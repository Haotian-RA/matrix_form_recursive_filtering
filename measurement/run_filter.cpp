// run_filter.cpp — End-to-end cascaded IIR filter benchmark
//
// Measures clock cycles per sample using wall-clock, TSC, and core cycle
// counters, then compares against PMCTest reference values.
//
// Compile-time constants (set via -D flags from run_filter.sh):
//   M_VAL  — number of biquads (filter_order / 2)
//   N_VAL  — block count per signal block group
//
// Runtime options:
//   --algo scalar|bf|ph|cr    (default: cr)
//   --arch haswell|skylake|meteorlake  (default: haswell)
//
// Example:
//   clang++ -std=c++20 -mavx2 -mfma -march=native -O2 -ffast-math
//           -DM_VAL=8 -DN_VAL=32 -lpthread run_filter.cpp -o run_filter
//   taskset -c 0 ./run_filter --algo cr --arch meteorlake

#include <chrono>
#include <cstdint>
#include <iostream>
#include <pthread.h>
#include <thread>
#include <sched.h>
#include <sys/syscall.h>
#include <unistd.h>
#include <linux/perf_event.h>
#include <sys/ioctl.h>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <fstream>
#include <string>
#include "../src/vcl/vectorclass.h"
#include "../include/filter.h"

// ─── Compile-time configuration ─────────────────────────────────────────────
#ifndef M_VAL
#define M_VAL 8
#endif

#ifndef N_VAL
#define N_VAL 8
#endif

#ifndef VECTOR_SIZE
#define VECTOR_SIZE 131072
#endif

#ifndef ITERS_VAL
#define ITERS_VAL 10000
#endif

using V = Vec8f;
using T = decltype(std::declval<V>().extract(0));
constexpr int L = V::size();
constexpr int N = N_VAL;
constexpr int M = M_VAL;
constexpr int vector_size = VECTOR_SIZE;
constexpr int ITERS = ITERS_VAL;
constexpr int WARMUP = 200;

// ─── Filter coefficients ────────────────────────────────────────────────────
constexpr T b1 = 2, b2 = 1, a1 = 1.3, a2 = -0.4;
constexpr T xi1 = 2, xi2 = 1, yi1 = -3, yi2 = -5;

// Generate M-element arrays at compile time
template<int K>
struct CoefArrays {
    T inits[K][4];
    T coefs[K][5];
    constexpr CoefArrays() : inits{}, coefs{} {
        for (int i = 0; i < K; ++i) {
            inits[i][0] = xi1; inits[i][1] = xi2;
            inits[i][2] = yi1; inits[i][3] = yi2;
            coefs[i][0] = 1;   coefs[i][1] = b1;
            coefs[i][2] = b2;  coefs[i][3] = a1;
            coefs[i][4] = a2;
        }
    }
};

static constexpr auto arrays = CoefArrays<M>();

// ─── TSC timing ─────────────────────────────────────────────────────────────
static inline uint64_t rdtsc_begin() {
    unsigned hi, lo;
    asm volatile("cpuid\n\trdtsc" : "=a"(lo), "=d"(hi) :: "rbx", "rcx");
    return (uint64_t(hi) << 32) | lo;
}

static inline uint64_t rdtsc_end() {
    unsigned hi, lo;
    asm volatile(
        "rdtscp\n\tmov %%edx,%0\n\tmov %%eax,%1\n\tcpuid"
        : "=r"(hi), "=r"(lo) :: "rax", "rbx", "rcx", "rdx");
    return (uint64_t(hi) << 32) | lo;
}

// ─── Core cycle counter (perf_event) ────────────────────────────────────────
class CoreCycleCounter {
    int fd;
    bool enabled;
public:
    CoreCycleCounter() : fd(-1), enabled(false) {
        struct perf_event_attr pe;
        memset(&pe, 0, sizeof(pe));
        pe.type = PERF_TYPE_HARDWARE;
        pe.config = PERF_COUNT_HW_CPU_CYCLES;
        pe.size = sizeof(pe);
        pe.disabled = 1;
        pe.exclude_kernel = 1;
        pe.exclude_hv = 1;
        pe.exclude_idle = 1;
        fd = syscall(__NR_perf_event_open, &pe, 0, -1, -1, 0);
        if (fd == -1) {
            std::cerr << "Warning: perf_event_open failed. "
                      << "Run: sudo sysctl -w kernel.perf_event_paranoid=0\n";
        } else {
            enabled = true;
        }
    }
    void start() {
        if (enabled) {
            ioctl(fd, PERF_EVENT_IOC_RESET, 0);
            ioctl(fd, PERF_EVENT_IOC_ENABLE, 0);
        }
    }
    void stop() {
        if (enabled) ioctl(fd, PERF_EVENT_IOC_DISABLE, 0);
    }
    uint64_t read() {
        if (!enabled) return 0;
        uint64_t count;
        ::read(fd, &count, sizeof(count));
        return count;
    }
    bool is_enabled() const { return enabled; }
    ~CoreCycleCounter() { if (fd != -1) close(fd); }
};

// ─── CPU frequency calibration ──────────────────────────────────────────────
double calibrate_tsc_frequency() {
    std::cout << "Calibrating TSC frequency... " << std::flush;
    double total = 0.0;
    for (int i = 0; i < 5; ++i) {
        auto t0 = std::chrono::high_resolution_clock::now();
        uint64_t c0 = rdtsc_begin();
        std::this_thread::sleep_for(std::chrono::milliseconds(200));
        uint64_t c1 = rdtsc_end();
        auto t1 = std::chrono::high_resolution_clock::now();
        double sec = std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count() / 1e9;
        total += (c1 - c0) / sec;
    }
    double freq = total / 5.0;
    std::cout << freq / 1e9 << " GHz\n";
    return freq;
}

// ─── PMCTest reference CPS data ─────────────────────────────────────────────
// Source: PMCTest measurements on three Intel micro-architectures (L=8, N=8)
//
// Layout: cps[order_idx] where order_idx: 0=order2, 1=order4, 2=order8, 3=order16
// These are PMCTest-predicted CPS (sum of per-stage measurements).

enum Arch { HASWELL = 0, SKYLAKE = 1, METEORLAKE = 2, NUM_ARCH = 3 };
enum Algo { SCALAR = 0, BF = 1, PH = 2, CR = 3, NUM_ALGO = 4 };

const char* arch_names[] = {"haswell", "skylake", "meteorlake"};
const char* algo_names[] = {"scalar", "bf", "ph", "cr"};

// Reference CPS: [algo][arch][order_idx]
// order_idx: 0=order2, 1=order4, 2=order8, 3=order16
const float ref_cps[NUM_ALGO][NUM_ARCH][4] = {
    // SCALAR
    {
        {27.614, 31.990, 47.428, 87.669},   // Haswell
        {23.268, 26.087, 41.483, 82.425},   // Skylake
        {17.279, 26.394, 36.262, 55.452},   // Meteor Lake
    },
    // BF (block filtering)
    {
        { 1.785,  5.580, 17.802, 49.280},   // Haswell
        { 1.418,  5.015, 13.330, 35.911},   // Skylake
        { 1.396,  2.890,  8.194, 23.865},   // Meteor Lake
    },
    // PH (N=8)
    {
        { 2.766,  4.447,  7.868, 14.509},   // Haswell
        { 2.479,  3.905,  6.808, 12.439},   // Skylake
        { 1.310,  2.433,  5.338, 11.008},   // Meteor Lake
    },
    // CR (N=8)
    {
        { 3.156,  5.320,  9.535, 18.128},   // Haswell
        { 2.862,  4.674,  8.424, 15.790},   // Skylake
        { 1.463,  2.793,  6.398, 13.526},   // Meteor Lake
    },
};

// Extended PH reference for different block sizes (Haswell only, N=8,16,32,64,128,256,512)
const float ref_ph_haswell_extended[4][7] = {
    { 2.766, 2.301, 2.214, 2.469, 2.713, 3.026, 3.752},  // order 2
    { 4.447, 3.666, 3.451, 3.990, 4.056, 4.354, 5.374},  // order 4
    { 7.868, 6.280, 5.890, 6.324, 6.577, 7.586, 9.064},  // order 8
    {14.509,11.505,10.631,11.273,11.721,15.154,16.179},   // order 16
};

// Extended CR reference for different block sizes (Haswell only, N=8,16,32,64,128,256,512)
const float ref_cr_haswell_extended[4][7] = {
    { 3.156, 2.324, 2.007, 2.101, 2.252, 2.501, 2.916},  // order 2
    { 5.320, 3.657, 3.011, 3.088, 3.270, 3.468, 4.389},  // order 4
    { 9.535, 6.407, 5.044, 5.067, 5.158, 5.426, 8.611},  // order 8
    {18.128,11.833, 9.065, 9.079,10.224,13.587,16.157},   // order 16
};

const int extended_blocksizes[] = {8, 16, 32, 64, 128, 256, 512};

// ─── Argument parsing ───────────────────────────────────────────────────────
struct Config {
    Algo algo = CR;
    Arch arch = HASWELL;
};

Config parse_args(int argc, char* argv[]) {
    Config cfg;
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--algo" && i + 1 < argc) {
            std::string val = argv[++i];
            if (val == "scalar") cfg.algo = SCALAR;
            else if (val == "bf")  cfg.algo = BF;
            else if (val == "ph")  cfg.algo = PH;
            else if (val == "cr")  cfg.algo = CR;
            else { std::cerr << "Unknown algo: " << val << "\n"; exit(1); }
        } else if (arg == "--arch" && i + 1 < argc) {
            std::string val = argv[++i];
            if (val == "haswell")       cfg.arch = HASWELL;
            else if (val == "skylake")  cfg.arch = SKYLAKE;
            else if (val == "meteorlake") cfg.arch = METEORLAKE;
            else { std::cerr << "Unknown arch: " << val << "\n"; exit(1); }
        } else if (arg == "--help" || arg == "-h") {
            std::cout << "Usage: " << argv[0]
                      << " [--algo scalar|bf|ph|cr] [--arch haswell|skylake|meteorlake]\n"
                      << "\nCompile-time constants (set via run_filter.sh):\n"
                      << "  --order 2|4|8|16      filter order (default: 16)\n"
                      << "  --blocksize N         block count (default: 8)\n";
            exit(0);
        }
    }
    return cfg;
}

// ─── Lookup reference CPS ───────────────────────────────────────────────────
float lookup_ref_cps(Config cfg) {
    int order_idx = -1;
    if      (M == 1) order_idx = 0;  // order 2
    else if (M == 2) order_idx = 1;  // order 4
    else if (M == 4) order_idx = 2;  // order 8
    else if (M == 8) order_idx = 3;  // order 16
    else return 0.0f;

    // For PH/CR on Haswell with extended block sizes, try exact match
    if (cfg.arch == HASWELL && (cfg.algo == PH || cfg.algo == CR)) {
        const auto& table = (cfg.algo == PH)
            ? ref_ph_haswell_extended : ref_cr_haswell_extended;
        for (int i = 0; i < 7; ++i) {
            if (extended_blocksizes[i] == N)
                return table[order_idx][i];
        }
    }

    // Fall back to N=8 reference
    return ref_cps[cfg.algo][cfg.arch][order_idx];
}

// ─── Main ───────────────────────────────────────────────────────────────────
int main(int argc, char* argv[]) {
    Config cfg = parse_args(argc, argv);

    std::cout << "========================================\n"
              << "  IIR Filter Benchmark\n"
              << "========================================\n\n";

    // Configuration
    std::cout << "Configuration:\n"
              << "  Algorithm:    " << algo_names[cfg.algo] << "\n"
              << "  Architecture: " << arch_names[cfg.arch] << "\n"
              << "  Filter order: " << 2 * M << " (" << M << " biquads)\n"
              << "  Block count:  " << N << " (NL = " << N * L << " samples/group)\n"
              << "  Vector size:  " << vector_size << " samples\n"
              << "  Iterations:   " << ITERS << "\n\n";

    // TSC calibration
    double TSC_HZ = calibrate_tsc_frequency();

    // Core cycle counter
    CoreCycleCounter core;
    if (core.is_enabled())
        std::cout << "Core cycle counter: enabled\n\n";
    else
        std::cout << "Core cycle counter: disabled (TSC only)\n\n";

    // Pin to core 0
    cpu_set_t cpus;
    CPU_ZERO(&cpus);
    CPU_SET(0, &cpus);
    pthread_setaffinity_np(pthread_self(), sizeof(cpus), &cpus);
    std::this_thread::sleep_for(std::chrono::milliseconds(100));

    // Build filter and data
    alignas(64) std::array<T, vector_size> in{}, out;
    in[0] = 1;
    auto _F = Filter<V, M, N>(arrays.coefs, arrays.inits);

    // Warm up
    for (int i = 0; i < WARMUP; ++i)
        _F(in.begin(), in.end(), out.begin());

    // Measure loop overhead
    uint64_t oh_tsc = 0, oh_core = 0;
    double oh_ns = 0;
    {
        auto t0 = std::chrono::high_resolution_clock::now();
        uint64_t c0 = rdtsc_begin();
        core.start();
        for (int i = 0; i < ITERS; ++i) { asm volatile(""); }
        core.stop();
        uint64_t c1 = rdtsc_end();
        auto t1 = std::chrono::high_resolution_clock::now();
        oh_ns = double(std::chrono::duration_cast<std::chrono::nanoseconds>(t1 - t0).count()) / ITERS;
        oh_tsc = (c1 - c0) / ITERS;
        if (core.is_enabled()) oh_core = core.read() / ITERS;
    }

    // ─── Benchmark ──────────────────────────────────────────────────────
    auto wc0 = std::chrono::high_resolution_clock::now();
    uint64_t tsc0 = rdtsc_begin();
    core.start();

    for (int i = 0; i < ITERS; ++i)
        _F(in.begin(), in.end(), out.begin());

    core.stop();
    uint64_t tsc1 = rdtsc_end();
    auto wc1 = std::chrono::high_resolution_clock::now();

    // ─── Results ────────────────────────────────────────────────────────
    double total_ns = double(std::chrono::duration_cast<std::chrono::nanoseconds>(wc1 - wc0).count());
    double avg_ns = total_ns / ITERS - oh_ns;

    double avg_tsc = double(tsc1 - tsc0) / ITERS - oh_tsc;
    double tsc_cps = avg_tsc / vector_size;

    double core_cps = 0;
    if (core.is_enabled()) {
        double avg_core = double(core.read()) / ITERS - oh_core;
        core_cps = avg_core / vector_size;
    }

    float predicted = lookup_ref_cps(cfg);

    std::cout << "========================================\n"
              << "  Results\n"
              << "========================================\n\n";

    std::cout << "Wall-clock:  " << avg_ns << " ns/iter\n";
    std::cout << "TSC:         " << tsc_cps << " cycles/sample\n";

    if (core.is_enabled()) {
        std::cout << "Core cycles: " << core_cps << " cycles/sample\n";
        std::cout << "Boost ratio: " << core_cps / tsc_cps
                  << "x (core/TSC, >1 = turbo active)\n";
    }

    std::cout << "\n";

    if (predicted > 0) {
        double best_cps = core.is_enabled() ? core_cps : tsc_cps;
        const char* method = core.is_enabled() ? "Core" : "TSC";
        double diff = best_cps - predicted;
        double pct = (diff / predicted) * 100.0;

        std::cout << "PMCTest prediction (" << algo_names[cfg.algo]
                  << ", " << arch_names[cfg.arch]
                  << ", N=" << N << "): " << predicted << " cycles/sample\n";
        std::cout << method << " measured: " << best_cps << " cycles/sample\n";
        std::cout << "Difference: " << std::abs(diff) << " ("
                  << (diff >= 0 ? "+" : "") << pct << "%)\n";
    } else {
        std::cout << "No PMCTest reference available for this configuration.\n";
        std::cout << "(N=" << N << " reference data may not be in the table.)\n";
    }

    std::cout << "\n========================================\n";

    return 0;
}
