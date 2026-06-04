#!/bin/bash

# Ultra-optimized Multi-run PMC benchmark script
# Minimal setup, compile once, execute many times
# Now includes double-filtering: 1st within runs, 2nd across collection
# Usage: ./multi_run_benchmark_ultra.sh [num_runs] [vector_type] [N] [function_name] [vector_size] [M] [second_sigma]

# Auto-detect parameter format: support both 
# ./script runs vector_type N function_name (old format)
# ./script runs vector_type function_name (new format - skip N)

NUM_RUNS=${1:-1000}
VECTOR_TYPE=${2:-"Vec8f"}

# Auto-detect if $3 is N (number) or FUNCTION_NAME (string)
if [[ "$3" =~ ^[0-9]+$ ]]; then
    # $3 is a number, use old format: runs vector_type N function_name
    N_VALUE=$3
    FUNCTION_NAME=${4:-"permute"}
    VECTOR_SIZE=${5:-128}  # Only used for filter functions
    M_VALUE=${6:-1}        # Filter order, only used for filter functions
    SECOND_SIGMA=${7:-2.0} # Second-level filtering threshold (default: 2.0σ)
else
    # $3 is not a number, use new format: runs vector_type function_name
    N_VALUE=8  # Use default
    FUNCTION_NAME=${3:-"permute"}
    VECTOR_SIZE=${4:-128}  # Only used for filter functions
    M_VALUE=${5:-1}        # Filter order, only used for filter functions
    SECOND_SIGMA=${6:-2.0} # Second-level filtering threshold (default: 2.0σ)
fi

echo "=== Ultra-Optimized Multi-Run PMC Benchmark ==="
echo "Configuration:"
echo "  Number of runs: $NUM_RUNS"
echo "  Vector Type: $VECTOR_TYPE"
echo "  N (block size): $N_VALUE"
echo "  Function: $FUNCTION_NAME"
if [[ "$FUNCTION_NAME" == *"Filter"* || "$FUNCTION_NAME" == *"_F"* ]]; then
    echo "  Vector Size: $VECTOR_SIZE"
    echo "  Filter Order (M): $M_VALUE"
fi
echo "  Each run: 100 repetitions → select repetition with median Core cyc"
echo "  Second-level filtering: ${SECOND_SIGMA}σ (filter Core cyc outliers from collection)"
if [[ "$3" =~ ^[0-9]+$ ]]; then
    echo "  Parameter format: Long (runs vector_type N function_name)"
else
    echo "  Parameter format: Short (runs vector_type function_name, N=default)"
fi
echo

# Create directory for storing individual run results
RESULTS_DIR="$HOME/matrix_form_recursive_filtering/measurement/multi_run_results"
mkdir -p "$RESULTS_DIR"

# Clear previous results
rm -f "$RESULTS_DIR"/*.txt
rm -f "$HOME/matrix_form_recursive_filtering/measurement/filtered_means.txt"

# Path to PMCTest directory and files
PMCTEST_DIR="$HOME/testp/PMCTest"
PMCTEST_FILE="$PMCTEST_DIR/PMCTestB.cpp"
EXECUTABLE="$PMCTEST_DIR/a.out"

echo "=== SETUP PHASE ==="

# Determine function category and calling pattern
case $FUNCTION_NAME in
    "permute")
        FUNCTION_TYPE="simple"
        TEST_CALL="M_out = permuteV(M_in);"
        SETUP_TYPE="vector"
        ;;
    "depermute")
        FUNCTION_TYPE="simple"
        TEST_CALL="M_out = depermuteV(M_in);"
        SETUP_TYPE="vector"
        ;;
    "F"|"FirCoreOrderTwo")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = F(M_in);"
        OBJECT_DECL="FirCoreOrderTwo<V,N> F(b1,b2,xi1,xi2);"
        SETUP_TYPE="vector"
        ;;
    "I"|"IirCoreOrderTwo")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = I(M_in);"
        OBJECT_DECL="IirCoreOrderTwo<V,N> I(b1,b2,a1,a2,xi1,xi2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "IM"|"IirCoreOrderTwo")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = I(M_in);"
        OBJECT_DECL="IirCoreOrderTwo<V,N> I(b1,b2,a1,a2,xi1,xi2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "IV"|"IirCoreOrderTwo")
        FUNCTION_TYPE="object"
        TEST_CALL="V_out = I(V_in);"
        OBJECT_DECL="IirCoreOrderTwo<V,N> I(b1,b2,a1,a2,xi1,xi2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "IS"|"IirCoreOrderTwo")
        FUNCTION_TYPE="object"
        TEST_CALL="S_out = I(S_in);"
        OBJECT_DECL="IirCoreOrderTwo<V,N> I(b1,b2,a1,a2,xi1,xi2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "BF"|"BlockFiltering")
        FUNCTION_TYPE="object"
        TEST_CALL="V_out = BF(V_in);"
        OBJECT_DECL="BlockFiltering<V> BF(b1,b2,a1,a2,xi1,xi2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "PS"|"PartSolutionV")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = PS(M_in);"
        OBJECT_DECL="PartSolutionV<V,N> PS(a1,a2);"
        SETUP_TYPE="vector"
        ;;
    "HS"|"HomoSolutionV")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = HS(M_in);"
        OBJECT_DECL="HomoSolutionV<V,N> HS(a1,a2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "PH"|"HomoSolutionV")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = HS(PS(M_in));"
        OBJECT_DECL="HomoSolutionV<V,N> HS(a1,a2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "CRF"|"CyclicReduction")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = CR.forward<0>(M_in);"
        OBJECT_DECL="CyclicReduction<V,N> CR(a1,a2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "CRB"|"CyclicReduction")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = CR.backward(M_in);"
        OBJECT_DECL="CyclicReduction<V,N> CR(a1,a2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "CR"|"CyclicReduction")
        FUNCTION_TYPE="object"
        TEST_CALL="M_out = CR(M_in);"
        OBJECT_DECL="CyclicReduction<V,N> CR(a1,a2,yi1,yi2);"
        SETUP_TYPE="vector"
        ;;
    "Filter"|"_F")
        FUNCTION_TYPE="filter"
        TEST_CALL="auto r = _F(in.begin(),in.end(),out.begin());"
        OBJECT_DECL="auto _F = Filter<V,M,N>(coefs,inits);"
        SETUP_TYPE="filter"
        ;;
    *)
        FUNCTION_TYPE="simple"
        TEST_CALL="M_out = ${FUNCTION_NAME}V(M_in);"
        SETUP_TYPE="vector"
        ;;
esac

echo "Function type: $FUNCTION_TYPE ($SETUP_TYPE setup)"

# Generate coefficient arrays for filter functions
generate_filter_coefficients() {
    local M=$1
    if [ $M -eq 1 ]; then
        INITS_ARRAY="{xi1,xi2,yi1,yi2}"
        COEFS_ARRAY="{1,b1,b2,a1,a2}"
    else
        local inits_inner="{xi1,xi2,yi1,yi2}"
        local coefs_inner="{1,b1,b2,a1,a2}"
        INITS_ARRAY="{"
        COEFS_ARRAY="{"
        for ((i=1; i<=M; i++)); do
            INITS_ARRAY+="$inits_inner"
            COEFS_ARRAY+="$coefs_inner"
            if [ $i -lt $M ]; then
                INITS_ARRAY+=","
                COEFS_ARRAY+=","
            fi
        done
        INITS_ARRAY+="}"
        COEFS_ARRAY+="}"
    fi
}

echo "Configuring PMCTestB.cpp..."

# Create backup if it doesn't exist
if [ -f "$PMCTEST_FILE" ] && [ ! -f "${PMCTEST_FILE}.backup" ]; then
    cp "$PMCTEST_FILE" "${PMCTEST_FILE}.backup"
    echo "  Created backup of original PMCTestB.cpp"
fi

# Quick configuration using sed (minimal changes)
echo "  Updating vector type and block size..."
sed -i "s/using V = Vec[0-9]*[fd];/using V = $VECTOR_TYPE;/" "$PMCTEST_FILE"
sed -i "s/constexpr static int N = [0-9]*;/constexpr static int N = $N_VALUE;/" "$PMCTEST_FILE"

if [ "$SETUP_TYPE" = "filter" ]; then
    echo "  Setting up filter configuration..."
    sed -i "s/constexpr static int M = [0-9]*;/constexpr static int M = $M_VALUE;/" "$PMCTEST_FILE"
    sed -i "s/constexpr static int vector_size = [0-9]*;/constexpr static int vector_size = $VECTOR_SIZE;/" "$PMCTEST_FILE"
    
    # Generate and insert coefficients
    generate_filter_coefficients $M_VALUE
    
    # Remove old coefficients and add new ones
    sed -i '/^\/\/ filter parameter/d' "$PMCTEST_FILE"
    sed -i '/^constexpr T inits\[M\]\[4\] = /d' "$PMCTEST_FILE"
    sed -i '/^constexpr T coefs\[M\]\[5\] = /d' "$PMCTEST_FILE"
    
    # Insert new coefficients after M definition (FIXED: removed extra blank line)
    sed -i "/^constexpr static int M = $M_VALUE;/a\\
// filter parameter\\
constexpr T inits[M][4] = $INITS_ARRAY;\\
constexpr T coefs[M][5] = $COEFS_ARRAY;" "$PMCTEST_FILE"
fi

echo "  Updating test call to: $TEST_CALL"

# FIXED: Update the test call to match ALL possible output variable patterns
sed -i "/Test code start/,/Test code end/{
    /\(M_out\|V_out\|S_out\) = /c\\        $TEST_CALL
    /auto r = /c\\        $TEST_CALL
}" "$PMCTEST_FILE"

echo "Configuration complete!"

echo
echo "=== COMPILATION PHASE ==="
echo "Compiling PMCTest (one-time compilation)..."

# Change to PMCTest directory
cd "$PMCTEST_DIR"

# Compile A file if needed
if [ PMCTestA.cpp -nt a64.o ] ; then
    echo "  Compiling PMCTestA.cpp..."
    g++ -O2 -c -m64 -oa64.o PMCTestA.cpp
    if [ $? -ne 0 ] ; then 
        echo "Error: Failed to compile PMCTestA.cpp"
        exit 1
    fi
fi

# Compile B file and link
echo "  Compiling and linking PMCTestB.cpp..."
clang++ -std=c++20 -mavx2 -mfma -O2 -ffast-math a64.o -march=native -lpthread PMCTestB.cpp

if [ $? -ne 0 ] ; then 
    echo "Error: Failed to compile PMCTestB.cpp"
    exit 1
fi

echo "Compilation successful!"
echo

echo "=== EXECUTION PHASE ==="
echo "Starting $NUM_RUNS benchmark runs (execution only, no recompilation)..."

# Progress tracking
progress_interval=$((NUM_RUNS / 20))
if [ $progress_interval -eq 0 ]; then
    progress_interval=1
fi

# Record start time
start_time=$(date +%s)
successful_runs=0
failed_runs=0

echo "Progress: "

# Run the benchmark multiple times (execution only)
for ((i=1; i<=NUM_RUNS; i++)); do
    # Show progress with timing estimates
    if [ $((i % progress_interval)) -eq 0 ] || [ $i -eq $NUM_RUNS ]; then
        current_time=$(date +%s)
        elapsed=$((current_time - start_time))
        if [ $i -gt 0 ]; then
            avg_time_per_run=$((elapsed * 1000 / i))  # milliseconds
            estimated_total=$((elapsed * NUM_RUNS / i))
            remaining=$((estimated_total - elapsed))
            remaining_min=$((remaining / 60))
            remaining_sec=$((remaining % 60))
        else
            avg_time_per_run=0
            remaining_min=0
            remaining_sec=0
        fi
        
        progress=$((i * 100 / NUM_RUNS))
        printf "\rProgress: [%-20s] %d%% (%d/%d) - %dms/run, ETA: %dm%ds    " \
               "$(printf '%*s' $((progress/5)) '' | tr ' ' '#')" \
               $progress $i $NUM_RUNS \
               $((avg_time_per_run)) \
               $remaining_min $remaining_sec
    fi
    
    # Execute the benchmark directly to the run file (taskset for CPU affinity, silent execution)
    taskset -c 0 "$EXECUTABLE" > "$RESULTS_DIR/run_${i}.txt" 2>/dev/null
    
    # Check if benchmark succeeded
    if [ $? -eq 0 ] && [ -s "$RESULTS_DIR/run_${i}.txt" ]; then
        ((successful_runs++))
    else
        ((failed_runs++))
        # Create empty file to maintain numbering
        touch "$RESULTS_DIR/run_${i}.txt"
    fi
done

# Record end time
end_time=$(date +%s)
total_time=$((end_time - start_time))

echo
echo
echo "Execution completed!"
echo "  Total time: ${total_time}s ($((total_time/60))m $((total_time%60))s)"
echo "  Average per run: $((total_time*1000/NUM_RUNS))ms"
echo "  Successful runs: $successful_runs"
echo "  Failed runs: $failed_runs"
echo

if [ $successful_runs -eq 0 ]; then
    echo "Error: No successful runs to analyze!"
    exit 1
fi

echo "=== ANALYSIS PHASE ==="
echo "Analyzing results..."

# Run the multi-run analysis
cd "$HOME/matrix_form_recursive_filtering/measurement"
python3 data_analyze.py "$RESULTS_DIR" "$NUM_RUNS" "$SECOND_SIGMA"

if [ $? -eq 0 ]; then
    echo
    echo "=== Analysis Complete ==="
    echo "Results saved to:"
    echo "  - Individual runs: $RESULTS_DIR/"
    echo "  - Filtered means: ~/matrix_form_recursive_filtering/measurement/filtered_means.txt"
else
    echo "Error: Multi-run analysis failed"
    exit 1
fi

# Usage examples
echo
echo "Usage examples:"
echo "  # Unit functions (both formats supported):"
echo "  ./pmc_measure.sh 1000 Vec8f BF            # Block Filtering (N=8 default)"
echo "  ./pmc_measure.sh 1000 Vec8f 16 BF         # Block Filtering (N=16 explicit)"
echo "  ./pmc_measure.sh 100 Vec8f F              # FIR Core (N=8 default)"
echo "  ./pmc_measure.sh 100 Vec8f 8 F            # FIR Core (N=8 explicit)"
echo "  ./pmc_measure.sh 200 Vec8f I              # IIR Core"  
echo "  ./pmc_measure.sh 500 Vec8f PS             # Particular Solution"
echo "  ./pmc_measure.sh 1000 Vec8f HS            # Homogeneous Solution"
echo "  ./pmc_measure.sh 1000 Vec8f PH            # PH Decomposition"
echo "  ./pmc_measure.sh 1000 Vec8f CRF           # Cyclic Reduction forward"
echo "  ./pmc_measure.sh 1000 Vec8f CRB           # Cyclic Reduction backward"
echo "  ./pmc_measure.sh 1000 Vec8f CR            # Cyclic Reduction"
echo ""
echo "  # Filter functions:"
echo "  ./pmc_measure.sh 1000 Vec8f Filter 128 1      # input size 128, M=1 (single filter), N=8 default"
echo "  ./pmc_measure.sh 1000 Vec8f 8 Filter 128 2    # input size 128, M=2 (two filters), N=8 explicit"
echo "  ./pmc_measure.sh 1000 Vec8f 16 Filter 256 4   # input size 256, M=4 (four filters), N=16"
