Set the CPU governor to “performance” (once in your shell):
echo performance | sudo tee /sys/devices/system/cpu/cpu*/cpufreq/scaling_governor

sudo update-initramfs -u

# Method 1: Temporary fix
sudo sysctl -w kernel.perf_event_paranoid=0

# Method 2: Permanent fix  
echo 'kernel.perf_event_paranoid = 0' | sudo tee -a /etc/sysctl.conf
