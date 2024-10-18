import os
import numpy as np
import matplotlib.pyplot as plt

def plot_peak_slip_rate_vs_time(base_dir, folders, label):
    # Initialize lists to hold the concatenated data
    times = []
    peak_slip_rates = []

    # Loop over each folder and read the data
    for folder in folders:
        file_path = os.path.join(base_dir, folder, 'global.dat')
        if os.path.exists(file_path):
            data = np.loadtxt(file_path)
            if times:
                # Adjust the time to be consecutive
                last_time = times[-1]
                data[:, 0] += last_time
            times.extend(data[:, 0])
            peak_slip_rates.extend(data[:, 1])
        else:
            print(f'File not found: {file_path}')

    # Convert lists to numpy arrays for plotting
    times = np.array(times)
    peak_slip_rates = np.array(peak_slip_rates)

    # Convert time from seconds to years
    times /= (60 * 60 * 24 * 365.25)

    # Plot the data
    plt.plot(times, peak_slip_rates, label=label)
    plt.yscale('log')
    plt.xlabel('Time (years)')
    plt.ylabel('Peak Slip Rate, m/s')
    plt.title('Peak Slip Rate vs Time')
    plt.grid(True)

# Example usage for three models
base_dir1 = '/Users/dliu/scratch/eqquasi.hbi/bp5.qdc.2000.dip50'
base_dir2 = '/Users/dliu/scratch/eqquasi.hbi/bp5.qdc.2000.dip70'
base_dir3 = '/Users/dliu/scratch/eqquasi.hbi/bp5.qdc.2000.dip90'
folders = [f'Q{i}' for i in range(10)]

plt.figure()
plot_peak_slip_rate_vs_time(base_dir1, folders, label='Dip 50')
plot_peak_slip_rate_vs_time(base_dir2, folders, label='70')
plot_peak_slip_rate_vs_time(base_dir3, folders, label='90')
plt.legend()
plt.show()
