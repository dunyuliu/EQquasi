import os, argparse
import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

def plot_profile(caseDir, cycleFolder, fileName, cellSize, dip=90, var='effective_normal', profile='vertical', label=''):

    file_path = os.path.join(caseDir, cycleFolder, fileName)
    if os.path.exists(file_path):
        f_tdyna = os.path.join(caseDir, cycleFolder, 'tdyna.txt')
        tdyna = np.loadtxt(f_tdyna)
        global_data = np.loadtxt(os.path.join(caseDir, cycleFolder, 'global.dat'))
        index = np.where(global_data[:, 0] == tdyna[0])
        # Find the closest fault file
        snapshotRuptInit = [f for f in os.listdir(os.path.join(caseDir, cycleFolder)) if f.startswith('fault.') and f.endswith('.nc')]
        snapshotRuptInit = [f for f in snapshotRuptInit if f != 'fault.r.nc']
        time_indices = [int(f.split('.')[1]) for f in snapshotRuptInit]
        closest_index = min(time_indices, key=lambda x: abs(x - index[0][0]))
        closest_snapshot = f'fault.{closest_index:05d}.nc'
        print(f'Loading closest fault.*.nc snapshot file to nucleation is : {closest_snapshot}')
        fault_file_path = os.path.join(caseDir, cycleFolder, closest_snapshot)
        with nc.Dataset(fault_file_path, 'r') as fault_ds:
            varAtInit = fault_ds.variables[var][:]
            # You can now use fault_var_distribution as needed

        with nc.Dataset(file_path, 'r') as ds:
            var_distribution = ds.variables[var][:]
            if var=='slips':
                var_distribution = var_distribution - varAtInit

            if profile == 'horizontal':
                cellSize = cellSize / np.sin(dip / 180. * np.pi)
                profileDepth = 10e3  # meters
                print(cellSize)
                index = var_distribution.shape[0] - round(profileDepth / cellSize)
                print(index)
                profile_data = var_distribution[index, :]
            elif profile == 'vertical':
                profile_data = var_distribution[:, var_distribution.shape[1] // 2]
            elif profile == 'sum_horizontal':
                profiles = []
                for i in range(var_distribution.shape[0]):
                    profile_data = var_distribution[i, :]
                    if np.mean(profile_data) > 0.5:
                        profiles.append(profile_data)
                if profiles:
                    profile_data = np.mean(profiles, axis=0)
                else:
                    profile_data = np.zeros(var_distribution.shape[1])
            elif profile == 'sum_vertical':
                profiles = []
                for i in range(var_distribution.shape[1]):
                    profile_data = var_distribution[:, i]
                    if np.mean(profile_data) > 0.5:
                        profiles.append(profile_data)
                if profiles:
                    profile_data = np.mean(profiles, axis=0)
                else:
                    profile_data = np.zeros(var_distribution.shape[0])
            else:
                raise ValueError("Profile must be 'horizontal' or 'vertical' or 'sum_horizontal' or 'sum_vertical'")
            
    else:
        print(f'File not found: {file_path}')

    # Plot the data
    scale = 1
    unit = 'm'
    if var=='effective_normal' or var=='shear_strike' or var=='shear_dip':
        scale = 1e6
        unit = 'MPa'
    elif var=='slip_rate':
        unit = 'm/s'
    elif var=='state_variable':
        unit = ''
    plt.plot(profile_data/scale, label=label)
    plt.title('EQ Cycle ID is ' + cycleFolder)
    plt.xlabel('Grids along ' + profile + ' profile')
    plt.ylabel(var + f' ({unit})')   

# Parse command line arguments
parser = argparse.ArgumentParser(description='Generate stress profile plots.')
parser.add_argument('cycleFolder', type=str, help='Cycle folder containing the data files')
parser.add_argument('--var', type=str, default='slips', choices=['shear_strike', 'shear_dip', 'effective_normal', 'slips', 'slipd', 'slip_rate'], help='Variable to plot (default: slips)')
parser.add_argument('--file', type=str, default='fault.r.nc', help='Name of the file to read data from (default: fault.r.nc)')
parser.add_argument('--profile_type', type=str, default='horizontal', choices=['horizontal', 'vertical', 'sum_horizontal', 'sum_vertical'], help='Type of profile to plot (default: horizontal)')
args = parser.parse_args()
cycleFolder = args.cycleFolder
var = args.var
fileName = args.file
profile_type = args.profile_type
cellSize = 2000

plt.figure()
plot_profile('bp5.qdc.2000.dip50', cycleFolder, fileName, cellSize, dip=50, var=var, profile=profile_type, label='Dip50')
plot_profile('bp5.qdc.2000.dip70', cycleFolder, fileName, cellSize, dip=70, var=var, profile=profile_type, label='Dip70')
plot_profile('bp5.qdc.2000.dip90', cycleFolder, fileName, cellSize, dip=90, var=var, profile=profile_type, label='Dip90')
plt.legend()
plt.grid(True)
plt.show()