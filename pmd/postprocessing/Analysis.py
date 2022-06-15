import math
import os
import warnings
from typing import List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize
from scipy.optimize import curve_fit
from sklearn.linear_model import LinearRegression

from pmd.util.Log import Pmdlogging

warnings.filterwarnings("ignore")


def calculate_Tg(result_fname: str,
                 make_plot: bool = True,
                 append_result_to_yaml: Optional[str] = None) -> int:
    '''Method to calculate glass transition temperature based on the
    result file obtained from TgMeasurement Procedure

    Parameters:
        result_fname (str): Name of the result file from TgMeasurement
                            Procedure
        make_plot (bool): Whether to make a plot to visualize the fitting
        append_result_to_yaml (str): YAML file name to append result value to
                                     ; default: `None`

    Returns:
        Tg (float): Glass transition temperature of the system
    '''

    # TODO: implement append_result_to_yaml with dict structure like
    # {result: {Tg: value}} or {result_{smiles}: {Tg: value}}

    def piecewise_linear(x, x0, y0, k1, k2):
        return np.piecewise(
            x, [x < x0],
            [lambda x: k1 * x + y0 - k1 * x0, lambda x: k2 * x + y0 - k2 * x0])

    # Read result file (should have first column as temperature,
    # and second column as the corresponding density)
    df = pd.read_csv(
        result_fname,
        header=1,
        usecols=[1, 2],  # TODO: relax hard-coded col index
        names=['Temp', 'Rho'],
        delim_whitespace=True)
    x = np.array(df['Temp'])
    y = np.array(df['Rho'])

    p0 = (250, 0.85, -0.0001, 0.0001)
    p, e = optimize.curve_fit(piecewise_linear, x, y, p0)
    xd = np.linspace(x.min(), x.max(), 100)

    if make_plot:
        plt.rcParams.update({'font.size': 12})
        plt.plot(x, y, "ko", fillstyle='none', markersize=8)
        plt.plot(xd, piecewise_linear(xd, *p), 'r')
        plt.ylabel('Density $(g/cm^3)$')
        plt.xlabel('Temperature ($K$)')
        plt.savefig('temp_vs_density.png', dpi=300)

    Tg = p[0]
    Pmdlogging.info(f'Glass transition temperature: {Tg}')

    return Tg


def calculate_diffusivity(result_folder: str = 'result',
                          block_list: List[int] = [1, 2, 4, 5, 10, 20],
                          time_array: List[int] = [
                              5000000, 10000000, 20000000, 50000000, 100000000,
                              200000000
                          ]):
    '''Method to calculate diffusivity based on the files in the
    result folder obtained from MSDMeasurement Procedure

    Parameters:
        result_folder (str): Name of the result folder from MSDMeasurement
                             Procedure
        block_list (list): A list of number of blocks to use
        time_array (list): A list of time as the start and end points of
                        fitting region

    Returns:
        D (float): Diffusivity of the system
    '''

    def read_files(dir: str, block_list: List[int]):
        process_dict = {}
        block_dict = {}
        file_counter = 0
        for root, dirs, files in os.walk(dir):
            for file in files:
                parts = file.split('.')[0].split('_')
                if file.endswith(".txt") and len(parts) == 3:
                    file_counter += 1
                    # read file and add to process_dict with starting times
                    # as keys and dfs as values
                    df = pd.read_csv(os.path.join(dir, file),
                                     header=1,
                                     usecols=[0, 1],
                                     delim_whitespace=True,
                                     names=['time', 'msd'])

                    process_dict[int(parts[1])] = df

        duration = int(parts[2])
        # validate all blocks
        for block in block_list:
            if file_counter % block != 0:
                Pmdlogging.warning(f'block numbers in the block list should '
                                   f'be a factor of {file_counter}')
                continue

            for b in range(block):
                start_time = int(duration / block * b)
                final_time = int(duration / block * (b + 1))
                df_block = process_dict[start_time].copy()
                df_block.columns = ['time', f'block{b}']
                final_index = df_block.index[df_block['time'] ==
                                             final_time].to_list()[0]

                df_block_sliced = df_block.loc[:final_index]
                df_block_sliced['time'] = df_block_sliced['time'] - \
                    df_block_sliced['time'][0]

                if b == 0:
                    df_final = df_block_sliced
                else:
                    df_final = df_final.merge(df_block_sliced, on='time')

            block_dict[block] = df_final
        return block_dict

    def find_nearest(array, value):
        array = np.asarray(array)
        nearest_index = (np.abs(array - value)).argmin()
        if array[nearest_index] > value:
            nearest_index -= 1
        return nearest_index

    def main_calculation(block_dict, t1, t2):
        dif_mean = np.zeros(len(block_dict))
        msd_slope_log_mean = np.zeros(len(block_dict))
        dif_std = np.zeros(len(block_dict))
        msd_slope_log_std = np.zeros(len(block_dict))
        block_list = np.zeros(len(block_dict))

        for index, num_blocks in enumerate(sorted(block_dict)):
            time = list(block_dict[num_blocks]['time'])

            cut_point = find_nearest(time, t1)
            cut_point2 = find_nearest(time, t2)

            if cut_point == cut_point2:
                cut_point = cut_point2 - 5

            dif = np.zeros(num_blocks)
            msd_slope_log = np.zeros(num_blocks)
            msd_slope_linear = np.zeros(num_blocks)

            for n in range(num_blocks):
                x = np.log10(
                    block_dict[num_blocks]['time'].loc[cut_point:cut_point2])
                y = np.log10(block_dict[num_blocks]
                             [f'block{n}'].loc[cut_point:cut_point2])
                x_sklearn = block_dict[num_blocks][[
                    'time'
                ]].loc[cut_point:cut_point2]
                y_sklearn = block_dict[num_blocks][f'block{n}'].loc[
                    cut_point:cut_point2]

                def fit_func(x, a, b):
                    return a * x + b

                msd_slope_log[n] = curve_fit(fit_func, x, y)[0][0]

                msd_slope_linear[n] = LinearRegression().fit(
                    x_sklearn, y_sklearn).coef_

                # 0.1 is for A^2/fs to cm^2/s and 6 is for 3 dimensions
                dif[n] = np.log10(msd_slope_linear[n] / 6 * 0.1)

            dif_mean[index] = dif.mean()
            msd_slope_log_mean[index] = msd_slope_log.mean()
            dif_std[index] = dif.std()
            msd_slope_log_std[index] = msd_slope_log.std()
            block_list[index] = num_blocks

        best_index = np.argmin(np.abs(msd_slope_log_mean - 1))
        return dif_mean[best_index], dif_std[best_index], msd_slope_log_mean[
            best_index], block_list[best_index]

    curr_best_dif_mean = 0
    curr_best_dif_std = 0
    curr_best_slope = 0
    curr_best_block = 0
    curr_best_t1 = 0
    curr_best_t2 = 0
    block_dict = read_files(result_folder, block_list)
    for t1 in range(len(time_array) - 1):
        for t2 in range(t1 + 1, len(time_array)):
            dif_mean, dif_std, slope, block = main_calculation(
                block_dict, time_array[t1], time_array[t2])
            if not math.isnan(dif_mean):
                if np.abs(slope - 1) < np.abs(curr_best_slope - 1):
                    curr_best_dif_mean = dif_mean
                    curr_best_dif_std = dif_std
                    curr_best_slope = slope
                    curr_best_block = block
                    curr_best_t1 = time_array[t1]
                    curr_best_t2 = time_array[t2]

    Pmdlogging.info(f'************** Best Result: '
                    f'{curr_best_t1}-{curr_best_t2}ns **************')
    Pmdlogging.info(curr_best_dif_mean, curr_best_dif_std, curr_best_slope,
                    curr_best_block)

    return (curr_best_dif_mean, curr_best_dif_std, curr_best_slope,
            curr_best_block)


# if __name__ == '__main__':
#     main()


def calculate_MSD(r, ir, box_bounds, id2type=[]):
    '''Method to calculate mean squared displacement for each type as given in
    `id2type`; NOTE: does NOT do any sort of block averaging; assumes mass = 1
    for all beads; does not account for changes in box size

    Parameters:
        r: unscaled (but wrapped) coordinates (format as read in from
           `read_lammpstrj`)

        ir: image flags (format as read in from `read_lammpstrj`)

        box_bounds: boundaries of the box (format as read in from
                    `read_lammpstrj`)

        id2type: array that maps atom id to type (format as read in from
                 `read_lammpstrj`)

    Returns:
        msd_dict: dict of the calculated MSDs for each type
    '''

    # set up some constants
    frames = len(r)
    box_size = np.array([
        box_bounds[0][0][1] - box_bounds[0][0][0],
        box_bounds[0][1][1] - box_bounds[0][1][0],
        box_bounds[0][2][1] - box_bounds[0][2][0]
    ])

    # allocate an array for the box center of mass which needs to be
    # subtracted off
    box_com = np.zeros([frames, 3], float)

    # preallocate msd vectors
    msd_dict = {}
    for type_id in set(id2type):
        msd_dict[type_id] = np.zeros(frames, float)

    # loop over frames
    for t in range(frames):
        # calculate the center of mass of the entire box
        for atom in range(1, len(r[0])):
            box_com[t] += r[t][atom] + ir[t][atom] * box_size
        box_com[t] = box_com[t] / (len(r[0]) - 1)

        # loop over atoms
        for atom in range(1, len(id2type)):
            # calculate how much the bead has moved reletive to the center of
            # mass (note that this is a vector equation)
            diff = (r[t][atom] + ir[t][atom] * box_size - box_com[t]) - (
                r[0][atom] + ir[0][atom] * box_size - box_com[0])
            # mean squared displacement is this difference dotted with itself
            msd_dict[id2type[atom]][t] += diff.dot(diff)

    # scale MSD by the number of beads of each type, to get the average MSD
    for type_id in set(id2type):
        msd_dict[type_id] = msd_dict[type_id] / sum(id2type == type_id)
    # this is needed since id2type has a dummy entry of 0 at index 0 so that
    # it is indexed by LAMMPS atom_id
    del msd_dict[0]

    return msd_dict
