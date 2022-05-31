from typing import Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import optimize

from pmd.util.Log import Pmdlogging


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
                                     ;default: `None`

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
    box_com = np.zeros([frames, 3], np.float)

    # preallocate msd vectors
    msd_dict = {}
    for type_id in set(id2type):
        msd_dict[type_id] = np.zeros(frames, np.float)

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
