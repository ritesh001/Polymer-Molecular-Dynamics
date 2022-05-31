import sys

import numpy as np


def _read_header(f):
    f.readline()  # ITEM: TIMESTEP
    timestep = int(f.readline())

    f.readline()  # ITEM: NUMBER OF ATOMS
    num_atoms = int(f.readline())

    f.readline()  # ITEM: BOX BOUNDS xx yy zz
    line = f.readline()
    line = line.split()
    xlo = float(line[0])
    xhi = float(line[1])
    line = f.readline()
    line = line.split()
    ylo = float(line[0])
    yhi = float(line[1])
    line = f.readline()
    line = line.split()
    zlo = float(line[0])
    zhi = float(line[1])

    return timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi


def read_lammpstrj(fname,
                   num_frames=float('inf'),
                   skip_beginning=0,
                   skip_between=0):
    '''Method to read in a lammps trajectory
    NOTE: assumes that the number of atoms in the simulation is fixed;
    NOTE: also assumes that the coordinates are wrapped, so x or xs type
          coordinates are allowed but not xu or xsu

    Parameters:
        fname (str): filename string or 'stdin' (or a value that evaluates
                     to false) for reading from standard in

        num_frames (int): optional number of frames to read before stopping,
                          defaults to reading in all frames

        skip_beginning (int): skip this many frames at the beginning of the
                              dump file

        skip_between (int): skip this many frames between saved frames

    Returns:
        r: num_frames by num_atoms+1 by 3 array of wrapped and unscaled
                      coordinates (indexed by frame number then atom id)

        ir: num_frames by num_atoms+1 by 3 array of image flags

        timestep: num_frames length array of timesteps

        box_bounds: 3D array to store boundaries of the box, indexed by frame,
                    x/y/z, then lower/upper

        id2type: num_atoms+1 length arrays to map atom id to type id
                 (if available, may be None)

        id2mol: num_atoms+1 length arrays to map atom id to molecule id
                  (if available, may be None)

        mol2ids: num_mols+1 length list of atom id arrays corresponding to the
                 molecules (if available, may be None)
    '''

    # allow reading from standard input
    if not fname or fname == 'stdin':
        f = sys.stdin
    else:
        f = open(fname, 'r')

    # read in the initial header
    frame = 0
    init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = _read_header(f)

    # skip the beginning frames, if requested
    for skippedframe in range(skip_beginning):
        f.readline()  # ITEM: ATOMS
        # loop over the atoms lines
        for atom in range(num_atoms):
            f.readline()
        init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = _read_header(
            f)

    # preallocate arrays, if possible
    if num_frames < float('inf'):
        alloc = num_frames
        inf_frames = False
    else:
        alloc = 1
        inf_frames = True
    # 1D array of timesteps
    timestep = np.zeros(alloc, np.int)
    # 3D array to store boundaries of the box, indexed by frame, x/y/z, then
    # lower/upper
    box_bounds = np.zeros([alloc, 3, 2], np.float)

    timestep[frame] = init_timestep
    box_bounds[frame][0][0] = xlo
    box_bounds[frame][0][1] = xhi
    box_bounds[frame][1][0] = ylo
    box_bounds[frame][1][1] = yhi
    box_bounds[frame][2][0] = zlo
    box_bounds[frame][2][1] = zhi

    # NOTE: using num_atoms+1 here so that the arrays are indexed by their
    # LAMMPS atom id
    # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    r = np.zeros([alloc, num_atoms + 1, 3], np.float)
    # 3D array of x, y, z image flags, r[frame][id][coordinate]
    ir = np.zeros([alloc, num_atoms + 1, 3], np.int)

    # array to map from atom id to molecule id, builds this from the first
    # frame, if available
    id2mol = np.zeros(num_atoms + 1, np.int)
    # array to map from atom id to type, builds this from the first frame,
    # if available
    id2type = np.zeros(num_atoms + 1, np.int)

    # separately do the first ATOMS section so that we can initialize things,
    # build the id2mol and id2type arrays, and so that the main loop starts
    # with reading in the header
    line = f.readline()
    line = line.split()
    id_index = line.index("id") - 2
    if "mol" in line:
        mol_index = line.index("mol") - 2
    else:
        mol_index = None
    if "type" in line:
        type_index = line.index("type") - 2
    else:
        type_index = None

    if "x" in line:
        scaled = False
        x_index = line.index("x") - 2
        y_index = line.index("y") - 2
        z_index = line.index("z") - 2
    elif "xs" in line:
        scaled = True
        x_index = line.index("xs") - 2
        y_index = line.index("ys") - 2
        z_index = line.index("zs") - 2
    else:
        print("ERROR: x coordinate not found in lammps trajectory",
              file=sys.stderr)
        return

    if "ix" in line:
        ix_index = line.index("ix") - 2
        iy_index = line.index("iy") - 2
        iz_index = line.index("iz") - 2
    else:
        print("ERROR: x image flag not found in lammps trajectory",
              file=sys.stderr)
        return

    # loop over the atoms lines for the first frame separately, the rest of the
    # frames will be read in below
    for atom in range(num_atoms):
        line = f.readline()
        line = line.split()

        # get the atom id
        my_id = int(line[id_index])

        # x, y, z coordinates
        r[frame][my_id][0] = float(line[x_index])
        r[frame][my_id][1] = float(line[y_index])
        r[frame][my_id][2] = float(line[z_index])

        # unscale, if necessary
        if scaled:
            r[frame][my_id][0] = r[frame][my_id][0] * (
                box_bounds[frame][0][1] -
                box_bounds[frame][0][0]) + box_bounds[frame][0][0]
            r[frame][my_id][1] = r[frame][my_id][1] * (
                box_bounds[frame][1][1] -
                box_bounds[frame][1][0]) + box_bounds[frame][1][0]
            r[frame][my_id][2] = r[frame][my_id][2] * (
                box_bounds[frame][2][1] -
                box_bounds[frame][2][0]) + box_bounds[frame][2][0]

        # x, y, z image flags
        ir[frame][my_id][0] = int(line[ix_index])
        ir[frame][my_id][1] = int(line[iy_index])
        ir[frame][my_id][2] = int(line[iz_index])

        # if available, buidl the i2mol and id2type arrays
        if mol_index is not None:
            id2mol[my_id] = int(line[mol_index])
        if type_index is not None:
            id2type[my_id] = int(line[type_index])

    # build the reverse of the id2mol array
    # this is a 2D array with rows of (potentially) varying length, so nest a
    # numpy array into a python list
    if mol_index is not None:
        num_mols = id2mol.max()
        mol2ids = [[]]
        for molid in range(1, num_mols + 1):
            mol2ids.append(np.where(id2mol == molid)[0])
    else:
        num_mols = None
        mol2ids = None

    # loop over number of num_frames frames, if num_frames is infinite, will
    # loop over all the frames in the file
    # this is the frame counter for frames actually read in
    frame = 1
    # this is the actual frame count in the file (not counting the ones skipped
    # in the beginning
    frame_attempt = 0
    while frame < num_frames:

        frame_attempt += 1

        # try to read in a new header
        try:
            (my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo,
             my_zhi) = _read_header(f)
        except Exception:
            print("WARNING: hit end of file when reading in {} at frame {}".
                  format(fname, skip_beginning + frame_attempt),
                  file=sys.stderr)
            break

        # skip the frame if between frames to be read in and restart the loop
        if frame_attempt % (skip_between + 1) > 0:
            f.readline()  # ITEM: ATOMS
            # loop over the atoms lines
            for atom in range(num_atoms):
                f.readline()
            continue

        # if we don't know how many frames to read in, have to allocate more
        # memeory for the arrays
        if inf_frames:
            timestep = np.append(timestep, 0)

            box_bounds = np.concatenate(
                (box_bounds, np.zeros([1, 3, 2], np.float)))

            r = np.concatenate((r, np.zeros([1, num_atoms + 1, 3], np.float)))
            ir = np.concatenate((ir, np.zeros([1, num_atoms + 1, 3],
                                              np.float)))

        # update the timestep and box size arrays
        timestep[frame] = my_timestep
        box_bounds[frame][0][0] = my_xlo
        box_bounds[frame][0][1] = my_xhi
        box_bounds[frame][1][0] = my_ylo
        box_bounds[frame][1][1] = my_yhi
        box_bounds[frame][2][0] = my_zlo
        box_bounds[frame][2][1] = my_zhi

        f.readline()  # ITEM: ATOMS
        # loop over the atoms lines
        for atom in range(num_atoms):
            line = f.readline()
            line = line.split()

            # get the atom id
            my_id = int(line[id_index])

            # x, y, z coordinates
            r[frame][my_id][0] = float(line[x_index])
            r[frame][my_id][1] = float(line[y_index])
            r[frame][my_id][2] = float(line[z_index])

            # unscale, if necessary
            if scaled:
                r[frame][my_id][0] = r[frame][my_id][0] * (
                    box_bounds[frame][0][1] -
                    box_bounds[frame][0][0]) + box_bounds[frame][0][0]
                r[frame][my_id][1] = r[frame][my_id][1] * (
                    box_bounds[frame][1][1] -
                    box_bounds[frame][1][0]) + box_bounds[frame][1][0]
                r[frame][my_id][2] = r[frame][my_id][2] * (
                    box_bounds[frame][2][1] -
                    box_bounds[frame][2][0]) + box_bounds[frame][2][0]

            # x, y, z image flags
            ir[frame][my_id][0] = int(line[ix_index])
            ir[frame][my_id][1] = int(line[iy_index])
            ir[frame][my_id][2] = int(line[iz_index])

        frame += 1

    return r, ir, timestep, box_bounds, id2type, id2mol, mol2ids


def read_lammpstrj_by_type(fname,
                           types,
                           num_frames=float('inf'),
                           skip_beginning=0,
                           skip_between=0):
    '''Method to read in a lammps trajectory, but only extract data of certain
    bead types;

    Parameters:
        fname (str): filename string or 'stdin' (or a value that evaluates
                     to false) for reading from standard in

        types (int[]): types of atom to be read in; defaults to all types
                       of beads

        num_frames (int): optional number of frames to read before stopping,
                          defaults to reading in all frames

        skip_beginning (int): skip this many frames at the beginning of the
                              dump file

        skip_between (int): skip this many frames between saved frames

    Returns:
        r: num_frames by num_atoms+1 by 3 array of wrapped and unscaled
                      coordinates (indexed by frame number then atom id)

        ir: num_frames by num_atoms+1 by 3 array of image flags

        timestep: num_frames length array of timesteps

        box_bounds: 3D array to store boundaries of the box, indexed by frame,
                    x/y/z, then lower/upper

        id2type: num_atoms+1 length arrays to map atom id to type id
                 (if available, may be None)

        id2index: num_atoms+1 length arrays to map atom id to index id
                  (if available, may be None)
    '''

    # allow reading from standard input
    if not fname or fname == 'stdin':
        f = sys.stdin
    else:
        f = open(fname, 'r')

    # read in the initial header
    frame = 0
    init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = _read_header(f)

    # skip the beginning frames, if requested
    for skippedframe in range(skip_beginning):
        print(f"Skipping {skippedframe}")
        f.readline()  # ITEM: ATOMS
        # loop over the atoms lines
        for atom in range(num_atoms):
            f.readline()
        init_timestep, num_atoms, xlo, xhi, ylo, yhi, zlo, zhi = _read_header(
            f)

    # preallocate arrays, if possible
    if num_frames < float('inf'):
        alloc = num_frames
        inf_frames = False
    else:
        alloc = 1
        inf_frames = True

    # 1D array of timesteps
    timestep = np.zeros(alloc, np.int)
    # 3D array to store boundaries of the box, indexed by frame, x/y/z,
    # then lower/upper
    box_bounds = np.zeros([alloc, 3, 2], np.float)

    timestep[frame] = init_timestep
    box_bounds[frame][0][0] = xlo
    box_bounds[frame][0][1] = xhi
    box_bounds[frame][1][0] = ylo
    box_bounds[frame][1][1] = yhi
    box_bounds[frame][2][0] = zlo
    box_bounds[frame][2][1] = zhi

    # separately do the first ATOMS section so that we can initialize things,
    # build the id2mol and id2type arrays, and so that the main loop starts
    # with reading in the header
    line = f.readline()
    line = line.split()
    id_index = line.index("id") - 2
    # if "mol" in line:
    #     mol_index = line.index("mol") - 2
    # else:
    #     mol_index = None
    if "type" in line:
        type_index = line.index("type") - 2
    else:
        type_index = None

    if "x" in line:
        scaled = False
        x_index = line.index("x") - 2
        y_index = line.index("y") - 2
        z_index = line.index("z") - 2
    elif "xs" in line:
        scaled = True
        x_index = line.index("xs") - 2
        y_index = line.index("ys") - 2
        z_index = line.index("zs") - 2
    else:
        print("ERROR: x coordinate not found in lammps trajectory",
              file=sys.stderr)

    if "ix" in line:
        ix_index = line.index("ix") - 2
        iy_index = line.index("iy") - 2
        iz_index = line.index("iz") - 2
    else:
        print("ERROR: x image flag not found in lammps trajectory",
              file=sys.stderr)
        return

    # NOTE: using num_atoms+1 here so that the arrays are indexed by their
    # LAMMPS atom id
    # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    r_temp = np.zeros([num_atoms + 1, 3], np.float)
    # 3D array of x, y, z image flags, r[frame][id][coordinate]
    ir_temp = np.zeros([num_atoms + 1, 3], np.int)

    # array to map from atom id to type, builds this from the first frame,
    # if available
    index2type = np.zeros(num_atoms + 1, np.int)
    # array to map from atom id to index, builds this from the first frame,
    # if available
    index2id = np.zeros(num_atoms + 1, np.int)

    num_types = 0
    # loop over the atoms lines for the first frame separately, the rest of the
    # frames will be read in below
    for atom in range(num_atoms):
        line = f.readline()
        line = line.split()

        # get the atom id
        my_id = int(line[id_index])

        # build the index2type array
        my_type = int(line[type_index])
        index2type[my_id] = my_type
        if my_type in types:
            num_types += 1

        # x, y, z coordinates
        r_temp[my_id][0] = float(line[x_index])
        r_temp[my_id][1] = float(line[y_index])
        r_temp[my_id][2] = float(line[z_index])

        # unscale, if necessary
        if scaled:
            r_temp[my_id][0] = r_temp[my_id][0] * (
                box_bounds[frame][0][1] -
                box_bounds[frame][0][0]) + box_bounds[frame][0][0]
            r_temp[my_id][1] = r_temp[my_id][1] * (
                box_bounds[frame][1][1] -
                box_bounds[frame][1][0]) + box_bounds[frame][1][0]
            r_temp[my_id][2] = r_temp[my_id][2] * (
                box_bounds[frame][2][1] -
                box_bounds[frame][2][0]) + box_bounds[frame][2][0]

        # x, y, z image flags
        ir_temp[my_id][0] = int(line[ix_index])
        ir_temp[my_id][1] = int(line[iy_index])
        ir_temp[my_id][2] = int(line[iz_index])

    # NOTE: using num_types+1 here so that the arrays are indexed by their
    # LAMMPS atom id
    # 3D array of x, y, z coordinates, r[frame][id][coordinate]
    r = np.zeros([alloc, num_types + 1, 3], np.float)
    # 3D array of x, y, z image flags, r[frame][id][coordinate]
    ir = np.zeros([alloc, num_types + 1, 3], np.int)

    # array to map from atom id to type, builds this from the first frame,
    # if available
    id2type = np.zeros(num_types + 1, np.int)
    # array to map from atom id to index, builds this from the first frame,
    # if available
    id2index = np.zeros(num_types + 1, np.int)

    # store the temporary data into real arrays
    my_id = 0
    for atom in range(num_atoms):
        index = atom + 1
        if index2type[index] in types:
            my_id += 1
            # x, y, z coordinates
            r[frame][my_id][0] = r_temp[index][0]
            r[frame][my_id][1] = r_temp[index][1]
            r[frame][my_id][2] = r_temp[index][2]

            # x, y, z image flags
            ir[frame][my_id][0] = ir_temp[index][0]
            ir[frame][my_id][1] = ir_temp[index][1]
            ir[frame][my_id][2] = ir_temp[index][2]

            id2type[my_id] = index2type[index]
            id2index[my_id] = index
            index2id[index] = my_id

    # loop over number of num_frames frames, if num_frames is infinite,
    # will loop over all the frames in the file
    # this is the frame counter for frames actually read in
    frame = 1
    # this is the actual frame count in the file (not counting the ones
    # skipped in the beginning)
    frame_attempt = 0
    while frame < num_frames:

        frame_attempt += 1

        # try to read in a new header
        try:
            (my_timestep, my_num_atoms, my_xlo, my_xhi, my_ylo, my_yhi, my_zlo,
             my_zhi) = _read_header(f)
        except Exception:
            print("WARNING: hit end of file when reading in {} at frame {}".
                  format(fname, skip_beginning + frame_attempt),
                  file=sys.stderr)
            break

        # skip the frame if between frames to be read in and restart the loop
        if frame_attempt % (skip_between + 1) > 0:
            f.readline()  # ITEM: ATOMS
            # loop over the atoms lines
            for atom in range(num_atoms):
                f.readline()
            continue

        # if we don't know how many frames to read in, have to allocate more
        # memeory for the arrays
        if inf_frames:
            timestep = np.append(timestep, 0)

            box_bounds = np.concatenate(
                (box_bounds, np.zeros([1, 3, 2], np.float)))

            r = np.concatenate((r, np.zeros([1, num_types + 1, 3], np.float)))
            ir = np.concatenate((ir, np.zeros([1, num_types + 1, 3],
                                              np.float)))

        # update the timestep and box size arrays
        timestep[frame] = my_timestep
        box_bounds[frame][0][0] = my_xlo
        box_bounds[frame][0][1] = my_xhi
        box_bounds[frame][1][0] = my_ylo
        box_bounds[frame][1][1] = my_yhi
        box_bounds[frame][2][0] = my_zlo
        box_bounds[frame][2][1] = my_zhi

        f.readline()  # ITEM: ATOMS
        # loop over the atoms lines
        for atom in range(num_atoms):
            line = f.readline()
            line = line.split()

            # get the atom id
            index = int(line[id_index])
            if index2type[index] in types:
                my_id = index2id[index]

                # x, y, z coordinates
                r[frame][my_id][0] = float(line[x_index])
                r[frame][my_id][1] = float(line[y_index])
                r[frame][my_id][2] = float(line[z_index])

                # unscale, if necessary
                if scaled:
                    r[frame][my_id][0] = r[frame][my_id][0] * (
                        box_bounds[frame][0][1] -
                        box_bounds[frame][0][0]) + box_bounds[frame][0][0]
                    r[frame][my_id][1] = r[frame][my_id][1] * (
                        box_bounds[frame][1][1] -
                        box_bounds[frame][1][0]) + box_bounds[frame][1][0]
                    r[frame][my_id][2] = r[frame][my_id][2] * (
                        box_bounds[frame][2][1] -
                        box_bounds[frame][2][0]) + box_bounds[frame][2][0]

                # x, y, z image flags
                ir[frame][my_id][0] = int(line[ix_index])
                ir[frame][my_id][1] = int(line[iy_index])
                ir[frame][my_id][2] = int(line[iz_index])

        print(f"Reading frame {frame}")
        frame += 1

    print('===============Summary===============')
    print('Total number of beads =', num_atoms)
    print('Total number of selected beads =', num_types)

    return r, ir, timestep, box_bounds, id2type, id2index
