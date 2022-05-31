---
sidebar_label: TrajectoryReader
title: postprocessing.TrajectoryReader
---

### read\_lammpstrj

```python
def read_lammpstrj(fname,
                   num_frames=float('inf'),
                   skip_beginning=0,
                   skip_between=0)
```

Method to read in a lammps trajectory
NOTE: assumes that the number of atoms in the simulation is fixed;
NOTE: also assumes that the coordinates are wrapped, so x or xs type
coordinates are allowed but not xu or xsu

**Arguments**:

- `fname` _str_ - filename string or &#x27;stdin&#x27; (or a value that evaluates
  to false) for reading from standard in

- `num_frames` _int_ - optional number of frames to read before stopping,
  defaults to reading in all frames

- `skip_beginning` _int_ - skip this many frames at the beginning of the
  dump file

- `skip_between` _int_ - skip this many frames between saved frames


**Returns**:

- `r` - num_frames by num_atoms+1 by 3 array of wrapped and unscaled
  coordinates (indexed by frame number then atom id)

- `ir` - num_frames by num_atoms+1 by 3 array of image flags

- `timestep` - num_frames length array of timesteps

- `box_bounds` - 3D array to store boundaries of the box, indexed by frame,
  x/y/z, then lower/upper

- `id2type` - num_atoms+1 length arrays to map atom id to type id
  (if available, may be None)

- `id2mol` - num_atoms+1 length arrays to map atom id to molecule id
  (if available, may be None)

- `mol2ids` - num_mols+1 length list of atom id arrays corresponding to the
  molecules (if available, may be None)

### read\_lammpstrj\_by\_type

```python
def read_lammpstrj_by_type(fname,
                           types,
                           num_frames=float('inf'),
                           skip_beginning=0,
                           skip_between=0)
```

Method to read in a lammps trajectory, but only extract data of certain
bead types;

**Arguments**:

- `fname` _str_ - filename string or &#x27;stdin&#x27; (or a value that evaluates
  to false) for reading from standard in

- `types` _int[]_ - types of atom to be read in; defaults to all types
  of beads

- `num_frames` _int_ - optional number of frames to read before stopping,
  defaults to reading in all frames

- `skip_beginning` _int_ - skip this many frames at the beginning of the
  dump file

- `skip_between` _int_ - skip this many frames between saved frames


**Returns**:

- `r` - num_frames by num_atoms+1 by 3 array of wrapped and unscaled
  coordinates (indexed by frame number then atom id)

- `ir` - num_frames by num_atoms+1 by 3 array of image flags

- `timestep` - num_frames length array of timesteps

- `box_bounds` - 3D array to store boundaries of the box, indexed by frame,
  x/y/z, then lower/upper

- `id2type` - num_atoms+1 length arrays to map atom id to type id
  (if available, may be None)

- `id2index` - num_atoms+1 length arrays to map atom id to index id
  (if available, may be None)
