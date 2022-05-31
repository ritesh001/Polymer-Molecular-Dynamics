---
sidebar_label: Analysis
title: postprocessing.Analysis
---

### calculate\_Tg

```python
def calculate_Tg(result_fname: str,
                 make_plot: bool = True,
                 append_result_to_yaml: Optional[str] = None) -> int
```

Method to calculate glass transition temperature based on the
result file obtained from TgMeasurement Procedure

**Arguments**:

- `result_fname` _str_ - Name of the result file from TgMeasurement
  Procedure
- `make_plot` _bool_ - Whether to make a plot to visualize the fitting
- `append_result_to_yaml` _str_ - YAML file name to append result value to
- `;default` - `None`


**Returns**:

- `Tg` _float_ - Glass transition temperature of the system

### calculate\_MSD

```python
def calculate_MSD(r, ir, box_bounds, id2type=[])
```

Method to calculate mean squared displacement for each type as given in
`id2type`; NOTE: does NOT do any sort of block averaging; assumes mass = 1
for all beads; does not account for changes in box size

**Arguments**:

- `r` - unscaled (but wrapped) coordinates (format as read in from
  `read_lammpstrj`)

- `ir` - image flags (format as read in from `read_lammpstrj`)

- `box_bounds` - boundaries of the box (format as read in from
  `read_lammpstrj`)

- `id2type` - array that maps atom id to type (format as read in from
  `read_lammpstrj`)


**Returns**:

- `msd_dict` - dict of the calculated MSDs for each type
