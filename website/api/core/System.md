---
sidebar_label: System
title: core.System
---

## System Objects

```python
class System()
```

Template object to contain System initialization settings

**Attributes**:

- `smiles` _str_ - SMILES string of the polymer (use * as connecting point)
  
- `density` _float_ - Density of the system
  
- `force_field` _str_ - Force field (One of `gaff2` and `opls`)
  
- `natoms_total` _int_ - Total number of atoms of the system
  
- `natoms_per_chain` _int_ - Number of atoms of the polymer, one of this
  attribute and `mw_per_chain` has to be
  provided but not both (providing both will
  result in an error); default: `None`
  
- `mw_per_chain` _int_ - Molecular weight of the polymer, one of this
  attribute and `natoms_per_chain` has to be
  provided but not both (providing both will
  result in an error); default: `None`
  
- `data_fname(str)` - File name of the data file, which will be read in by LAMMPS
  [read_data](https://docs.lammps.org/read_data.html) command,
  only enter it at instantialization if you will not use the
  `write_data` function

### write\_data

```python
def write_data(output_dir: str, data_fname: str = 'amor_data.lmps') -> None
```

Method to make LAMMPS data file (which contains coordinates and force
field parameters)

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS data file
  
- `data_fname` _str_ - File name of the data file, which will be read in by LAMMPS
  [read_data](https://docs.lammps.org/read_data.html) command
  

**Returns**:

  None

