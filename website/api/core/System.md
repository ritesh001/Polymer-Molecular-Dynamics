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
  
- `data_fname` _str_ - File name of the output data file, which will be read in by
  LAMMPS [read_data](https://docs.lammps.org/read_data.html)
  command

### update\_smiles

```python
def update_smiles(smiles: str) -> System
```

Method to update the SMILES of this system

**Arguments**:

- `smiles` _str_ - Updated SMILES string
  

**Returns**:

- `system` _System_ - System instance itself (builder design pattern)

### update\_natoms\_total

```python
def update_natoms_total(natoms_total: int) -> System
```

Method to update total number of atoms of the system

**Arguments**:

- `natoms_total` _int_ - Updated total number of atoms
  

**Returns**:

- `system` _System_ - System instance itself (builder design pattern)

### update\_force\_field

```python
def update_force_field(force_field: Str) -> System
```

Method to update force field of the system

**Arguments**:

- `force_field` _str_ - Updated force field
  

**Returns**:

- `system` _System_ - System instance itself (builder design pattern)

### write\_data

```python
def write_data(output_dir: str, cleanup: bool = True) -> None
```

Method to make LAMMPS data file (which contains coordinates and force
field parameters)

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS data file
  
- `cleanup` _bool_ - Whether to clean up files other than the LAMMPS data file PSP
  generated
  

**Returns**:

  None

