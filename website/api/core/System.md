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

- `force_field` _ForceField_ - Force field (One of `GAFF2` and `OPLS`)

- `natoms_total` _int_ - Total number of atoms of the system

- `natoms_per_chain` _int_ - Number of atoms per polymer chain, one of this
  attribute, `mw_per_chain`, and `ru_per_chain`
  has to be provided but not both (providing both
  will result in an error); default: `None`

- `mw_per_chain` _int_ - Molecular weight of the polymer, one of this
  attribute, `natoms_per_chain`, and `ru_per_chain`
  has to be provided but not both (providing both
  will result in an error); default: `None`

- `ru_per_chain` _int_ - Number of repeating unit per polymer chain, one of
  this attribute, `natoms_per_chain`, and
  `mw_per_chain` has to be provided but not both
  (providing both will result in an error); default:
  `None`

- `data_fname` _str_ - File name of the output data file, which will be read in by
  LAMMPS [read_data](https://docs.lammps.org/read_data.html)
  command; default: `"data.lmps"`

### write\_data

```python
def write_data(output_dir: str = '.', cleanup: bool = True) -> None
```

Method to make LAMMPS data file (which contains coordinates and force
field parameters)

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS data file
  ; default: `"."`

- `cleanup` _bool_ - Whether to clean up files other than the LAMMPS data
  file PSP generated


**Returns**:

  None

## SolventSystem Objects

```python
class SolventSystem(System)
```

Template object to contain System with solvents initialization settings

**Attributes**:

- `smiles` _str_ - SMILES string of the polymer (use * as connecting point)

- `solvent_smiles` _str_ - SMILES string of the solvent

- `ru_nsolvent_ratio` _float_ - The ratio of total number of repeating units
  in the system and total number of solvent
  molecules

- `density` _float_ - Density of the system

- `force_field` _Force Field_ - Force field (One of `GAFF2` and `OPLS`)

- `natoms_total` _int_ - Total number of atoms of the system

- `natoms_per_chain` _int_ - Number of atoms per polymer chain, one of this
  attribute, `mw_per_chain`, and `ru_per_chain`
  has to be provided but not both (providing both
  will result in an error); default: `None`

- `mw_per_chain` _int_ - Molecular weight of the polymer, one of this
  attribute, `natoms_per_chain`, and `ru_per_chain`
  has to be provided but not both (providing both
  will result in an error); default: `None`

- `ru_per_chain` _int_ - Number of repeating unit per polymer chain, one of
  this attribute, `natoms_per_chain`, and
  `mw_per_chain` has to be provided but not both
  (providing both will result in an error); default:
  `None`

- `data_fname` _str_ - File name of the output data file, which will be read in by
  LAMMPS [read_data](https://docs.lammps.org/read_data.html)
  command; default: `"data.lmps"`

### write\_data

```python
def write_data(output_dir: str = '.', cleanup: bool = True) -> None
```

Method to make LAMMPS data file (which contains coordinates and force
field parameters)

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS data file
  ; default: `"."`

- `cleanup` _bool_ - Whether to clean up files other than the LAMMPS data
  file PSP generated


**Returns**:

  None
