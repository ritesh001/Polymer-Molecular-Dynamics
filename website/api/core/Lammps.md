---
sidebar_label: Lammps
title: core.Lammps
---

## Lammps Objects

```python
class Lammps()
```

Template object to contain LAMMPS initialization settings

**Attributes**:

- `data_fname` _str_ - File name of the data file, which will be read in by LAMMPS
  [read_data](https://docs.lammps.org/read_data.html) command
  
- `force_field` _str_ - Force field (`GAFF2` or `OPLS`)
  
- `atom_style` _str_ - LAMMPS [atom_style](https://docs.lammps.org/atom_style.html)
  to use during simulation; default=full
  
- `units` _str_ - LAMMPS [units](https://docs.lammps.org/units.html) to use during
  simulation; default=real
  
- `timestep` _float_ - LAMMPS [timestep](https://docs.lammps.org/timestep.html) to
  use during simulation; default=1 fs
  
- `neighbor_skin` _float_ - LAMMPS [neighbor](https://docs.lammps.org/neighbor.html)
  skin size to use during simulation; default=2.0 Angstrom
  
- `neighbor_every` _int_ - LAMMPS [neighbor](https://docs.lammps.org/neighbor.html)
  list checking frequency to use during simulation; default=1 fs
  
- `thermo` _int_ - LAMMPS [thermo](https://docs.lammps.org/thermo.html) to use during
  simulation; default=1000 timestep

### add\_procedure

```python
def add_procedure(procedure: Procedure)
```

Method to add simulation procedure

**Arguments**:

- `procedure` _Procedure_ - One of `minimization`, `equilibration`, or `Tg_measurement`
  

**Returns**:

- `Lammps` _Lammps_ - Lammps instance itself (builder design pattern)

### write\_input

```python
def write_input(output_dir: str, lmp_input_fname: str)
```

Method to make LAMMPS input files

**Arguments**:

- `output_dir` _str_ - Directory for all the generated LAMMPS input files
  

**Returns**:

  None

