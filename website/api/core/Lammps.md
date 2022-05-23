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

- `read_data_from` _System | str_ - System object that the data file will
  be read from. This can also be your data
  file name string if you do not generate your
  system via PMD; one of this attribute and
  `read_restart_from` has to be provided but not
  both (providing both will result in an error)
  ; default: `None`
- `read_restart_from` _Lammps | str_ - Lammps object that the last restart
  file created will be read from. This can
  also be your previous lammps file name
  string if you do not have previous Lammps
  object; one of this attribute and
  `read_data_from` has to be provided but
  not both (providing both will result in an
  error); default: `None`
- `force_foeld` _ForceField_ - Only needed if `read_data_from` or
  `read_restart_from` is provided as a file
  name string. This is needed for specifying
  potential form for LAMMPS input file
  ; default: `None`
- `atom_style` _str_ - LAMMPS
  [atom_style](https://docs.lammps.org/atom_style.html)
  to use during simulation; default: `"full"`
- `units` _str_ - LAMMPS [units](https://docs.lammps.org/units.html) to use
  during simulation; default: `"real"`
- `timestep` _float_ - LAMMPS
  [timestep](https://docs.lammps.org/timestep.html) to
  use during simulation; default: `1.0` (in unit of fs
  if `units` is `"real"`)
- `neighbor_skin` _float_ - LAMMPS
  [neighbor](https://docs.lammps.org/neighbor.html)
  skin size to use during the simulation; default:
  `2.0 Angstrom`
- `neighbor_every` _int_ - LAMMPS
  [neighbor](https://docs.lammps.org/neighbor.html)
  list checking frequency to use during the
  simulation; default: `1`
- `thermo` _int_ - LAMMPS [thermo](https://docs.lammps.org/thermo.html)
  to use during simulation; default: `1000 timestep`
- `lmp_input_fname` _str_ - Name of the LAMMPS input file; default:
  `"lmp.in"`

### add_procedure

```python
def add_procedure(procedure: Union[Procedure, List[Procedure]]) -> Lammps
```

Method to add simulation procedure

**Arguments**:

- `procedure` _Procedure_ - One of `Minimization`, `Equilibration`,
  `NPT`, `NVT`, `MSDMeasurement`, `TgMeasurement`, and `Deformation`

**Returns**:

- `Lammps` _Lammps_ - Lammps instance itself (builder design pattern)

### write_lammps

```python
@build_dir
def write_lammps(output_dir: str = '.') -> None
```

Method to make LAMMPS input files

**Arguments**:

- `output_dir` _str_ - Directory for the generated LAMMPS input file
  ; default: `"."`

**Returns**:

None
