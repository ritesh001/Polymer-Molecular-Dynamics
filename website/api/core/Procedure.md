---
sidebar_label: Procedure
title: core.Procedure
---

## Minimization Objects

```python
class Minimization(Procedure)
```

Perform an energy minimization of the system, by iteratively adjusting
atom coordinates. Iterations are terminated when one of the stopping
criteria is satisfied. At that point the configuration will hopefully be in
local potential energy minimum.

**Attributes**:

- `min_style` _str_ - Minimization algorithm, see
  [here](https://docs.lammps.org/min_style.html) for all
  options; default: `"cg"`

- `etol` _float_ - Stopping tolerance for energy (unitless); default:
  `10**(-6)`

- `ftol` _float_ - Stopping tolerance for force (force units); default:
  `10**(-8)`

- `maxiter` _int_ - Max iterations of minimizer; default: `10**5`

- `maxeval` _int_ - Max number of force/energy evaluations; default: `10**7`

## Equilibration Objects

```python
class Equilibration(Procedure)
```

Perform a 21-step amorphous polymer equilibration process. Ref: Abbott,
Hart, and Colina, Theoretical Chemistry Accounts, 132(3), 1-19, 2013.

**Attributes**:

- `Teq` _float_ - Target equilibration temperature; default: `300`

- `Peq` _float_ - Target equilibration pressure; default: `1`

- `Tmax` _float_ - Maximum temperature during the equilibration; default:
  `600`

- `Pmax` _float_ - Maximum pressure during the equilibration; default:
  `50000`

- `Tdamp` _str_ - Damping parameter for the thermostat; default:
  `"$(100.0*dt)"`

- `Pdamp` _str_ - Damping parameter for the barostat; default:
  `"$(100.0*dt)"`

- `dump_fname` _str_ - Name of the dump file; default: `"equil.lammpstrj"`

- `dump_every` _int_ - Dump every this many timesteps; default: `10000`

- `dump_image` _bool_ - Whether to dump a image file at the end of the run
  ; default: `False`

- `reset_timestep_before_run` _bool_ - Whether to reset timestep after the
  procedure; default: `True`

## NPT Objects

```python
class NPT(Procedure)
```

Perform the simulation under NPT ensemble (via Nose-Hoover thermostat
and barostat).

**Attributes**:

- `duration` _int_ - Duration of this NPT procedure (timestep unit)

- `Tinit` _float_ - Initial temperature

- `Tfinal` _float_ - Final temperature

- `Pinit` _float_ - Initial pressure

- `Pfinal` _float_ - Final pressure

- `Tdamp` _str_ - Damping parameter for the thermostat; default:
  `"$(100.0*dt)"`

- `Pdamp` _str_ - Damping parameter for the barostat; default:
  `"$(100.0*dt)"`

- `dump_fname` _str_ - Name of the dump file; default: `"npt.lammpstrj"`

- `dump_every` _int_ - Dump every this many timesteps; default: `10000`

- `dump_image` _bool_ - Whether to dump a image file at the end of the run
  ; default: `False`

- `reset_timestep_before_run` _bool_ - Whether to reset timestep after the
  procedure; default: `False`

## NVT Objects

```python
class NVT(Procedure)
```

Perform the simulation under NVT ensemble (via Nose-Hoover thermostat).

**Attributes**:

- `duration` _int_ - Duration of this NVT procedure (timestep unit)

- `Tinit` _float_ - Initial temperature

- `Tfinal` _float_ - Final temperature

- `Tdamp` _str_ - Damping parameter for thermostats; default:
  `"$(100.0*dt)"`

- `dump_fname` _str_ - Name of the dump file; default: `"nvt.lammpstrj"`

- `dump_every` _int_ - Dump every this many timesteps; default: `10000`

- `dump_image` _bool_ - Whether to dump a image file at the end of the run
  ; default: `False`

- `reset_timestep_before_run` _bool_ - Whether to reset timestep after the
  procedure; default: `False`

## MSDMeasurement Objects

```python
class MSDMeasurement(Procedure)
```

Perform the simulation under NVT ensemble (via Nose-Hoover thermostat).

**Attributes**:

- `duration` _int_ - Duration of this NVT procedure (timestep unit)

- `Tinit` _float_ - Initial temperature

- `Tfinal` _float_ - Final temperature

- `group` _str_ - The group of atoms that will be considered for MSD
  calculation. This has to be a string that matches the
  syntax of [group](https://docs.lammps.org/group.html)
  LAMMPS command (e.g. `"molecule <=50"`, `"type 1 2"`, etc)

- `create_block_every` _int_ - The time interval that new MSD calculation
  starting point will be created (e.g. for a
  1000 fs run, a `create_block_every` value of
  100fs would result in 10 blocks with 10
  different MSD starting point and length)
  ; default: `None`

- `result_folder_name` _str_ - The name of the folder that PMD creates and
  put result files in; default: `"result"`

- `Tdamp` _str_ - Damping parameter for thermostats; default:
  `"$(100.0*dt)"`

- `dump_fname` _str_ - Name of the dump file; default: `"nvt.lammpstrj"`

- `dump_every` _int_ - Dump every this many timesteps; default: `10000`

- `dump_image` _bool_ - Whether to dump a image file at the end of the run
  ; default: `False`

- `reset_timestep_before_run` _bool_ - Whether to reset timestep after the
  procedure; default: `False`

## TgMeasurement Objects

```python
class TgMeasurement(Procedure)
```

Perform glass transition temperature measurement of the system,
by iteratively cooling the system and equilibrate.

**Attributes**:

- `Tinit` _float_ - Initial temperature of the cooling process; default:
  `500`

- `Tfinal` _float_ - Final temperature of the cooling process; default:
  `100`

- `Tinterval` _float_ - Temperature interval of the cooling process
  ; default: `20`

- `step_duration` _int_ - Duration of each temperature step
  (timestep unit); default: `1000000`

- `pressure` _float_ - Pressure during the cooling process; default: `1`

- `Tdamp` _str_ - Damping parameter for the thermostat; default:
  `"$(100.0*dt)"`

- `Pdamp` _str_ - Damping parameter for the barostat; default:
  `"$(100.0*dt)"`

- `dump_fname` _str_ - Name of the dump file; default:
  `"Tg_measurement.lammpstrj"`

- `dump_every` _int_ - Dump every this many timesteps; default: `10000`

- `dump_image` _bool_ - Whether to dump a image file at the end of the run
  ; default: `False`

- `result_fname` _str_ - Name of the result file; default:
  `"temp_vs_density.txt"`

## Deformation Objects

```python
class Deformation(Procedure)
```

Perform a uniaxial tensile deformation in the x direction.
This can be used to calculate modulus and tensile strengths.

**Attributes**:

- `duration` _int_ - Duration of the deformation procedure (timestep unit)

- `erate` _float_ - Engineering strain rate. The units of the specified
  strain rate are 1/time

- `Tinit` _float_ - Initial temperature

- `Tfinal` _float_ - Final temperature

- `Tdamp` _str_ - Damping parameter for thermostats; default:
  `"$(100.0*dt)"`

- `print_every` _int_ - Print result to the result file every this many
  timesteps; default: `1000`

- `dump_fname` _str_ - Name of the dump file; default:
  `"deformation.lammpstrj"`

- `dump_every` _int_ - Dump every this many timesteps; default: `10000`

- `dump_image` _bool_ - Whether to dump a image file at the end of the run
  ; default: `False`

- `reset_timestep_before_run` _bool_ - Whether to reset timestep after the
  procedure; default: `False`

- `result_fname` _str_ - Name of the result file; default:
  `"stress_vs_strain.txt"`
