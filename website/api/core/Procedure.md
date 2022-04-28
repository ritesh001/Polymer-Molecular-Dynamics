---
sidebar_label: Procedure
title: core.Procedure
---

## Minimization Objects

```python
class Minimization(Procedure)
```

Perform an energy minimization of the system, by iteratively adjusting atom coordinates.
Iterations are terminated when one of the stopping criteria is satisfied. At that point the
configuration will hopefully be in local potential energy minimum.

**Attributes**:

- `min_style` _str_ - Minimization algorithm, see [here](https://docs.lammps.org/min_style.html)
  for all options; default: `cg`
- `etol` _float_ - Stopping tolerance for energy (unitless); default: `10**(-6)`
- `ftol` _float_ - Stopping tolerance for force (force units); default: `10**(-8)`
- `maxiter` _int_ - Max iterations of minimizer; default: `10**5`
- `maxeval` _int_ - Max number of force/energy evaluations; default: `10**7`

## Equilibration Objects

```python
class Equilibration(Procedure)
```

Perform a 21-step amorphous polymer equilibration process.
Ref: Abbott, Hart, and Colina, Theoretical Chemistry Accounts, 132(3), 1-19, 2013.

**Attributes**:

- `Teq` _float_ - Target equilibration temperature; default: `300`
- `Peq` _float_ - Target equilibration pressure; default: `1`
- `Tmax` _float_ - Maximum temperature during the equilibration; default: `600`
- `Pmax` _float_ - Maximum pressure during the equilibration; default: `50000`
- `Tdamp` _str_ - Damping parameter for thermostats; default: `&#x27;$(100.0*dt)&#x27;`
- `Pdamp` _str_ - Damping parameter for barostats; default: `&#x27;$(100.0*dt)&#x27;`
- `dump_fname` _str_ - Name of the dump file; default: `&#x27;equil.lammpstrj&#x27;`
- `reset_timestep` _bool_ - Whether to reset timestep after the procedure; default: `True`

## TgMeasurement Objects

```python
class TgMeasurement(Procedure)
```

Perform glass transition temperature measurement of the system,
by iteratively cooling the system and equilibrate.

## NPT Objects

```python
class NPT(Procedure)
```

Perform the simulation under NPT ensemble.

## NVT Objects

```python
class NVT(Procedure)
```

Perform the simulation under NVT ensemble.
