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

## Equilibration Objects

```python
class Equilibration(Procedure)
```

Perform a 21-step amorphous polymer equilibration process.
Ref: Abbott, Hart, and Colina, Theoretical Chemistry Accounts, 132(3), 1-19, 2013.

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

