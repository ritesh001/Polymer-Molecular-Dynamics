---
sidebar_label: Builder
title: core.Builder
---

## EMC Objects

```python
class EMC(Builder)
```

Object to perform system structure generation using
[EMC](http://montecarlo.sourceforge.net/): Enhanced Monte Carlo package.
This object should be used as input argument of `System` or `Lammps`
objects

**Attributes**:

- `force_field` _str_ - Force field, options are `"pcff"`, `"opls-aa"`,
  `"opls-ua"`, and `"trappe"`

## PSP Objects

```python
class PSP(Builder)
```

Object to perform system structure generation using
[PSP](https://github.com/Ramprasad-Group/PSP): Polymer Structure Predictor
package. This object should be used as input argument of `System` or
`Lammps` objects

**Attributes**:

- `force_field` _str_ - Force field, options are `"opls-lbcc"`,
  `"opls-cm1a"`, `"gaff2-gasteiger"`, and `"gaff2-am1bcc"`

