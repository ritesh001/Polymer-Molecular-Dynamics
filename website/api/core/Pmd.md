---
sidebar_label: Pmd
title: core.Pmd
---

## Pmd Objects

```python
class Pmd()
```

Template object to perform tasks for Systems, Lammps, and Jobs
altogether (e.g. create data files, lammps input files, job scheduler
files, or config files)

**Attributes**:

- `system` _System_ - a System object

- `lammps` _Lammps or list[Lammps]_ - one or a list of Lammps objects

- `job` _Job or list[Job]_ - one or a list of Job objects

### create

```python
@build_dir
def create(output_dir: str = '.',
           save_config: bool = False,
           config_fname: str = 'config.yaml') -> None
```

Method to create files from all the pmd objects. This method can
can also automatically generate a config file if `save_config` input
argument is set to True.

**Arguments**:

- `output_dir` _str_ - Directory for all the generated files; default:
  `"."`

- `save_config` _bool_ - Whether to save a config file; default: `False`

- `config_fname` _str_ - Name of the config file; default:
  `"config.yaml"`


**Returns**:

  None

### save\_config

```python
@build_dir
def save_config(output_dir: str, config_fname: str = 'config.yaml')
```

Method to create a config file with all the details of the System,
Lammps, or Job settings. This method only creates the config file.

**Arguments**:

- `output_dir` _str_ - Directory for all the generated files; default:
  `"."`

- `config_fname` _str_ - Name of the config file; default:
  `"config.yaml"`


**Returns**:

  None

### load\_config

```python
@staticmethod
def load_config(config_file: str, output_dir: str = '.')
```

Method to load a config file and create all the objects listed in
the config file

**Arguments**:

- `config_file` _str_ - Config file to load

- `output_dir` _str_ - Directory for all the generated files; default:
  `"."`


**Returns**:

  None
