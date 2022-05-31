---
sidebar_label: Job
title: core.Job
---

## Torque Objects

```python
class Torque(Job)
```

Template Torque object to contain Torque job initialization settings

**Attributes**:

- `jobname` _str_ - Job name
- `project` _str_ - Project name
- `nodes` _int_ - Number of nodes
- `ppn` _int_ - Number of processors (CPU)
- `walltime` _str_ - Job walltime (in format of HH:MM:SS)
- `gpus` _int_ - Number of processors (GPU)
- `job_fname` _str_ - Name of the PBS input file; default: `"job.pbs"`

### write\_job

```python
@build_dir
def write_job(output_dir: str = '.') -> None
```

Method to make the Torque job scheduler input file

**Arguments**:

- `output_dir` _str_ - Directory for the Torque input file; default: `"."`


**Returns**:

  None

## Slurm Objects

```python
class Slurm(Job)
```

Template Slurm object to contain Slurm job initialization settings

**Attributes**:

- `jobname` _str_ - Job name
- `nodes` _int_ - Number of nodes
- `ntasks_per_node` _int_ - Number of processors (CPU)
- `time` _str_ - Job time
- `gpus` _int_ - Number of processors (GPU)
- `job_fname` _str_ - Name of the Slurm input file; default: `"job.sh"`

### write\_job

```python
@build_dir
def write_job(output_dir: str = '.') -> None
```

Method to make the Slurm job scheduler input file

**Arguments**:

- `output_dir` _str_ - Directory for the Slurm input file; default: `.`


**Returns**:

  None
