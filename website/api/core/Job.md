---
sidebar_label: Job
title: core.Job
---

## Job Objects

```python
class Job()
```

Template object to contain job initialization settings

**Attributes**:

- `jobname` _str_ - Job name
- `project` _str_ - Project name
- `nodes` _int_ - Number of nodes
- `ppn` _int_ - Number of processors (CPU)
- `walltime` _str_ - Job walltime
- `gpus` _int_ - Number of processors (GPU)

### write\_pbs

```python
def write_pbs(output_dir: str, pbs_fname: str) -> None
```

Method to make the PBS job scheduler input file

**Arguments**:

- `output_dir` _str_ - Directory for the PBS input file
- `pbs_fname` _str_ - Name of the PBS input file
  

**Returns**:

  None

### write\_slurm

```python
def write_slurm(output_dir: str, slurm_fname: str) -> None
```

Method to make the Slurm job scheduler input file

**Arguments**:

- `output_dir` _str_ - Directory for the Slurm input file
- `slurm_fname` _str_ - Name of the Slurm input file
  

**Returns**:

  None

