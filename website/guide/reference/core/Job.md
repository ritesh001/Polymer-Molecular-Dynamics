---
sidebar_label: Job
title: core.Job
---

## Job Objects

```python
class Job()
```

pmd.core.Job.Job
Template object to contain job initialization settings

**Attributes**:

- `jobname` - str
  Job name
- `project` - str
  Project name
- `nodes` - int
  Number of nodes
- `ppn` - int
  Number of processors (CPU)
- `gpus` - int
  Number of processors (GPU)
- `walltime` - str
  Job walltime
- `LAMMPS_EXEC` - str
  Directory of the LAMMPS executable file

### write_pbs

```python
def write_pbs(output_dir)
```

Method to make the PBS job scheduler input file

**Arguments**:

- `output_dir` _str_ - Directory for the PBS job scheduler input file

**Returns**:

None
