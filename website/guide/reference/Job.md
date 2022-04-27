---
sidebar_label: Job
title: Job
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

