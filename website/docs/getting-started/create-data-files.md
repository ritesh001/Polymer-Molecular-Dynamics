---
sidebar_position: 2
title: Create data files
---

To create a data file, first we need to define our system.

Instantiate a System object with your polymer SMILES string and other specs:

```python
import pmd

s = pmd.System(smiles=smiles,
               force_field='gaff2',
               density=0.8,
               natoms_total=5000,
               natoms_per_chain=150)
```

Once the System object is created, you can call the `write_data` function:

```python
s.write_data(output_dir="your-output-dir", data_fname = 'amor_data.lmps')
```
