---
sidebar_position: 2
title: Creating data files
---

To create a data file, first we need to define our system.

## Define the system

Instantiate a System object with your polymer SMILES string and other specs:

```python
import pmd

s = pmd.System(smiles=smiles,
               force_field='gaff2',
               density=0.8,
               natoms_total=5000,
               natoms_per_chain=150)
```

### Polymer Chain Length Options

You can define the polymer chain length by `natoms_per_chain` (number of atoms in a polymer chain) as shown above. Alternatively, `mw_per_chain` (molecular weight of the polymer) or `ru_per_chain` (number of repeating units in a polymer chain) can be used, for example:

```python
s = pmd.System(smiles=smiles,
               force_field='gaff2',
               density=0.8,
               natoms_total=5000,
               #highlight-next-line
               mw_per_chain=1000)
```

```python
s = pmd.System(smiles=smiles,
               force_field='gaff2',
               density=0.8,
               natoms_total=5000,
               #highlight-next-line
               ru_per_chain=25)
```

### Force field Options

For the force field, other than the `gaff2` option shown above, we can also use the `opls` option, for example:

```python
s = pmd.System(smiles=smiles,
               #highlight-next-line
               force_field='opls',
               density=0.8,
               natoms_total=5000,
               mw_per_chain=1000)
```

## Generate LAMMPS data file

Once the System object is created, you can call the `write_data` function:

```python
s.write_data()
```

This outputs your data file with a default name of `'data.lmps'` at the directory that you run the script

### Customize data file name and directory

```python
import pmd

s = pmd.System(smiles=smiles,
               force_field='gaff2',
               density=0.8,
               natoms_total=5000,
               natoms_per_chain=150,
               # highlight-next-line
               data_fname='custom_data_file.data')

# highlight-next-line
s.write_data(output_dir='./custom_output_dir')
```
