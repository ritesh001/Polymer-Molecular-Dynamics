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
               density=0.8,
               natoms_total=5000,
               natoms_per_chain=150,
               builder=pmd.EMC(force_field='pcff'))
```

### Polymer Chain Length Options

You can define the polymer chain length by `natoms_per_chain` (number of atoms in a polymer chain) as shown above. Alternatively, `mw_per_chain` (molecular weight of the polymer) or `ru_per_chain` (number of repeating units in a polymer chain) can be used, for example:

```python
s = pmd.System(smiles=smiles
               density=0.8,
               natoms_total=5000,
               #highlight-next-line
               mw_per_chain=1000,
               builder=pmd.EMC(force_field='pcff'))
```

```python
s = pmd.System(smiles=smiles,
               density=0.8,
               natoms_total=5000,
               #highlight-next-line
               ru_per_chain=25,
               builder=pmd.EMC(force_field='pcff'))
```

### Builder Options

For the builder, other than the `EMC` option shown above, we can also use the `PSP` option, for example:

```python
s = pmd.System(smiles=smiles,
               force_field='opls',
               density=0.8,
               natoms_total=5000,
               mw_per_chain=1000,
               #highlight-next-line
               builder=pmd.PSP(force_field='gaff2-gasteiger'))
```

The `EMC` builder has `force_field` options of **`pcff`**, **`opls-aa`**, **`opls-ua`**, and **`trappe`**.

And the `PSP` builder has `force_field` options of **`opls-cm1a`**, **`opls-lbcc`**, **`gaff2-gasteiger`**, and **`gaff2-am1bcc`**.

:::caution
If you want to use `PSP` as the builder, [PSP](https://github.com/Ramprasad-Group/PSP) and its dependencies have to be installed manually.
:::

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
               density=0.8,
               natoms_total=5000,
               natoms_per_chain=150,
               builder=pmd.EMC(force_field='pcff')
               # highlight-next-line
               data_fname='custom_data_file.data')

# highlight-next-line
s.write_data(output_dir='./custom_output_dir')
```
