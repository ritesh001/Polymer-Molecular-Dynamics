# High Throughput Polymer MD Simulations [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Python toolkit and guides for molecular dynamics prediction of polymer properties.

You can find the documentation [on the website](https://high-throughput-pmd.netlify.app/api/overview).

Check out the [Getting Started](https://high-throughput-pmd.netlify.app/docs/getting-started/installation) page for a quick overview.

This package allows you to calculate polymer properties including:

- [Glass transition temperature (Tg)](http://high-throughput-pmd.netlify.app/docs/guides/glass-transition-temperature)
- [Gas diffusivity](http://high-throughput-pmd.netlify.app/docs/guides/gas-diffusivity)
- [Solvent diffusivity](http://high-throughput-pmd.netlify.app/docs/guides/solvent-diffusivity)
- Thermal conductivity - in-progress
- Mechanical properties - in-progress
- Viscosity - planned
- Melting temperature (Tm)
- Solubility (MC)

## Installation

```bash
pip install pmd
```

### Prerequisites

- [PSP](https://github.com/Ramprasad-Group/PSP) - for generating polymer topology
- [RDKit](https://www.rdkit.org/) - for getting molecule info such as molecular weight

Note that, [PSP](https://github.com/Ramprasad-Group/PSP) and [RDKit](https://www.rdkit.org/) are required to be installed manually.

## Example

Below is an example python script where we use PMD to generate LAMMPS data and input files for Tg measurement with a list of SMILES strings.

```python
import pmd

for smiles in ['*CC*', '*CC(*)CC','*CC(*)c1ccccc1']:
    # Define system specs and make the data file
    s = pmd.System(smiles=smiles, force_field='opls', density=0.8, natoms_total=5000, natoms_per_chain=150)
    s.write_data(output_dir=smiles)

    # Customize LAMMPS simulation and make the input file
    lmp = pmd.Lammps(s)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(pmd.Equilibration())
    lmp.add_procedure(pmd.TgMeasurement())
    lmp.write_input(output_dir=smiles)

    # Create job scheduler file
    job = pmd.Job(jobname=smiles, project='Your-project-id', nodes=2, ppn=24, walltime='48:00:00')
    job.write_pbs(output_dir=smiles)
```
