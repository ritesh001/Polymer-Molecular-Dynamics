# High Throughput Polymer MD Simulations [![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

Python toolkit and guides for molecular dynamics prediction of polymer properties.

- Quick overview - Check out the [Introduction page](https://high-throughput-pmd.netlify.app/docs/intro)
- Documentation - You can find it [on the website](https://high-throughput-pmd.netlify.app/api/overview).

## Guides and scripts

You can calculate various polymer properties from MD simulations with this package, including:

- Glass transition temperature (Tg) - [[Guide](http://high-throughput-pmd.netlify.app/docs/guides/glass-transition-temperature)] [[Scripts](https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations/tree/main/scripts/Tg)]
- Gas diffusivity - [[Guide](http://high-throughput-pmd.netlify.app/docs/guides/gas-diffusivity)] [[Scripts](https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations/tree/main/scripts/Gas_diffusivity)]
- Solvent diffusivity - [[Guide](http://high-throughput-pmd.netlify.app/docs/guides/solvent-diffusivity)] [[Scripts](https://github.com/Ramprasad-Group/High-Throughput-Polymer-MD-Simulations/tree/main/scripts/Solvent_diffusivity)]
- Thermal conductivity - In-progress
- Mechanical properties - In-progress
- Solubility - Planned
- Viscosity - Planned
- Melting temperature (Tm) - Planned

## Example

Below is an example python script where we use PMD to generate LAMMPS data and input files for Tg measurement with a list of SMILES strings.

```python
import pmd

for smiles in ['*CC*', '*CC(*)CC','*CC(*)c1ccccc1']:
    # Define system specs and make the data file
    s = pmd.System(smiles=smiles, force_field='opls', density=0.8,
                   natoms_total=5000, natoms_per_chain=150)
    s.write_data(output_dir=smiles)

    # Customize LAMMPS simulation and make the input file
    lmp = pmd.Lammps(read_data_from=s)
    lmp.add_procedure(pmd.Minimization())
    lmp.add_procedure(pmd.Equilibration())
    lmp.add_procedure(pmd.TgMeasurement())
    lmp.write_lammps(output_dir=smiles)

    # Create job scheduler file
    job = pmd.Torque(run_lammps=lmp, jobname=smiles, project='Your-project-id',
                     nodes=2, ppn=24, walltime='48:00:00')
    job.write_job(output_dir=smiles)
```

## Installation

```bash
pip install pmd
```

### Prerequisites

- [PSP](https://github.com/Ramprasad-Group/PSP) - for generating polymer topology
- [RDKit](https://www.rdkit.org/) - for getting molecule info such as molecular weight

PSP and RDKit are required to be installed for the package to fully function.
