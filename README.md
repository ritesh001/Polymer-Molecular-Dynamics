# High-throughput-MD-simulations
This repo include Polymer Molecular Dynamics (PMD) package and examples for high-throughput MD simulations 

Properties
- Glass transition temperature (Tg) [[Tutorial]](https://github.com/Ramprasad-Group/High-throughput-MD-simulations/tree/main/tutorials/Tg)
- Gas diffusivity [[Tutorial]](https://github.com/Ramprasad-Group/High-throughput-MD-simulations/tree/main/tutorials/Gas_diffusivity)
- Solvent diffusivity [[Tutorial]](https://github.com/Ramprasad-Group/High-throughput-MD-simulations/tree/main/tutorials/Solvent_diffusivity)
- Thermal conductivity - upcoming
- Mechanical properties - upcoming
- Melting temperature (Tm)
- Viscosity
- Solubility (MC)

## Polymer Molecular Dynamics (PMD) package
[![License](https://img.shields.io/badge/license-MIT-blue.svg)](http://opensource.org/licenses/MIT)

PMD automatically builds LAMMPS data and input files for molecular dynamics simulations

## Example
Below is an example python script where we use PMD to generate LAMMPS data and input files for Tg measurement with a list of SMILES strings.
```python
import pmd

smiles = ['*CC*', '*CC(*)C', '*CC(*)CC', '*CC(*)CCC', '*CC(*)CCCC','*CC(*)c1ccccc1']

for s in smiles:
    lmp = pmd.Lammps(s, force_field='gaff2')
    lmp.add_procedure('minimization', min_style='cg')
    lmp.add_procedure('equilibration', Tfinal=600, Pfinal=1, Tmax=800, Pmax=49346.163)
    lmp.add_procedure('Tg_measurement', Tinit=600, Tfinal=100, Tinterval=25, step=1000000)
    lmp.write_input(output_dir=s)
                       
    job = pmd.Job(jobname=s, project='GT-rramprasad3-CODA20', nodes=2, ppn=24, walltime='48:00:00')
    job.write_pbs(output_dir=s)
```

A tutorial on using PMD to create simulations of a list of polymers can be found [here](https://github.com/Ramprasad-Group/High-throughput-MD-simulations/tree/main/Equilibration).

## Installation

```bash
git clone https://github.com/Ramprasad-Group/High-throughput-MD-simulations.git
pip install .
```

### Requirements
* [PSP](https://github.com/Ramprasad-Group/PSP) - for generating polymer topology
* [RDKit](https://www.rdkit.org/)

Note that, [PSP](https://github.com/Ramprasad-Group/PSP) and [RDKit](https://www.rdkit.org/) are required to be installed manually. NNMDKit requires PSP to create polymer structures. To configure the integration with PSP, environment variables are required to be addded:

For Ramprasad Group users on Gaanams and Tyrions, simply type in the command line:

```
echo "# PACKMOL PATH
export PACKMOL_EXEC='/home/appls/utility/packmol/packmol'

# pysimm PATHs
export PYTHONPATH=$PYTHONPATH:'/home/appls/utility/pysimm'
export PATH=$PATH:'/home/appls/utility/pysimm/bin'

# AmberTools PATH
export ANTECHAMBER_EXEC='/home/modules/anaconda3/bin/antechamber'" >> ~/.bashrc
```

