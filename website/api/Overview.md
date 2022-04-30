---
sidebar_position: 1
title: Overview
---

## Core module

- [**System**](/api/core/System): Objects that stores system specs such as polymer SMILES and force field of the simulation
  - [write_data](/api/core/System#write_data): Method to make LAMMPS data file (which contains coordinates and force field parameters)
- [**Lammps**](/api/core/Lammps): Objects related to configure LAMMPS settings and writing LAMMPS input files
  - [add_procedure](/api/core/Lammps#add_procedure): Method to add simulation procedure
  - [write_input](/api/core/Lammps#write_input): Method to make LAMMPS input files
- [**Procedure**](/api/core/Procedure): A group of objects that are simulation procedures that can be linked to a Lammps object
  - [Minimization](/api/core/Procedure#minimization-objects): Procedure object that enables performing energy minimization of the system
  - [Equilibrium](/api/core/Procedure#equilibration-objects): Procedure object that enables performing a 21-step equilibration process
  - [NPT](/api/core/Procedure#npt-objects): Procedure object that enables performing a simulation under NPT ensemble
  - [NVT](/api/core/Procedure#nvt-objects): Procedure object that enables performing a simulation under NVT ensemble
  - [TgMeasurement](/api/core/Procedure#tgmeasurement-objects): Procedure object that enables performing a cooling simulation for Tg measurement
- [**Job**](/api/core/Job): Objects related to creating job scheduler-relevant files for submitting jobs on supercomputers
  - [write_pbs](/api/core/Job#write_pbs): Method to make the PBS job scheduler input file
  - [write_slurm](/api/core/Job#write_slurm): Method to make the Slurm job scheduler input file

## Postprocess module

- [**Analysis**](/api/postprocessing/Analysis): A suite of methods for calculating system properties
  - [calculate_Tg](/api/postprocessing/Analysis#calculate_tg): Method to calculate glass transition temperature based on the result from the TgMeasurement Procedure
- [**TrajectoryReader**](/api/postprocessing/TrajectoryReader): Methods to read in a lammps trajectory file