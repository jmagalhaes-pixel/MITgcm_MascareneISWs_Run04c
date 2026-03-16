# MITgcm_MascareneISWs_Run04c
This repository contains the MITgcm configuration used to simulate internal solitary waves (ISWs) generated at the Mascarene Ridge in the Indian Ocean.
See **da Silva et al. (2015)** for more details.
The configuration is part of the **My Leap** project investigating the formation and evolution of ISWs using nonhydrostatic simulations.

---

## Repository Structure

code/
Model configuration and package options used during compilation.

input/
Runtime configuration files including physical parameters, grid, bathymetry, and forcing fields.

---

## Model

Simulations are performed using the **MIT General Circulation Model (MITgcm)** in nonhydrostatic mode.
Relevant MITgcm documentation:  
https://mitgcm.readthedocs.io

---

## Compilation

Compile MITgcm following the standard workflow:

./genmake2 -mods=code -optfile=OPTFILE
make depend
make


This produces the executable:

mitgcmuv


---

## Running the Simulation

Run the model using MPI:

mpirun -np 15 ./mitgcmuv > output.log 2>&1 &

---

## Experiment Description

Run04c explores the generation of internal solitary waves forced by tidal flow over realistic Mascarene Ridge bathymetry.

Key elements include:

- realistic ridge topography
- stratified ocean configuration
- barotropic tidal forcing and a barotropic body force in external_forcinf.F
- nonhydrostatic dynamics

---

## Notes

This repository contains **only the configuration files required to reproduce the experiment**.

Not included:

- compiled executable (`mitgcmuv`)
- build directory
- model output files
- restart files

---

## Author

Jorge M. R. Magalhães  
CIIMAR – Interdisciplinary Centre of Marine and Environmental Research  
Portugal


## Acknowledgments
Code originally developed by Professor Maarten Buijsman from the University of Southern Mississippi (USA) and used to simulate ISWs in the Mascarene Ridge of the Indian Ocean together with Professor José da Silva and Prof. Jorge Magalhaes from the University of Porto (Portugal).

Contact: jmagalhaes@fc.ul.pt

## References
da Silva, J. C. B., Buijsman, M. C., & Magalhães, J. M. R. (2015).
Internal solitary waves in the Mascarene Ridge region of the Indian Ocean.
Deep-Sea Research Part I: Oceanographic Research Papers, 99, 146–163.
https://doi.org/10.1016/j.dsr.2015.01.002



