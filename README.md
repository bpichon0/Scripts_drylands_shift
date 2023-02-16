# Code for running the model and analyses of vegetation dynamics in drylands

This repository contains the code used to perform simulations and the postprocessing for both main text and supplementary informations.

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**

All the code was made on R (v4.1.0) and Julia (1.7.3).

Here is the different steps to reproduce the figures:

**1.** The file `Dryland_shift_function.R` and `N_Species_Sim_functions.jl` contain all functions necessary to run the code for ODE and spatial model respectively. 

**2.** To run all analyses made with two species, you cade run the codes:
    - `Dryland_shift_main.R` for the non-spatial simulations
    - For the spatially explicit simulations done with 2 species, run the code `2_species.jl`. 

**3.** To run the analysis for N-species, you can run the file `N_Species_Sim_main.jl`. Note that simulations with 15 and 25 species using pair approximation takes some time run. The code is nevertheness parallalized. It took about **1.5 day on a 25 cores CPU**. 

**4.** Once simulations are made, all figures are made using the same file `Make_figures.R`.

