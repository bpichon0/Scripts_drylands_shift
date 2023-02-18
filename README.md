# Code for running the model and analyses of vegetation dynamics in drylands

This repository contains the code used to perform simulations and the postprocessing for both main text and supplementary informations.

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**

All the code was made on R (v4.1.0) and Julia (1.7.3).

Here is the different steps to reproduce the figures:

## `Replicating the figures`

**0.** The file `Dryland_shift_function.R` and `Dryland_shift_Nspecies_function.jl` contain all functions necessary to run the code for the ODEs models and the spatially explicit one. 

**1.** To run all analyses made with two species, you cade run the codes:
    - `Dryland_shift_main.R` for the non-spatial simulations
    - For the spatially explicit simulations done with 2 species, run the code `Dryland_shift_Nspecies_main.jl` (**regions 1 to 5**). 

**3.** To run the analysis for N-species, you can run the file `Dryland_shift_Nspecies_main.jl` (**regions 6 to 10**). Note that simulations with 15 and 25 species using pair approximation takes some time run. It took about **1.5 days on a 25 cores CPU**. The code is nevertheness parallalized: you can change the number of cores used in the file. 
Simulatations are then post-processed using `Dryland_shift_main.R` script.

**4.** Once simulations are made and post-processed, all figures can be generated using the same file `Make_figures.R`. The code is organized by figures: a figure per chunk of code.



## `Exploring the model`

<p align="center">
    <img src="https://github.com/BPerezLamarque/ABDOMEN/blob/main/example/ABDOMEN.png" width="600">
</p>

<p align="center">
    <b>Figure 1: ABDOMEN: A comparative phylogenetic model for the dynamics of microbiota composition during host diversification.</b>
</p>


