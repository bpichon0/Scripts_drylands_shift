# Code for running the model and analyses of vegetation dynamics in drylands

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**

This folder contains all code necessary to replicate the analysis in the main text and in supplementary. The file Dryland_shift_function.R and N_Species_Sim_functions.jl contain all functions necessary to run the code for ODE and spatial model respectively.

All *non spatial* simulations done with 2 species can be found in Dryland_shift_main.R and the figures can be generated using Make_figures.R file.

*Spatially explicit* simulations with 2 species can be found in 2_speicies.jl and the figures can be generated using Make_figures.R file.

Finally, all simulations (with explicit space or not) are done in the file N_Species_Sim_main.jl. Note that bifurcation diagrams with 15 and 25 species using pair approximation takes a lot of time (~2 weeks on 15 cores). 
