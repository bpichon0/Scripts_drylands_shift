# Code for running the model and analyses of vegetation dynamics in drylands

Contact: Benoît Pichon, **benoit.pichon0@gmail.com**


This repository contains the code used to perform simulations and the postprocessing for both main text and supplementary informations.
All the code was made on R (*v4.1.0*) and Julia (*1.7.3*).

Here is the different steps to reproduce the figures:


## `Installing Julia dependancies`

To install julia dependancies go into your working directory, press "]" and enter "activate .". Once activating the local folder, all dependancies can be loaded using "] instantiate".

## `Installing R dependancies`

To install the pacakges needed for the analyses and create the folder architecture used to save the data sets, please load the file `Dryland_shift_function.R` using: 

```R
source(Dryland_shift_function.R)
```

## `Replicating the figures`

**0.** The file `Dryland_shift_function.R` and `Dryland_shift_Nspecies_function.jl` contain all functions necessary to run the code for the ODEs models and the spatially explicit one. 

**1.** To run all analyses made with two species, you can run the codes:
    - `Dryland_shift_main.R` for the non-spatial simulations
    - For the spatially explicit simulations done with 2 species, run the code `Dryland_shift_Nspecies_main.jl` (**regions 1 to 5**). 

**3.** To run the analysis for N-species, you can run the file `Dryland_shift_Nspecies_main.jl` (**regions 6 to 9**). Note that simulations with 15 and 25 species using pair approximation takes some time run. It took about **1.5 days on a 25 cores CPU**. The code is nevertheness parallalized: you can change the number of cores used in the file. 
Simulatations are then post-processed using `Dryland_shift_main.R` script.

**4.** Once simulations are made and post-processed, all figures can be generated using the same file `Make_figures.R`. The code is organized by figures: a figure per chunk of code. If you want to generate only a part of the results, we indicate in each chunk in Make_figure.R file the script and step to simulate in order to replicate the given figure. 



## `Exploring the model`


### Model description

<p align="center">
    <img src="https://github.com/bpichon0/Scripts_drylands_shift/blob/master/Example/Model.jpg" width="1000">
</p>

<p align="center">
    <b>Figure 1: Trait based plant community model in arid ecosystems.</b>
</p>

The model describes the temporal changes of a landscape of $n \times n$ sites.
Each site can be in one of $N+2$ states: colonized by the plant $i$ ($+_i$), with a fertile soil ($0$) or in a degraded soil state ($-$).
The transitions between the different states are probabilistic and depend on the local neighborhood of a focal site (*e.g.* the local vegetation density) and the trait of the plants (Fig. 1a).

Each plant species $i$ is characterized by a strategy modelled as a position $\psi_i$ along the competitivity <-> stress-tolerance trade-off.
This strategy determines whether the species $i$ facilitates its local environment as well as its tolerance to competition and abiotic stress (*e.g.,* aridity).




### A toy example


Let's take as an example a community of 5 plant species equally distributed along the strategy trade-off axis.
We begin by the spatially explicit dynamics:

```julia

include("Dryland_shift_Nspecies_function.jl") #loading the functions


Nsp = 5 #number of species
p = Get_classical_param_dict(N_species=Nsp, alpha_e=0.1, #level of interspecific competition
            scenario_trait="spaced", #species traits are equally spaced along the trade-off axis. Other methods are possible: you can give the composition of the community in term of strategies or consider randomly distributed species ("random") along this axis
             cintra=0.3)

state = Get_initial_lattice_Nspecies(param=p, size_mat=100, branch="Degradation", type_ini="equal") #intial lattice

d2, landscape = Gillespie_CA_N_species(; param=copy(p), landscape=copy(state), tmax=1000, type_competition="global") #Simulate the dynamics

Plot_dynamics_Nspecies(d=d2, Nsp=5) # temporal dynamics

Plot_landscape_Nspecies(landscape) # spatial organization of the vegetation

```

<p align="center">
    <img src="https://github.com/bpichon0/Scripts_drylands_shift/blob/master/Example/Temporal_dynamics.svg" width="600">
</p>

<p align="center">
    <b>Figure 2: Temporal dynamics of the 5 species community.</b>
    Green lines are the more competitive species while blue lines are the most stress-tolerant (and therefore facilitating) species 
</p>



<p align="center">
    <img src="https://github.com/bpichon0/Scripts_drylands_shift/blob/master/Example/Landscape_asymp_state.svg" width="600">
</p>

<p align="center">
    <b>Figure 3: Spatial organization of the vegetation at asymptotic state.</b>
</p>



We can also simulate the temporal dynamics of this 5-species community with the pair-approximation model (deterministic system of ODEs).



```julia


Nsp = 5

p = Get_classical_param_Nspecies(N_species=Nsp,
    alpha_e=0.1, scenario_trait="spaced", cintra=0.3)

state = Get_initial_state(Nsp=Nsp, type="equal", branch="Degradation", PA=true)
prob = ODEProblem(PA_N_species, state, (0, 5000), p)
sol = solve(prob, callback=TerminateSteadyState(1e-8)) #solving the pair approximation model

dyn = Reorder_dynamics(sol) #to reorder the output of DifferentialEquation package
Plot_dynamics_Nspecies(d=dyn, Nsp=Nsp) 

# Note that you can also run the mean-field model by changin PA = false in the Get_initial_state function and
# calling MF_N_species in the ODEProblem instead of PA_N_species. 
```



<p align="center">
    <img src="https://github.com/bpichon0/Scripts_drylands_shift/blob/master/Example/PA_dynamics.svg" width="600">
</p>

<p align="center">
    <b>Figure 3: Temporal dynamics with the pair-approximation model.</b>
</p>



