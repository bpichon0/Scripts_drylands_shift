include("Network_shift_functions.jl")


#region :  Step 1 : Exploration

#niche of species
param = Get_params(cg=0.2, cl=0.2, fs=0.5, S=0.5, z=4, nb_sp=4)

landscape = Get_initial_lattice(param=param, size_mat=25)
Plot_landscape(landscape)
