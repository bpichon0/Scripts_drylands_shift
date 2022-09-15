include("Network_shift_functions.jl")


#region :  Step 1 : Exploration

#niche of species
param = Get_params(cg=0.2, cl=0.2, fs=0.5, S=0.5, z=4, nb_sp=4)

landscape = Get_initial_lattice(param=param, size_mat=25)
Plot_landscape(landscape)


time_step=1:1000

param = [0.02,0.1,.9,0.8, 0.05,.1,1,0.1,.4,.1,0,.1,1,4]
function Gillepsie_tau_leeping_Julia(param,landscape,time_step)
    
    r,d,f,beta,m,e,emax,cintra,cinter1,cinter2,S,delta,dt,z=param

  
  
  
end
  
  
  