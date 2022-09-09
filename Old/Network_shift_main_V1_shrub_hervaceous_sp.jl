cd("C:\\Users\\Benoi\\OneDrive\\Documents\\Phd\\Network_shift\\Scripts")
include("Network_shift_functions.jl")


#region :  Step 1 : Exploration 

#niche of species
param=Get_params(cg=0.2,cl=0.2,fs=0.5,S=2,z=4,frac_shrub=0.5)
Plot_species_niche(param)
savefig("../Figures/niche_species.pdf")


land_ini=Get_initial_lattice(frac=[0.9,0.1,0],size_mat=25)
Plot_landscape(land_ini)


#bifu with all parameters equal except niche center
A_seq=0:0.1:1.2  ; timeseq=1:1:1000
d=Array{Float64}(undef,length(A_seq),4)
param=Get_params(cg=0.2,cl=0.2,fs=0.5,S=2,z=4,frac_shrub=0.5)
param=Set_equal_param(param)

for (i,Aloop) in enumerate(A_seq)
    land_ini=Get_initial_lattice()
    merge!(param,Dict("A"=>Aloop))
    d2,landscape=Run_CA_shift(timeseq,copy(param),copy(land_ini))
    mean_densities=[mean(d2[900:size(d2,1),x]) for x in 2:size(d2,2)]
    d[i,:]=pushfirst!(mean_densities,Aloop)
    #Plot_dynamics(d2,"time (t)")
    #Plot_landscape(landscape)
end

Plot_dynamics(d,"Aridity (A)")
savefig("../Figures/Comunity_gradient.pdf")


CSV.write("./Figures/test.csv",  Tables.table(d), writeheader=false)

param=Get_params(cg=0.2,cl=0.2,fs=0.5,S=2,z=4,frac_shrub=0.5)
land_ini=Get_initial_lattice()
merge!(param,Dict("A"=>Aloop))
d2,landscape=Run_CA_shift(timeseq,copy(param),copy(land_ini))
Plot_dynamics(d2,"time (t)")
Plot_landscape(landscape)







u0 = Float64[0.4 , 0.4*0.4 , 0.4 , 0.4*0.4 , 0.1 , 0.4*0.4 ,  0.1*0.4 , 0.1*0.4 , 0.1*0.1]
tspan = (0.0,1000.0)
t=0:0.1:1000
param=Get_params(cg=0.2,cl=0.2,fs=0.5,S=2,z=4,frac_shrub=0.5)
prob = ODEProblem(PA_shift,u0,tspan,param)
sol = solve(prob,Rosenbrock23(),saveat=t)
sol_dyn=Reorder_dynamics(sol)
sol_dyn=sol_dyn[:,1:5]
Plot_dynamics(sol_dyn,"time t")





