include("Network_shift_functions.jl")


# region : Step 0 : Comparing computational speed of two species functions 
p = Get_classical_param(N_species=2, type_interaction="low_inter", relative_competition=3, scenario_trait="spaced")
nb_cell = 25
landscape = Get_initial_lattice(param=p, size_mat=nb_cell)

@time d2, landscape = Gillespie_CA_N_species_V1(param=p, landscape=landscape, tmax=2000)
Plot_dynamics(d=d2, Nsp=4, name_x_axis="time")



#endregion

p = Get_classical_param()
landscape = Get_initial_lattice()

@time d2, landscape = Gillespie_tau_leeping(param=p, landscape=landscape, time=2000)
Plot_dynamics(d2)



N_sim = 80
S_seq = collect(range(0, 1, length=N_sim))
dt = 1
p = Get_classical_param(N_species=4, type_interaction="low_inter", relative_competition=3)
state = Get_initial_state(param=p)
tspan = (0.0, 2000)

d = zeros(N_sim, 7)
for stress in S_seq

    p["S"] = stress
    p["cg"] = 0.01
    prob = ODEProblem(MF_N_species, state, tspan, p)
    sol = solve(prob, Tsit5())
    d2 = Reorder_dynamics(sol)
    d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress)
    dt = dt + 1

end

Plot_dynamics(d=d, Nsp=p["Nsp"], name_x_axis="Stress (S)")
