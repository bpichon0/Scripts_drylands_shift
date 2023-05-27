include("Dryland_shift_Nspecies_function.jl")



#region : 1-- Illustration competitive exclusion 2 species (Fig 2) 




S_seq = collect(range(0, stop=0.2, length=2))[2] #the point where there is bistability in PA
c_seq = collect(range(0, stop=0.4, length=10))[10]
intra_comp_seq = [0.3]
param = Get_classical_param_2species()
size_landscape = 100
trajec_seq = ["Degradation" "Restoration"]
count = 1
scale_competition = ["global"]
disp_seq = [0.1]
d2 = zeros(length(S_seq) * length(c_seq) * length(scale_competition) * length(intra_comp_seq) * 2, 4) #Allocating

for disp in disp_seq
    for scale_comp in scale_competition
        for aii in intra_comp_seq
            for stress in S_seq
                param["S"] = stress
                for alpha_e in c_seq
                    param["alpha_e"] = alpha_e
                    for traj in trajec_seq
                        if traj == "Degradation"
                            ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])
                        else
                            ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.05, 0.05, 0.49, 0.5])
                        end


                        d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=5000, type_competition=scale_comp, save=false, burning=15000, N_snap=40,
                            name_save="")

                        #display(Plot_dynamics(d))
                        CSV.write("../Table/2_species/CA/Illustration/Dynamics_stress_" * repr(round(stress, digits=3)) *
                                  "_trajectory_" * traj * "_alpha0_" * repr(alpha_e) * ".csv", Tables.table(d), writeheader=false)

                        CSV.write("../Table/2_species/CA/Illustration/Landscape_stress_" * repr(round(stress, digits=3)) *
                                  "_trajectory_" * traj * "_alpha0_" * repr(alpha_e) * ".csv", Tables.table(state), writeheader=false)

                    end
                end
            end
        end
    end
end




#endregion


#region : 2-- Convergence CA & PA (SI fig)



S_seq = collect(range(0, stop=0.8, length=10)) #the point where there is bistability in PA
c_seq = collect(range(0, stop=0.4, length=10))
intra_comp_seq = [0.3]
param = Get_classical_param_2species()
size_landscape = 100
trajec_seq = ["Degradation" "Restoration"]
count = 1
scale_competition = ["global"]
disp_seq = [0.1]
d2 = zeros(length(S_seq) * length(c_seq) * length(scale_competition) * length(intra_comp_seq) * 2, 4) #Allocating

for disp in disp_seq
    for scale_comp in scale_competition
        for aii in intra_comp_seq
            for stress in S_seq
                param["S"] = stress
                for alpha_e in c_seq
                    param["alpha_e"] = alpha_e
                    for traj in trajec_seq
                        if traj == "Degradation"
                            ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])
                        else
                            ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.05, 0.05, 0.49, 0.5])
                        end


                        d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=3000, type_competition=scale_comp, save=false, burning=15000, N_snap=40,
                            name_save="")

                        #display(Plot_dynamics(d))
                        CSV.write("../Table/2_species/CA/Justifying_use_PA/Dynamics_stress_" * repr(round(stress, digits=3)) *
                                  "_trajectory_" * traj * "_alpha0_" * repr(alpha_e) * ".csv", Tables.table(d), writeheader=false)

                    end
                end
            end
        end
    end
end



#endregion


#region : 3-- Species pairs (Fig 3b)


param = Get_classical_param_2species()
size_landscape = 100
trajec_seq = ["Degradation"]
count = 1
scale_competition = ["global"]
disp_seq = collect(range(0, stop=1, length=12))

for disp in disp_seq
    param["S"] = 0
    param["alpha_e"] = 0.2
    param["delta"] = disp
    ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])


    d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=5000,
        type_competition="global", save=false, burning=15000, N_snap=40,
        name_save="")
    CSV.write("../Table/2_species/CA/Pairs" * "_delta_" * repr(disp) * ".csv", Tables.table(state), writeheader=false)

end




#endregion


#region : 4-- Clustering species dispersal (Fig 3c)


param = Get_classical_param_2species()
size_landscape = 100
trajec_seq = ["Degradation"]
count = 1
scale_competition = ["global"]
disp_seq = collect(range(0, stop=1, length=12))[[2 11]]

for disp in disp_seq
    param["S"] = 0.73
    param["alpha_e"] = 0.2
    param["delta"] = disp
    ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])


    d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=5000,
        type_competition="global", save=false, burning=15000, N_snap=40,
        name_save="")
    CSV.write("../Table/2_species/CA/Landscape_clustering" * "_delta_" * repr(disp) * ".csv", Tables.table(state), writeheader=false)

end




#endregion


#region : 5-- Vegetation along dispersal gradient (SI fig)



param = Get_classical_param_2species()
size_landscape = 100
scale_comp = ["global"]
param["alpha_e"] = 0.4
disp_seq = [0 0.1 0.2 0.3 0.7 1]

for disp in disp_seq
    param["delta"] = disp
    ini = Get_initial_lattice_2species(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])


    d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=3000, type_competition=scale_comp, save=false, burning=15000, N_snap=40,
        name_save="")

    #display(Plot_dynamics(d))
    CSV.write("../Table/2_species/CA/Dispersal_gradient_delta_" * repr(disp) * ".csv", Tables.table(state), writeheader=false)


end



#endregion


#region : 6-- 5species equal initial conditions PA (Fig 5)

N_sim = 100
Nsp = 5
a0_seq = [0 0.15 0.3]
tspan = (0.0, 200000)
branches = ["Degradation", "Restoration"]


N_tot_sim = N_sim * length(a0_seq) * 2

d = zeros(N_tot_sim, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 3)
dt = 1

#frac = rand(Nsp) #taking relative proportion of species
frac = zeros(Nsp) .+ 1 #same initial proportion

for a0 in a0_seq


    global p = Get_classical_param_Nspecies(N_species=Nsp,
        alpha_e=a0, scenario_trait="spaced", cintra=0.3)


    for branch_bifu in eachindex(1:2)

        if branch_bifu == 1 # Degradation
            state = Get_initial_state(Nsp=Nsp, type="equal", branch="Degradation", PA=true)
            S_seq = collect(range(0, 1, length=N_sim))

        else #Restoration
            state = Get_initial_state(Nsp=Nsp, type="equal", branch="Restoration", PA=true)
            S_seq = reverse(collect(range(0, 1, length=N_sim))) #to gradually decrease the level of stress

        end

        for stress in S_seq

            p[9] = stress

            prob = ODEProblem(PA_N_species, state, tspan, p)
            sol = solve(prob, callback=TerminateSteadyState(1e-8))
            d2 = Reorder_dynamics(sol)


            d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], a0, branch_bifu, stress)
            dt = dt + 1
        end
        print("Stress finished")
    end



end


CSV.write("../Table/N_species/PA/PA_equal_ini.csv", Tables.table(d), writeheader=false)

#endregion


#region : 7-- 5species equal initial conditions MF (Fig SI)

N_sim = 100
N_random_ini = 1

Nsp = 5
a0_seq = collect(range(0, 0.3, length=N_sim2))
tspan = (0.0, 200000)
branches = ["Degradation", "Restoration"]


N_tot_sim = N_sim * length(a0_seq) * 2

d = zeros(N_tot_sim, Nsp + 2 + 3)
dt = 1

#frac = rand(Nsp) #taking relative proportion of species
frac = zeros(Nsp) .+ 1 #same initial proportion

for a0 in a0_seq


    global p = Get_classical_param_Nspecies(N_species=Nsp,
        alpha_e=a0, scenario_trait="spaced", cintra=0.3)





    for branch_bifu in eachindex(1:2)

        if branch_bifu == 1 # Degradation
            state = Get_initial_state(Nsp=Nsp, type="equal", branch="Degradation", PA=false)
            S_seq = collect(range(0, 1, length=N_sim))

        else #Restoration
            state = Get_initial_state(Nsp=Nsp, type="equal", branch="Restoration", PA=false)
            S_seq = reverse(collect(range(0, 1, length=N_sim))) #to gradually decrease the level of stress

        end

        for stress in S_seq

            p[9] = stress

            prob = ODEProblem(MF_N_species, state, tspan, p)
            sol = solve(prob, callback=TerminateSteadyState(1e-8))
            d2 = Reorder_dynamics(sol)


            d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], a0, branch_bifu, stress)
            dt = dt + 1
        end
    end




end


CSV.write("../Table/N_species/MF/MF_equal_ini.csv", Tables.table(d), writeheader=false)





#endregion


#region : 8-- Nspecies example spatially explicit (Fig 5)


N_sim = 40
N_sim2 = 3
N_random_ini = 1

Nsp = 5
S_seq = collect(range(0, 0.6, length=2))
a0_seq = [0.15 0.3]
f_seq = collect(range(0, 0.9, length=3))[3]
tmax = 20000
branches = ["Restoration"]
tradeoff_seq = collect(range(0.25, 1, length=4))[4]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(tradeoff_seq) * 2 * N_random_ini


for random_ini in eachindex(1:N_random_ini)


    for a0 in a0_seq

        for tradeoff in tradeoff_seq

            global p = Get_classical_param_dict(N_species=Nsp,
                alpha_e=a0, scenario_trait="spaced", cintra=0.3)

            p["tau_leap"] = 0.5

            for facil in f_seq

                p["f"] = facil

                for branch_bifu in branches

                    if branch_bifu == 1 # Degradation
                        state = Get_initial_lattice_Nspecies(param=p, size_mat=100, branch="Degradation", type_ini="equal")

                    else #Restoration
                        state = Get_initial_lattice_Nspecies(param=p, size_mat=100, branch="Restoration", type_ini="equal")
                    end

                    for stress in S_seq

                        p["S"] = stress

                        d2, landscape = Gillespie_CA_N_species(; param=copy(p), landscape=copy(state), tmax=tmax, type_competition="global")

                        CSV.write("../Table/N_species/CA/Sim_Nrandom_" * repr(random_ini) *
                                  "_a0_" * repr(a0) * "_tradeoff_" * repr(tradeoff) * "_Nsp_" * repr(Nsp) * "_stress_" * repr(stress) *
                                  "_f_" * repr(facil) * ".csv", Tables.table(d2), writeheader=false)

                        CSV.write("../Table/N_species/CA/Landscape_Nrandom_" * repr(random_ini) *
                                  "_a0_" * repr(a0) * "_tradeoff_" * repr(tradeoff) * "_Nsp_" * repr(Nsp) * "_stress_" * repr(stress) *
                                  "_f_" * repr(facil) * ".csv", Tables.table(landscape), writeheader=false)

                    end

                end



            end
        end
    end
end

#endregion


#region : 9-- Full simulations for Fig 6 & SI 




using Distributed


addprocs(25, exeflags="--project=$(Base.active_project())")


@everywhere include("N_Species_Sim_functions.jl")

@everywhere begin
    using StatsBase, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings
    using BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames, SharedArrays
end







@everywhere function Run_sim_MF(N_random_ini)

    Nsp = 5
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = cat(collect(range(0.225, 0.35, length=length_aij_seq)), [0, 0.075, .15],dims=1)
    facil = 0.9
    tspan = (0.0, 1000000)

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Nsp + 2 + 5)
            dt = 1

            for direction in 1

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param_Nspecies(N_species=Nsp,
                    alpha_e=a0, scenario_trait="spaced", cintra=0.3, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=false)

                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, callback=TerminateSteadyState(1e-10))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("../Table/N_species/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_MF, 1:250)

@everywhere function Run_sim_MF(N_random_ini)

    Nsp = 15
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = cat(collect(range(0.225, 0.35, length=length_aij_seq)), [0, 0.075, .15],dims=1)
    facil = 0.9
    tspan = (0.0, 1000000)

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Nsp + 2 + 5)
            dt = 1

            for direction in 1

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param_Nspecies(N_species=Nsp,
                    alpha_e=a0, scenario_trait="spaced", cintra=0.3, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=false)


                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, callback=TerminateSteadyState(1e-10))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/N_species/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_MF, 1:250)




@everywhere function Run_sim_MF(N_random_ini)

    Nsp = 25
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = cat(collect(range(0.225, 0.35, length=length_aij_seq)), [0, 0.075, .15],dims=1)
    facil = 0.9
    tspan = (0.0, 1000000)

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Nsp + 2 + 5)
            dt = 1

            for direction in 1

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param_Nspecies(N_species=Nsp,
                    alpha_e=a0, scenario_trait="spaced", cintra=0.3, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=false)


                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, callback=TerminateSteadyState(1e-10))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/N_species/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_MF, 1:250)




@everywhere function Run_sim_PA(N_random_ini)

    Nsp = 5
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = a0_seq = cat(collect(range(0.225, 0.35, length=length_aij_seq)), [0, 0.075, .15],dims=1)
    facil = 0.9
    tspan = (0.0, 1000000)

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 5)
            dt = 1

            for direction in 1

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param_Nspecies(N_species=Nsp,
                    alpha_e=a0, scenario_trait="spaced", cintra=0.3, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=true)


                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(PA_N_species, state, tspan, p)
                    sol = solve(prob, callback=TerminateSteadyState(1e-10))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/N_species/PA/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_PA, 1:250)



@everywhere function Run_sim_PA(N_random_ini)

    Nsp = 15
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = a0_seq = cat(collect(range(0.225, 0.35, length=length_aij_seq)), [0, 0.075, .15],dims=1)
    facil = 0.9
    tspan = (0.0, 1000000)

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 5)
            dt = 1

            for direction in 1
                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param_Nspecies(N_species=Nsp,
                    alpha_e=a0, scenario_trait="spaced", cintra=0.3, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=true)

                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(PA_N_species, state, tspan, p)
                    sol = solve(prob, callback=TerminateSteadyState(1e-10))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/N_species/PA/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_PA, 1:250)




@everywhere function Run_sim_PA(N_random_ini)

    Nsp = 25
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = a0_seq = cat(collect(range(0.225, 0.35, length=length_aij_seq)), [0, 0.075, .15],dims=1)
    facil = 0.9
    tspan = (0.0, 1000000)

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 5)
            dt = 1

            for direction in 1

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param_Nspecies(N_species=Nsp,
                    alpha_e=a0, scenario_trait="spaced", cintra=0.3, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=true)

                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(PA_N_species, state, tspan, p)
                    sol = solve(prob, callback=TerminateSteadyState(1e-10))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/N_species/PA/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_PA, 1:250)












#endregion
