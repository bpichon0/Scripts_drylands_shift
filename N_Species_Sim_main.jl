include("N_Species_Sim_functions.jl")



#region : testing PA/MF N species 

Nsp = 25
p = Get_classical_param(N_species=Nsp, type_interaction="classic",
    alpha_0=0.3, scenario_trait="spaced", cintra=0.3, h=1)



@time for i in 1:10
    Random.seed!(i)
    state = Get_initial_state(Nsp=Nsp, type="random", branch="Degradation", PA=false)

    p[9] = 0
    p[10] = 0.1
    tspan = (0, 30000)
    prob = ODEProblem(MF_N_species, state, tspan, p)

    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
    d2 = Reorder_dynamics(sol)
    display(Plot_dynamics(d=d2, Nsp=Nsp, PA=false))
end;

Nsp = 25
p = Get_classical_param(N_species=Nsp, type_interaction="classic",
    alpha_0=0.3, scenario_trait="spaced", cintra=0.3, h=1)
tspan = (0, 30000)

Random.seed!(213)
state = Get_initial_state(Nsp=Nsp, type="random", branch="Degradation", PA=true)
prob = ODEProblem(PA_N_species, state, tspan, p)

@time sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7));
d2 = Reorder_dynamics(sol)
Plot_dynamics(d=d2, Nsp=Nsp, PA=true)



#endregion


#region : 15sp equal initial conditions PA (Fig 5)

N_sim = 100
N_sim2 = 3
N_random_ini = 1

Nsp = 15
a0_seq = collect(range(0, 0.3, length=N_sim2))
f_seq = collect(range(0, 0.9, length=3))[3]
tspan = (0.0, 35000)
branches = ["Degradation", "Restoration"]
delta_seq = collect(range(0.1, 0.9, length=2))[1]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(delta_seq) * 2 * N_random_ini

d = zeros(N_tot_sim, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 3)
dt = 1

for random_ini in 1:1
    #frac = rand(Nsp) #taking relative proportion of species
    frac = zeros(Nsp) .+ 1 #same initial proportion

    for a0 in a0_seq

        for disp in delta_seq

            global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1)

            for facil in f_seq



                for branch_bifu in 1:2

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
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)


                        d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], a0, branch_bifu, stress)
                        dt = dt + 1
                    end
                end


            end
        end
    end
end

CSV.write("../Table/N_species/PA/PA_equal_ini.csv", Tables.table(d), writeheader=false)

#endregion

#region : 15sp equal initial conditions MF (Fig SI)

N_sim = 100
N_sim2 = 3
N_random_ini = 1

Nsp = 15
a0_seq = collect(range(0, 0.3, length=N_sim2))
f_seq = collect(range(0, 0.9, length=3))[3]
tspan = (0.0, 25000)
branches = ["Degradation", "Restoration"]
delta_seq = collect(range(0.1, 0.9, length=2))[1]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(delta_seq) * 2 * N_random_ini

d = zeros(N_tot_sim, Nsp + 2 + 3)
dt = 1

for random_ini in 1:1
    #frac = rand(Nsp) #taking relative proportion of species
    frac = zeros(Nsp) .+ 1 #same initial proportion

    for a0 in a0_seq

        for disp in delta_seq

            global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1)


            for facil in f_seq



                for branch_bifu in 1:2

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
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)


                        d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], a0, branch_bifu, stress)
                        dt = dt + 1
                    end
                end


            end
        end
    end
end

CSV.write("../Table/N_species/MF/MF_equal_ini.csv", Tables.table(d), writeheader=false)





#endregion



#region: 5 species example spatially explicit (Fig 5)


N_sim = 40
N_sim2 = 3
N_random_ini = 1

Nsp = 5
S_seq = collect(range(0, 0.7, length=3))
a0_seq = collect(range(0.1, 0.4, length=3))
f_seq = collect(range(0, 0.9, length=3))[3]
tmax = 50000
branches = ["Degradation"]
tradeoff_seq = collect(range(0.25, 1, length=4))[4]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(tradeoff_seq) * 2 * N_random_ini


for random_ini in 1:N_random_ini


    for a0 in a0_seq

        for tradeoff in tradeoff_seq

            global p = Get_classical_param_dict(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1)


            for facil in f_seq

                p["f"] = facil

                for branch_bifu in 1:2

                    if branch_bifu == 1 # Degradation
                        state = Get_initial_lattice(param=p, size_mat=100, branch="Degradation", type_ini="equal")

                    else #Restoration
                        state = Get_initial_lattice(param=p, size_mat=100, branch="Restoration", type_ini="equal")
                    end

                    for stress in S_seq

                        p["S"] = stress

                        d2, landscape = Gillespie_CA_N_species(; param=copy(p), landscape=copy(state), tmax=tmax, type_competition="global")
                        display(Plot_dynamics(d=d2, Nsp=5))
                        CSV.write("../Table/N_species/CA/Sim_Nrandom_" * repr(random_ini) *
                                  "_a0_" * repr(a0) * "_tradeoff_" * repr(tradeoff) * "_Nsp_" * repr(Nsp) * "_stress_" * repr(stress) *
                                  "_f_" * repr(facil) * ".csv", Tables.table(d2), writeheader=false)

                        CSV.write("../Table/N_species/CA/Landscape_Nrandom_" * repr(random_ini) *
                                  "_a0_" * repr(a0) * "_tradeoff_" * repr(tradeoff) * "_Nsp_" * repr(Nsp) * "_stress_" * repr(stress) *
                                  "_f_" * repr(facil) * ".csv", Tables.table(landscape), writeheader=false)


                        display(Plot_landscape(landscape))

                    end
                end



            end
        end
    end
end

#endregion


#region: 15 species example spatially explicit (Fig 5)


Nsp = 15
S_seq = collect(range(0, 0.7, length=3))[1]
a0_seq = collect(range(0, 0.4, length=3))[1]
tmax = 20000
branches = ["Degradation"]

global p = Get_classical_param_dict(N_species=Nsp, type_interaction="classic",
    alpha_0=0, scenario_trait="spaced", cintra=0.3, h=1)

state = Get_initial_lattice(param=p, size_mat=100, branch="Degradation", type_ini="equal")


d2, landscape = Gillespie_CA_N_species(; param=copy(p), landscape=copy(state), tmax=tmax, type_competition="global")
display(Plot_dynamics(d=d2, Nsp=15))

display(Plot_landscape(landscape))






#endregion



#region: Full simulations for Fig 6 & SI 




using Distributed


mkdir("./Table")

mkdir("./Table/MF")
mkdir("./Table/PA")

species_seq = [5 15 25]
for k in species_seq
    mkdir("./Table/MF/" * repr(k) * "_sp")
    mkdir("./Table/PA/" * repr(k) * "_sp")
end

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
    a0_seq = collect(range(0.225, 0.35, length=length_aij_seq))
    facil = 0.9

    if Nsp == 5
        tspan = (0.0, 20000)
    elseif Nsp == 15
        tspan = (0.0, 30000)
    elseif Nsp == 25
        tspan = (0.0, 40000)
    else
        tspan = (0.0, 40000)
    end

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Nsp + 2 + 5)
            dt = 1

            for direction in 1:2

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=false)

                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_MF, 1:150)

@everywhere function Run_sim_MF(N_random_ini)

    Nsp = 15
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = collect(range(0.225, 0.35, length=length_aij_seq))
    facil = 0.9

    if Nsp == 5
        tspan = (0.0, 20000)
    elseif Nsp == 15
        tspan = (0.0, 30000)
    elseif Nsp == 25
        tspan = (0.0, 40000)
    else
        tspan = (0.0, 40000)
    end

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Nsp + 2 + 5)
            dt = 1

            for direction in 1:2

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=false)


                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_MF, 1:150)




@everywhere function Run_sim_MF(N_random_ini)

    Nsp = 25
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = collect(range(0.225, 0.35, length=length_aij_seq))
    facil = 0.9

    if Nsp == 5
        tspan = (0.0, 20000)
    elseif Nsp == 15
        tspan = (0.0, 30000)
    elseif Nsp == 25
        tspan = (0.0, 40000)
    else
        tspan = (0.0, 40000)
    end

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Nsp + 2 + 5)
            dt = 1

            for direction in 1:2

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=false)


                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Nsp+2)], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_MF, 1:150)




@everywhere function Run_sim_PA(N_random_ini)

    Nsp = 5
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = collect(range(0.225, 0.35, length=length_aij_seq))
    facil = 0.9

    if Nsp == 5
        tspan = (0.0, 20000)
    elseif Nsp == 15
        tspan = (0.0, 30000)
    elseif Nsp == 25
        tspan = (0.0, 40000)
    else
        tspan = (0.0, 40000)
    end

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 5)
            dt = 1

            for direction in 1:2

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=true)


                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(PA_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/PA/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_PA, 1:150)



@everywhere function Run_sim_PA(N_random_ini)

    Nsp = 15
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = collect(range(0.225, 0.35, length=length_aij_seq))
    facil = 0.9

    if Nsp == 5
        tspan = (0.0, 20000)
    elseif Nsp == 15
        tspan = (0.0, 30000)
    elseif Nsp == 25
        tspan = (0.0, 40000)
    else
        tspan = (0.0, 40000)
    end

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 5)
            dt = 1

            for direction in 1:2
                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=true)

                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(PA_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/PA/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_PA, 1:150)




@everywhere function Run_sim_PA(N_random_ini)

    Nsp = 25
    N_sim_S = 100
    length_aij_seq = 6


    S_seq = collect(range(0, 0.9, length=N_sim_S))
    a0_seq = collect(range(0.225, 0.35, length=length_aij_seq))
    facil = 0.9

    if Nsp == 5
        tspan = (0.0, 20000)
    elseif Nsp == 15
        tspan = (0.0, 30000)
    elseif Nsp == 25
        tspan = (0.0, 40000)
    else
        tspan = (0.0, 40000)
    end

    for random_ini in N_random_ini

        for a0 in a0_seq
            d = zeros(N_sim_S * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 5)
            dt = 1

            for direction in 1:2

                if direction == 1
                    branches = "Degradation"
                else
                    branches = "Restoration"
                end

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=1)


                Random.seed!(random_ini)
                state = Get_initial_state(Nsp=Nsp, type="random", branch=branches, PA=true)

                for stress in S_seq

                    p[9] = stress

                    prob = ODEProblem(PA_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5(), callback=TerminateSteadyState(1e-7))
                    d2 = Reorder_dynamics(sol)


                    d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, facil, direction, stress)
                    dt = dt + 1
                end
            end
            CSV.write("./Table/PA/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                      "_a0_" * repr(a0) * "_Nsp_" * repr(Nsp) *
                      "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)

        end
    end
end


pmap(Run_sim_PA, 1:150)












#endregion
