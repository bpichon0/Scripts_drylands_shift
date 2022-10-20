include("N_Species_Sim_functions.jl")



#region : Step 1 : Testing MF model on 15 species MF
#Testing with 15 species : functional diversity, species diversity, stability, community weighted mean
N_sim = 100
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
tspan = (0.0, 10000)
d = zeros(N_sim * 3 * 2, Nsp + 5)
dt = 1

for a0 in [0 0.15 0.3]
    for h in [1 1.5]

        global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
            alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=h)
        state = Get_initial_state(param=p, type="equal", branch="Degradation")

        for stress in S_seq
            p["S"] = stress
            prob = ODEProblem(MF_N_species, state, tspan, p)
            sol = solve(prob, Tsit5())
            d2 = Reorder_dynamics(sol)
            d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, a0, h)
            dt = dt + 1
        end

    end
end
CSV.write("../Table/N_species/MF/15_species.csv", Tables.table(d), writeheader=false)
CSV.write("../Table/N_species/MF/Trait_15sp.csv", Tables.table(p["trait"]), writeheader=false)

#endregion

#region : Step 2 : 15 species functional diversity for varying competition and facilitation strength ----
N_sim = 100
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
tspan = (0.0, 10000)
d = zeros(N_sim * 10, Nsp + 6)
dt = 1

for a_ii in [0.1 0.3]
    if a_ii == 0.1
        a0_seq = [0 0.1]
    else
        a0_seq = [0 0.1 0.2]
    end
    for a0 in a0_seq
        for f0 in [0.5 0.9]

            global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=a_ii)
            state = Get_initial_state(param=p, type="equal", branch="Degradation")
            p["f"] = f0
            for stress in S_seq
                p["S"] = stress
                prob = ODEProblem(MF_N_species, state, tspan, p)
                sol = solve(prob, Tsit5())
                d2 = Reorder_dynamics(sol)
                d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, a_ii, a0, f0)
                dt = dt + 1
            end
        end
    end
end
CSV.write("../Table/N_species/MF/15_species_random_ini.csv", Tables.table(d), writeheader=false)





#endregion
