include("N_Species_Sim_functions.jl")






#region : Step 1 : Testing MF model on 4 species MF
N_sim = 100
S_seq = collect(range(0, 1, length=N_sim))
dt = 1
tspan = (0.0, 2000)
Scenario_comp = ["low_inter" "nested"]
Rela_comp = [1 2.5 4]
d = zeros(length(Scenario_comp) * length(Rela_comp) * N_sim, 9)
u = 0
for scena in Scenario_comp
    u += 1
    for rela_c in Rela_comp

        p = Get_classical_param(N_species=4, type_interaction=scena, relative_competition=rela_c,
            scenario_trait="spaced")
        state = Get_initial_state(param=p, type="equal", branch="Degradation")

        for stress in S_seq

            p["S"] = stress
            p["cg"] = 0.01

            prob = ODEProblem(MF_N_species, state, tspan, p)
            sol = solve(prob, Tsit5())
            d2 = Reorder_dynamics(sol)
            d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, rela_c, u)
            dt = dt + 1

        end
        p
    end
end

CSV.write("../Table/N_species/MF/4_species.csv", Tables.table(d), writeheader=false)


#endregion

#region : Step 2 : 15 species with varying community composition MF ----

N_sim = 100
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
tspan = (0.0, 20000)

Scenario_comp = ["low_inter" "nested"]
Rela_comp = [1 2.5 4]
d = zeros(length(Scenario_comp) * length(Rela_comp) * N_sim, Nsp + 5)
branch_seq = ["Degradation" "Restoration"]


for branch_sim in branch_seq
    for community_compo in [0.2 0.5 0.8] #initial community composition
        d = zeros(length(Scenario_comp) * length(Rela_comp) * N_sim, Nsp + 5)
        dt = 1
        u = 0
        for scena in Scenario_comp
            u += 1
            for rela_c in Rela_comp
                Random.seed!(Int64(100 * community_compo))
                global p = Get_classical_param(N_species=Nsp, type_interaction=scena, relative_competition=rela_c,
                    scenario_trait=community_compo)
                state = Get_initial_state(param=p, type="equal", branch=branch_sim)

                for stress in S_seq

                    p["S"] = stress
                    p["cg"] = 0.01

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5())
                    d2 = Reorder_dynamics(sol)
                    d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, rela_c, u)
                    dt = dt + 1

                end
            end
        end
        CSV.write("../Table/N_species/MF/15_species_" * branch_sim * "_frac_facilitator_" * repr(100 * community_compo) * ".csv", Tables.table(d), writeheader=false)
        CSV.write("../Table/N_species/MF/TRAITS_15_species_" * branch_sim * "_frac_facilitator_" * repr(100 * community_compo) * ".csv", Tables.table(p["trait"]), writeheader=false)
        CSV.write("../Table/N_species/MF/Aij_15_species_frac_" * branch_sim * "_facilitator_" * repr(100 * community_compo) * ".csv", Tables.table(p["alpha"]), writeheader=false)
    end
end
#endregion

#region : Step 2bis : Role of facilitation 
N_sim = 60
Nsp = 15
S_seq = collect(range(0, 0.9, length=N_sim))
tspan = (0.0, 20000)

Scenario_comp = ["nested"]
Rela_comp = [1 2.5 4]
d = zeros(length(Scenario_comp) * length(Rela_comp) * N_sim, Nsp + 5)
branch_seq = ["Degradation"]


for branch_sim in branch_seq
    for community_compo in [0.5] #initial community composition
        d = zeros(length(Scenario_comp) * length(Rela_comp) * N_sim, Nsp + 5)
        dt = 1
        u = 0
        for scena in Scenario_comp
            u += 1
            for facil in [0.5 0.75 1]
                rela_c = 2.5
                Random.seed!(Int64(100 * community_compo))
                global p = Get_classical_param(N_species=Nsp, type_interaction=scena, relative_competition=rela_c,
                    scenario_trait=community_compo)
                state = Get_initial_state(param=p, type="equal", branch=branch_sim)

                p["f"] = facil

                for stress in S_seq

                    p["S"] = stress
                    p["cg"] = 0.01

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5())
                    d2 = Reorder_dynamics(sol)
                    d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, facil, u)
                    dt = dt + 1

                end
            end
        end
    end
end
CSV.write("../Table/N_species/MF/Role_facilitation/Role_facilitation_dynamics.csv", Tables.table(d), writeheader=false)
CSV.write("../Table/N_species/MF/Role_facilitation/Role_facilitation_traits.csv", Tables.table(p["trait"]), writeheader=false)


#endregion

#region : Step 3 : 15 species : alternative stable states ----

N_sim = 80
N_initial_state = 20
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
tspan = (0.0, 20000)
Scenario_comp = ["low_inter" "nested"]
Rela_comp = [1 2.5 4]
d = zeros(N_initial_state * N_sim, Nsp + 6)

for community_compo in [0.2 0.5 0.8] #initial community composition
    u = 0
    for scena in Scenario_comp
        u += 1
        for rela_c in Rela_comp
            dt = 1
            Random.seed!(Int64(100 * community_compo))
            global p = Get_classical_param(N_species=Nsp, type_interaction=scena, relative_competition=rela_c,
                scenario_trait=community_compo)

            for stress in S_seq
                for n_ini in 1:N_initial_state

                    state = Get_initial_state(param=p, type="random", branch="Degradation")
                    p["S"] = stress
                    p["cg"] = 0.01

                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5())
                    d2 = Reorder_dynamics(sol)
                    d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, rela_c, u, n_ini)
                    dt = dt + 1

                end #end loop random initial states 
            end #stress loop

            #write the tables

            name_sim = "15_species_frac_facilitator_" * repr(100 * community_compo) * "_RelaC_" * repr(rela_c) * "_Scena_" * scena * ".csv"
            CSV.write("../Table/N_species/MF/ASS_15_species/" * name_sim, Tables.table(d), writeheader=false)
            CSV.write("../Table/N_species/MF/ASS_15_species/TRAITS_" * name_sim, Tables.table(p["trait"]), writeheader=false)
            CSV.write("../Table/N_species/MF/ASS_15_species/Aij_name_sim_" * name_sim, Tables.table(p["alpha"]), writeheader=false)
            d = zeros(N_initial_state * N_sim, Nsp + 6)

        end #competition loop 
    end # inter/intra loop
end
#endregion

#region : Step 4 : 15 species with varying community composition CA ----


Nsp = 15
tsim = 10000
community_compo = 0.2
Random.seed!(Int64(100 * community_compo))
p = Get_classical_param(N_species=Nsp, type_interaction="low_inter", relative_competition=1,
    scenario_trait=community_compo)
p["tau_leap"] = 1
Ncells = 100

for stress in [0 0.25 0.5 0.75]
    p["S"] = stress
    ini = Get_initial_lattice(param=copy(p), size_mat=Ncells, branch="Degradation")
    @time d, state = Gillespie_CA_N_species(param=copy(p), landscape=copy(ini), tmax=tsim)
    #Plot_landscape(state)
    #Plot_dynamics(d=d, Nsp=Nsp, name_x_axis="Time (t)")

    CSV.write("../Table/N_species/CA/Final_state_IBM_S_" * repr(stress) * ".csv", Tables.table(state), writeheader=false)
    CSV.write("../Table/N_species/CA/Dynamics_IBM_S_" * repr(stress) * ".csv", Tables.table(d), writeheader=false)
    CSV.write("../Table/N_species/CA/IBM_traits_sp_S_" * repr(stress) * ".csv", Tables.table(p["trait"]), writeheader=false)
    z_score = Compute_z_score_pairs(; final_state=state, Nsp=15, N_null=500)[:, :, 1]
    CSV.write("../Table/N_species/CA/Z_score_IBM_occurrence_S_" * repr(stress) * ".csv", Tables.table(z_score), writeheader=false)

end



#endregion
