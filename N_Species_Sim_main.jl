include("N_Species_Sim_functions.jl")



#region : Step 1 : Testing MF model on 15 species MF
N_sim = 100
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
N_random_ini = 10
tspan = (0.0, 10000)
d = zeros(N_sim * N_random_ini * 1 * 1 * 2, Nsp + 7)
dt = 1
branches = ["Restoration", "Degradation"]

for branch_bifu in 1:2
    for a0 in [0.3]#[0 0.15 0.3]
        for h in [1]
            for community_compo in ["spaced"]#[0.2 0.8 "spaced"]
                for random_ini in 1:N_random_ini

                    global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                        alpha_0=a0, scenario_trait=community_compo, cintra=0.3, h=h)
                    state = Get_initial_state(param=p, type="random", branch=branches[branch_bifu])

                    for stress in S_seq
                        p["S"] = stress
                        prob = ODEProblem(MF_N_species, state, tspan, p)
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)
                        if community_compo == "spaced"
                            community_compo = 0
                        end

                        d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, a0, community_compo, random_ini, branch_bifu)
                        dt = dt + 1
                    end
                    print(random_ini)
                end
            end
        end
    end
end
CSV.write("../Table/N_species/MF/15_species.csv", Tables.table(d), writeheader=false)

#endregion
