include("N_Species_Sim_functions.jl")




#region : Step 1 : Varying the facilitation and competition strength


N_sim = 100
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
N_random_ini = 1
tspan = (0.0, 10000)
d = zeros(N_sim * N_random_ini * 3 * 3, Nsp + 5)
dt = 1
branches = ["Restoration", "Degradation"]
facilitation_seq = [0 0.3 0.9]

for facil in facilitation_seq
    for a0 in [0 0.15 0.3]#[0 0.15 0.3]
        for h in [1]
            for community_compo in ["spaced"]#[0.2 0.8 "spaced"]
                for random_ini in 1:N_random_ini

                    global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                        alpha_0=a0, scenario_trait=community_compo, cintra=0.3, h=h)
                    state = Get_initial_state(Nsp=Nsp, type="equal", branch="Degradation")

                    for stress in S_seq
                        p["f"] = facil
                        p["S"] = stress
                        prob = ODEProblem(MF_N_species, state, tspan, p)
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)
                        if community_compo == "spaced"
                            community_compo = 0
                        end

                        d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, a0, facil)
                        dt = dt + 1
                    end
                    print(random_ini)
                end
            end
        end
    end
end
CSV.write("../Table/N_species/MF/Varying_competition_facilitation.csv", Tables.table(d), writeheader=false)

#endregion

#region : Step 2 : Number of tipping points 


N_sim = 50
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
N_random_ini = 50
tspan = (0.0, 10000)
d = zeros(N_sim * N_random_ini * 4 * 4, Nsp + 6)
dt = 1
branches = ["Degradation", "Restoration"]
facilitation_seq = collect(range(0, 0.9, length=4))
competition_seq = collect(range(0, 0.3, length=4))


for facil in facilitation_seq
    for a0 in competition_seq
        for h in [1]
            #for branch in branches[1]
            for random_ini in 1:N_random_ini

                global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=h)
                state = Get_initial_state(Nsp=Nsp, type="random", branch="Degradation")

                for stress in S_seq
                    p["f"] = facil
                    p["S"] = stress
                    prob = ODEProblem(MF_N_species, state, tspan, p)
                    sol = solve(prob, Tsit5())
                    d2 = Reorder_dynamics(sol)

                    d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], stress, a0, facil, random_ini)
                    dt = dt + 1
                end
                print(random_ini)
            end
            #end
        end
    end
end
CSV.write("../Table/N_species/MF/Number_tipping.csv", Tables.table(d), writeheader=false)


#endregion

#region : testing PA N species 

p = Get_classical_param(N_species=10, type_interaction="classic",
    alpha_0=0.3, scenario_trait="spaced", cintra=0.3, h=1)
state = Get_initial_state(Nsp=10, type="equal", branch="Degradation", PA=true)
p[9] = 0.25
p[10] = 0.1
tspan = (0, 20000)
prob = ODEProblem(PA_N_species, state, tspan, p)
sol = solve(prob, Tsit5())

d2 = Reorder_dynamics(sol)
Plot_dynamics(d=d2, Nsp=15, PA=true)


#endregion



#region : 15sp Final simulation 2 = equal ini

N_sim = 100
N_sim2 = 3
N_random_ini = 1

Nsp = 15
a0_seq = collect(range(0, 0.3, length=N_sim2))
f_seq = collect(range(0, 0.9, length=3))
tspan = (0.0, 15000)
branches = ["Degradation", "Restoration"]
delta_seq = collect(range(0.1, 0.9, length=2))[1]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(delta_seq) * 2 * N_random_ini


for random_ini in 1:N_random_ini

    #frac = rand(Nsp) #taking relative proportion of species
    frac = zeros(Nsp) .+ 1 #same initial proportion

    for a0 in a0_seq

        for disp in delta_seq

            global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1)

            p["delta"] = disp

            for facil in f_seq

                p["f"] = facil

                d = zeros(N_sim * 2, Nsp + 8)
                dt = 1

                for branch_bifu in 1:2

                    if branch_bifu == 1 # Degradation
                        state = (frac / sum(frac)) * 0.8 #normalizing the sum to 80% of vegetation cover
                        push!(state, 0.1)
                        push!(state, 0.1)
                        S_seq = collect(range(0, 1, length=N_sim))

                    else #Restoration
                        state = (frac / sum(frac)) * 0.01 #normalizing the sum to 1% of vegetation cover
                        push!(state, 0.49)
                        push!(state, 0.5)
                        S_seq = reverse(collect(range(0, 1, length=N_sim))) #to gradually decrease the level of stress

                    end

                    for stress in S_seq

                        p["S"] = stress

                        prob = ODEProblem(MF_N_species, state, tspan, p)
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)


                        d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], random_ini, a0, disp, facil, branch_bifu, stress)
                        dt = dt + 1
                    end
                end
                CSV.write("../Table/N_species/MF/Big_sim_Nsp/Constant_ini/Sim_Nrandom_" * repr(random_ini) *
                          "_a0_" * repr(a0) *
                          "_delta_" * repr(disp) *
                          "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)



            end
        end
    end
end






#endregion


#region : Nsp Final simulation 1 = random ini 






N_sim = 50
N_sim2 = 3
N_random_ini = 150


a0_seq = collect(range(0, 0.3, length=N_sim2))[3]

f_seq = collect(range(0, 0.9, length=3))[3]
branches = ["Degradation", "Restoration"]
delta_seq = collect(range(0.1, 0.9, length=2))[1]
Nsp_seq = [5, 10, 15, 20, 25, 30, 35]
runining_time = [10000, 15000, 20000, 20000, 25000, 25000, 30000]

N_tot_sim = N_sim * length(a0_seq) * length(f_seq) * length(delta_seq) * 2 * N_random_ini

for Nsp_index in 1:Nsp_seq


    Nsp = Nsp_seq[Nsp_index]

    tspan = (0.0, runining_time[Nsp_index])

    for random_ini in 1:N_random_ini

        frac = rand(Nsp) #taking relative proportion of species
        #frac = zeros(Nsp) .+ 1

        for a0 in a0_seq

            for disp in delta_seq

                global param = Get_classical_param(N_species=Nsp, type_interaction="classic",
                    alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1)


                for facil in f_seq


                    d = zeros(N_sim * 2, Nsp + 8)
                    dt = 1

                    for branch_bifu in 1:length(branches)

                        if branch_bifu == 1 # Degradation
                            state = (frac / sum(frac)) * 0.8 #normalizing the sum to 80% of vegetation cover
                            push!(state, 0.1)
                            push!(state, 0.1)
                            S_seq = collect(range(0, 1, length=N_sim))

                        else #Restoration
                            state = (frac / sum(frac)) * 0.01 #normalizing the sum to 1% of vegetation cover
                            push!(state, 0.49)
                            push!(state, 0.5)
                            S_seq = reverse(collect(range(0, 1, length=N_sim))) #to gradually decrease the level of stress

                        end

                        for stress in S_seq

                            param[9] = stress

                            prob = ODEProblem(MF_N_species, state, tspan, param)
                            sol = solve(prob, Tsit5())
                            d2 = Reorder_dynamics(sol)
                            #display(Plot_dynamics(d=d2, Nsp=Nsp))

                            d[dt, :] = push!(d2[size(d2)[1], 1:(param[1]+2)], random_ini, a0, disp, facil, branch_bifu, stress)
                            dt = dt + 1
                        end
                    end
                    CSV.write("../Table/N_species/MF/" * repr(Nsp) * "_sp/Sim_Nrandom_" * repr(random_ini) *
                              "_a0_" * repr(a0) *
                              "_delta_" * repr(disp) *
                              "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)



                end
            end
        end
    end
end









#endregion



#region : Pa N_species 




N_sim = 40
N_sim2 = 3
N_random_ini = 1

Nsp = 15
S_seq = collect(range(0, 0.82, length=N_sim))
a0_seq = collect(range(0.15, 0.4, length=30))
f_seq = collect(range(0, 0.9, length=3))[3]
tspan = (0.0, 30000)
branches = ["Degradation"]
tradeoff_seq = collect(range(0.25, 1, length=4))[4]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(tradeoff_seq) * 2 * N_random_ini


for random_ini in 1:N_random_ini

    #frac = rand(Nsp) #taking relative proportion of species
    frac = zeros(Nsp) .+ 1 #same initial proportion

    for a0 in a0_seq

        for tradeoff in tradeoff_seq

            global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1, trade_off=tradeoff)

            p[16] = tradeoff

            for facil in f_seq

                p[4] = facil

                d = zeros(N_sim * 2, Int(((Nsp^2 + 7 * Nsp) / 2) + 5) + 6)
                dt = 1

                for branch_bifu in 1:2

                    if branch_bifu == 1 # Degradation
                        state = Get_initial_state(Nsp=Nsp, type="equal", branch="Degradation", PA=true)

                    else #Restoration
                        state = Get_initial_state(Nsp=Nsp, type="equal", branch="Restoration", PA=true)
                    end

                    for stress in S_seq

                        p[9] = stress

                        prob = ODEProblem(MF_N_species, state, tspan, p)
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)


                        d[dt, :] = push!(d2[size(d2)[1], 1:(Int(((Nsp^2 + 7 * Nsp) / 2) + 5))], random_ini, a0, tradeoff, facil, branch_bifu, stress)
                        dt = dt + 1
                    end
                end
                CSV.write("../Table/N_species/PA/Minimal_comp/Sim_Nrandom_" * repr(random_ini) *
                          "_a0_" * repr(a0) * "_tradeoff_" * repr(tradeoff) * "_Nsp_" * repr(Nsp) *
                          "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)



            end
        end
    end
end


#


# Nspecies cellular automaton


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
