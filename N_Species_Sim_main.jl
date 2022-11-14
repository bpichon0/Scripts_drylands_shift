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
                    state = Get_initial_state(Nsp=Nsp, type="random", branch=branches[branch_bifu])

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


#region : Step 2 : Varying the facilitation and competition strength


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

#region : Step 3 : Changing the community composition


N_sim = 100
Nsp = 15
S_seq = collect(range(0, 1, length=N_sim))
N_random_ini = 1
tspan = (0.0, 10000)
d = zeros(N_sim * N_random_ini * 3 * 3, Nsp + 5)
dt = 1
facilitation_seq = [0 0.3 0.9]
community_comp_seq = collect(range(0, 1, length=10))


for facil in facilitation_seq
    for a0 in [0 0.15 0.3]#[0 0.15 0.3]
        for h in [1]
            for community_compo in community_comp_seq
                for random_ini in 1:N_random_ini

                    global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                        alpha_0=a0, scenario_trait=community_compo, cintra=0.3, h=h)
                    state = Get_initial_state(Nsp=Nsp, type="equal", branch="Degradation")

                    for stress in S_seq
                        p["S"] = stress
                        prob = ODEProblem(MF_N_species, state, tspan, p)
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)

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

#region :  Number of tipping points 


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

p = Get_classical_param(N_species=2, type_interaction="classic",
    alpha_0=0.1, scenario_trait="spaced", cintra=0.3, h=1)
state = Get_initial_state(Nsp=2, type="equal", branch="Degradation", PA=true)
p["S"] = 0.1
tspan = (0, 1000)
prob = ODEProblem(PA_N_species, state, tspan, p)
sol = solve(prob, Tsit5())
d2 = Reorder_dynamics(sol)
Plot_dynamics(d=d2, Nsp=2)

#endregion



#region : 15sp Final simulation 1 = random ini

N_sim = 50
N_sim2 = 3
N_random_ini = 40

Nsp = 15
a0_seq = collect(range(0, 0.3, length=N_sim2))
f_seq = collect(range(0, 0.9, length=3))
tspan = (0.0, 20000)
branches = ["Degradation", "Restoration"]
delta_seq = collect(range(0.1, 0.9, length=2))[1]


N_tot_sim = N_sim * N_sim2 * length(f_seq) * length(delta_seq) * 2 * N_random_ini


for random_ini in 1:N_random_ini

    frac = rand(Nsp) #taking relative proportion of species
    #frac = zeros(Nsp) .+ 1

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
                CSV.write("../Table/N_species/MF/Big_sim_Nsp/Random_ini/Sim_Nrandom_" * repr(random_ini) *
                          "_a0_" * repr(a0) *
                          "_delta_" * repr(disp) *
                          "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)



            end
        end
    end
end






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


#region : 30sp Testing random ini 




#endregion



N_sim = 40
N_sim2 = 3
N_random_ini = 100

Nsp = 35
a0_seq = collect(range(0, 0.3, length=N_sim2))[3]
f_seq = collect(range(0, 0.9, length=3))[3]
tspan = (0.0, 20000)
branches = ["Degradation", "Restoration"]
delta_seq = collect(range(0.1, 0.9, length=2))[1]


N_tot_sim = N_sim * length(a0_seq) * length(f_seq) * length(delta_seq) * 2 * N_random_ini


for random_ini in 1:100

    frac = rand(Nsp) #taking relative proportion of species
    #frac = zeros(Nsp) .+ 1

    for a0 in a0_seq

        for disp in delta_seq

            global p = Get_classical_param(N_species=Nsp, type_interaction="classic",
                alpha_0=a0, scenario_trait="spaced", cintra=0.3, h=1)

            p["delta"] = disp

            for facil in f_seq

                p["f"] = facil

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

                        p["S"] = stress

                        prob = ODEProblem(MF_N_species, state, tspan, p)
                        sol = solve(prob, Tsit5())
                        d2 = Reorder_dynamics(sol)


                        d[dt, :] = push!(d2[size(d2)[1], 1:(p["Nsp"]+2)], random_ini, a0, disp, facil, branch_bifu, stress)
                        dt = dt + 1
                    end
                end
                CSV.write("../Table/N_species/MF/35_sp/Sim_Nrandom_" * repr(random_ini) *
                          "_a0_" * repr(a0) *
                          "_delta_" * repr(disp) *
                          "_f_" * repr(facil) * ".csv", Tables.table(d), writeheader=false)



            end
        end
    end
end
