using StatsBase, RCall, Plots, Random, Agents, DifferentialEquations, BenchmarkTools, CSV, Tables, Distributions

#region : 0-- Functions, CA, Gillespie


function Get_classical_param()

    r = 0.01
    d = 0.025
    f = 0.9
    beta = 1
    m = 0.15
    e = 0.1
    cg = 0
    alpha_0 = 0.3
    cintra = 0.3
    S = 0
    delta = 0.1
    z = 4
    tau_leap = 0.5
    h = 1

    return Dict{String,Any}(
        "r" => r,
        "d" => d,
        "f" => f,
        "beta" => beta,
        "m" => m,
        "e" => e,
        "cg" => cg,
        "alpha_0" => alpha_0,
        "cintra" => cintra,
        "S" => S,
        "delta" => delta,
        "z" => z,
        "tau_leap" => tau_leap,
        "h" => h
    )
end

function Get_initial_lattice(; frac=[0.4, 0.4, 0.1, 0.1], size_mat=25)
    ini_vec = sample([1, 2, 0, -1], Weights(frac), size_mat * size_mat)
    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end







"Function for faster computing of CA using Gillespie Tau leeping method"
function Gillespie_tau_leeping(; landscape, param, time, type_competition, save=false, burning=1500, N_snap=25, name_save="")
    n_save = 1

    d2 = zeros(time, 5) #Allocating
    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    cg = param["cg"]
    cintra = param["cintra"]
    alpha_0 = param["alpha_0"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]
    tau_leap = param["tau_leap"]
    h = param["h"]


    rules_change = transpose([0 0 1 2 -1 0; 1 2 0 0 0 -1])

    #Global densities
    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    nb_cell = size(landscape)[1]

    #Allocating 
    Rate_landscape = zeros(nb_cell, nb_cell, 6)


    if type_competition == "global" #competition occurs globally

        for t in 1:time


            @rput landscape
            R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
            R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

            @rget neigh_1
            @rget neigh_2

            Rate_landscape[:, :, 1] .= @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - ((cintra * rho_1 + alpha_0 * h * (1 + exp(-1)) * rho_2)))
            Rate_landscape[:, :, 1] .= Rate_landscape[:, :, 1] .* (landscape .== 0)

            Rate_landscape[:, :, 2] .= @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - ((cintra * rho_2 + alpha_0 * rho_1)))
            Rate_landscape[:, :, 2] .= Rate_landscape[:, :, 2] .* (landscape .== 0)

            # calculate regeneration, degradation & mortality rate
            Rate_landscape[:, :, 3] .= m .* (landscape .== 1)
            Rate_landscape[:, :, 4] .= m .* (landscape .== 2)
            Rate_landscape[:, :, 5] .= @.(r + f * neigh_1 / z) .* (landscape .== -1)
            Rate_landscape[:, :, 6] .= d .* (landscape .== 0)

            Rate_landscape[findall(Rate_landscape .< 0)] .= 0 #to avoid problems with propensity

            #calculate propensity

            propensity = [sum(Rate_landscape[:, :, k]) for k in 1:size(Rate_landscape)[3]]

            nb_events = map(x -> rand(Poisson(x)), propensity * tau_leap)

            for event in 1:length(nb_events) #for each type of events
                patches = findall(landscape .== rules_change[event, 1])

                if nb_events[event] != 0 && length(patches) > nb_events[event]
                    landscape[wsample(patches, Rate_landscape[patches, event], nb_events[event])] .= rules_change[event, 2]
                end
            end

            rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
            rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
            rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
            rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

            @views d2[t, :] = [t rho_1 rho_2 rho_f rho_d]
            Rate_landscape = zeros(nb_cell, nb_cell, 6)

            if save && t > burning && t % ((time - burning) / N_snap) == 0
                CSV.write(name_save * "_nsave_" * repr(n_save) * ".csv", Tables.table(landscape), writeheader=false)
                n_save += 1
            end


        end


    else #competition occurs locally


        for t in 1:time

            @rput landscape
            R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
            R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

            @rget neigh_1
            @rget neigh_2

            Rate_landscape[:, :, 1] .= @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - (cg * (rho_1 + rho_2) + (cintra * (neigh_1 / z) + alpha_0 * (1 + h * exp(-1)) * (neigh_2 / z))))
            Rate_landscape[:, :, 1] .= Rate_landscape[:, :, 1] .* (landscape .== 0)

            Rate_landscape[:, :, 2] .= @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - (cg * (rho_1 + rho_2) + (cintra * (neigh_2 / z) + alpha_0 * (neigh_1 / z))))
            Rate_landscape[:, :, 2] .= Rate_landscape[:, :, 2] .* (landscape .== 0)

            # calculate regeneration, degradation & mortality rate
            Rate_landscape[:, :, 3] .= m .* (landscape .== 1)
            Rate_landscape[:, :, 4] .= m .* (landscape .== 2)
            Rate_landscape[:, :, 5] .= @.(r + f * neigh_1 / z) .* (landscape .== -1)
            Rate_landscape[:, :, 6] .= d .* (landscape .== 0)

            Rate_landscape[findall(Rate_landscape .< 0)] .= 0 #to avoid problems with propensity

            #calculate propensity

            propensity = [sum(Rate_landscape[:, :, k]) for k in 1:size(Rate_landscape)[3]]

            nb_events = map(x -> rand(Poisson(x)), propensity * tau_leap)

            for event in 1:length(nb_events) #for each type of events
                patches = findall(landscape .== rules_change[event, 1])

                if nb_events[event] != 0 && length(patches) > nb_events[event]
                    landscape[wsample(patches, Rate_landscape[patches, event], nb_events[event])] .= rules_change[event, 2]
                end
            end

            rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
            rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
            rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
            rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

            @views d2[t, :] = [t rho_1 rho_2 rho_f rho_d]
            Rate_landscape = zeros(nb_cell, nb_cell, 6)

            if save && t > burning && t % ((time - burning) / N_snap) == 0
                CSV.write(name_save * "_nsave_" * repr(n_save) * ".csv", Tables.table(landscape), writeheader=false)
                n_save += 1
            end



        end
    end
    return d2, landscape
end


function Ca_2_species_global(; landscape, param)


    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded


    # Neighbors :

    #using simcol package from R 
    @rput landscape
    R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

    @rget neigh_1
    @rget neigh_2

    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    cintra = param["cintra"]
    alpha_0 = param["alpha_0"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]

    colonization1 = @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - (cintra * rho_1 + alpha_0 * (1 + exp(-1)) * rho_2))
    colonization2 = @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - (cintra * rho_2 + alpha_0 * rho_1))
    # calculate regeneration, degradation & mortality rate
    death = m
    regeneration = (@.(r + f * neigh_1 / z)) .* (landscape .== -1)
    degradation = d

    # Apply rules
    rnum = reshape(rand(length(landscape)), Int64(sqrt(length(landscape))), Int64(sqrt(length(landscape))))# one random number between 0 and 1 for each cell
    landscape_update = copy(landscape)

    ## New vegetation
    landscape_update[findall((landscape .== 0) .& (rnum .<= colonization1))] .= 1
    landscape_update[findall((landscape .== 0) .& (rnum .> colonization1) .& (rnum .<= colonization1 .+ colonization2))] .= 2

    ## New fertile
    landscape_update[findall((landscape .== 1) .& (rnum .<= death))] .= 0
    landscape_update[findall((landscape .== 2) .& (rnum .<= death))] .= 0
    landscape_update[findall((landscape .== -1) .& (rnum .<= regeneration))] .= 0

    ## New degraded 
    landscape_update[findall((landscape .== 0) .& (rnum .> colonization1 .+ colonization2) .& (rnum .<= (colonization1 .+ colonization2 .+ degradation)))] .= -1

    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    return rho_1, rho_2, rho_f, rho_d, landscape_update

end


function Ca_2_species_local(; landscape, param)


    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded


    # Neighbors :

    #using simcol package from R 
    @rput landscape
    R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

    @rget neigh_1
    @rget neigh_2

    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    cg = param["cg"]
    cintra = param["cintra"]
    alpha_0 = param["alpha_0"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]

    colonization1 = @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - (cg * (rho_1 + rho_2) + (cintra * (neigh_1 / z) + alpha_0 * (1 + exp(-1)) * (neigh_2 / z))))
    colonization2 = @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - (cg * (rho_1 + rho_2) + (cintra * (neigh_2 / z) + alpha_0 * (neigh_1 / z))))
    # calculate regeneration, degradation & mortality rate
    death = m
    regeneration = @.(r + f * neigh_1 / z)
    degradation = d

    # Apply rules
    rnum = reshape(rand(length(landscape)), Int64(sqrt(length(landscape))), Int64(sqrt(length(landscape))))# one random number between 0 and 1 for each cell
    landscape_update = copy(landscape)

    ## New vegetation
    landscape_update[findall((landscape .== 0) .& (rnum .<= colonization1))] .= 1
    landscape_update[findall((landscape .== 0) .& (rnum .> colonization1) .& (rnum .<= colonization1 .+ colonization2))] .= 2

    ## New fertile
    landscape_update[findall((landscape .== 1) .& (rnum .<= death))] .= 0
    landscape_update[findall((landscape .== 2) .& (rnum .<= death))] .= 0
    landscape_update[findall((landscape .== -1) .& (rnum .<= regeneration))] .= 0

    ## New degraded 
    landscape_update[findall((landscape .== 0) .& (rnum .> colonization1 .+ colonization2) .& (rnum .<= (colonization1 .+ colonization2 .+ degradation)))] .= -1

    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    return rho_1, rho_2, rho_f, rho_d, landscape_update

end






function Run_CA_2_species(; time, param, landscape, save, name_save="", burning=1500, N_snap=25, type_competition)

    d = Array{Float64}(undef, time + 1, 5) #Allocating

    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    d[1, :] = [1 rho_1 rho_2 rho_f rho_d] #the dataframe 

    n_save = 1
    if type_competition == "global" #global competition
        @inbounds for k in 1:time

            rho_1, rho_2, rho_f, rho_d, landscape = Ca_2_species_global(landscape=copy(landscape), param=param)
            @views d[k+1, :] = [k + 1 rho_1 rho_2 rho_f rho_d]

            if save && k > burning && k % round(((time - burning) / N_snap)) == 0
                CSV.write(name_save * "_nsave_" * repr(n_save) * ".csv", Tables.table(landscape), writeheader=false)
                n_save += 1

            end
        end
    else #local competition

        @inbounds for k in 1:time

            rho_1, rho_2, rho_f, rho_d, landscape = Ca_2_species_local(landscape=copy(landscape), param=param)
            @views d[k+1, :] = [k + 1 rho_1 rho_2 rho_f rho_d]

            if save && k > burning && k % round(((time - burning) / N_snap)) == 0
                CSV.write(name_save * "_nsave_" * repr(n_save) * ".csv", Tables.table(landscape), writeheader=false)
                n_save += 1

            end
        end

    end


    return d, landscape

end
function Plot_dynamics(d)

    plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
    plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
    plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
    plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")


end

function Plot_landscape(landscape)
    colGRAD = cgrad([colorant"#696969", colorant"#D8CC7B", colorant"blue", colorant"#ACD87B"])
    heatmap(landscape, yflip=true, fill=true, c=colGRAD)
end


#endregion



#


#region : 1-- Illustration competitive exclusion 2 species (Fig 2) 




S_seq = collect(range(0, stop=0.2, length=2))[2] #the point where there is bistability in PA
c_seq = collect(range(0, stop=0.4, length=10))[10]
intra_comp_seq = [0.3]
param = Get_classical_param()
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
                    param["alpha_0"] = alpha_e
                    for traj in trajec_seq
                        if traj == "Degradation"
                            ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])
                        else
                            ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.05, 0.05, 0.49, 0.5])
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
param = Get_classical_param()
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
                    param["alpha_0"] = alpha_e
                    for traj in trajec_seq
                        if traj == "Degradation"
                            ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])
                        else
                            ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.05, 0.05, 0.49, 0.5])
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


param = Get_classical_param()
size_landscape = 100
trajec_seq = ["Degradation"]
count = 1
scale_competition = ["global"]
disp_seq = collect(range(0, stop=1, length=12))[[2 11]]

for disp in disp_seq
    param["S"] = 0
    param["alpha_0"] = 0.2
    param["delta"] = disp
    ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])


    d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=5000,
        type_competition="global", save=false, burning=15000, N_snap=40,
        name_save="")
    CSV.write("../Table/2_species/CA/Pairs" * "_delta_" * repr(disp) * ".csv", Tables.table(state), writeheader=false)

end




#endregion



#region : 4-- Clustering species dispersal (Fig 3c)


param = Get_classical_param()
size_landscape = 100
trajec_seq = ["Degradation"]
count = 1
scale_competition = ["global"]
disp_seq = collect(range(0, stop=1, length=12))[[2 11]]

for disp in disp_seq
    param["S"] = 0.73
    param["alpha_0"] = 0.2
    param["delta"] = disp
    ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])


    d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=5000,
        type_competition="global", save=false, burning=15000, N_snap=40,
        name_save="")
    CSV.write("../Table/2_species/CA/Landscape_clustering" * "_delta_" * repr(disp) * ".csv", Tables.table(state), writeheader=false)

end




#endregion



#region : 5-- Vegetation along dispersal gradient (SI fig)



param = Get_classical_param()
size_landscape = 100
scale_comp = ["global"]
param["alpha_0"] = 0.4
disp_seq = [0 0.1 0.2 0.3 0.7 1]

for disp in disp_seq
    param["delta"] = disp
    ini = Get_initial_lattice(size_mat=size_landscape, frac=[0.4, 0.4, 0.1, 0.1])


    d, state = Run_CA_2_species(; landscape=copy(ini), param=copy(param), time=3000, type_competition=scale_comp, save=false, burning=15000, N_snap=40,
        name_save="")

    #display(Plot_dynamics(d))
    CSV.write("../Table/2_species/CA/Dispersal_gradient_delta_" * repr(disp) * ".csv", Tables.table(state), writeheader=false)


end



#endregion
