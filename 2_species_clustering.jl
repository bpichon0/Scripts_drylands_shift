using StatsBase, RCall, Plots, Random, Agents, DifferentialEquations, BenchmarkTools, CSV, Tables, Distributions

#region : Step 0-- Functions, CA, Gillespie

# function Get_params(; r, d, f, beta, m, e, emax, cg, alpha11, alpha12, alpha21, alpha22, S, delta, z)
#     return vec(Float64[r d f beta m e emax cg alpha11 alpha12 alpha21 alpha22 S delta z])
# end

# function Get_classical_param()
#     return Get_params(r=0.05, d=0.1, f=0.9, beta=0.8, m=0.1, e=0, emax=1.2,
#         cg=0.1, alpha11=0.1, alpha12=0.1, alpha21=0.1, alpha22=0.1, S=0, delta=0.1, z=4)
# end

function Get_classical_param()

    r = 0.05
    d = 0.1
    f = 0.9
    beta = 0.8
    m = 0.1
    e = 0
    emax = 1.2
    cg = 0.1
    alpha11 = 0.1
    alpha12 = 0.1
    alpha21 = 0.1
    alpha22 = 0.1
    S = 0
    delta = 0.1
    z = 4
    tau_leap = 0.05

    return Dict{String,Any}(
        "r" => r,
        "d" => d,
        "f" => f,
        "beta" => beta,
        "m" => m,
        "e" => e,
        "emax" => emax,
        "cg" => cg,
        "alpha11" => alpha11,
        "alpha12" => alpha12,
        "alpha21" => alpha21,
        "alpha22" => alpha22,
        "S" => S,
        "delta" => delta,
        "z" => z,
        "tau_leap" => tau_leap
    )
end

function Get_initial_lattice(; frac=[0.4, 0.4, 0.1, 0.1], size_mat=25)

    ini_vec = sample([1, 2, 0, -1], Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end


function Ca_2_species(; landscape, param)


    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded


    # Neighbors :

    #using simcol package from R 
    @rput landscape
    R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

    @rget neigh_1
    @rget neigh_2
    @rget neigh_f
    @rget neigh_d

    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    emax = param["emax"]
    cg = param["cg"]
    alpha11 = param["alpha11"]
    alpha12 = param["alpha12"]
    alpha21 = param["alpha21"]
    alpha22 = param["alpha22"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]


    colonization1 = @. beta * (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(emax * (1 - S * (1 - e)) - (cg * (rho_1 + rho_2) + (alpha11 * (neigh_1 / z) + alpha21 * (neigh_2 / z))))
    colonization2 = @. beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * (neigh_2 / z) + alpha12 * (neigh_1 / z))))

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


function Run_CA_2_species(; time_sim, param, landscape, save, name_save, burning, N_snap)

    d = Array{Float64}(undef, length(time_sim) + 1, 5) #Allocating

    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    d[1, :] = [1 rho_1 rho_2 rho_f rho_d] #the dataframe 

    n_save = 1
    @inbounds for k = Base.OneTo(length(time_sim))

        rho_1, rho_2, rho_f, rho_d, landscape = Ca_2_species(landscape=landscape, param=param)
        @views d[k+1, :] = [k + 1 rho_1 rho_2 rho_f rho_d]

        if save && k > burning && k % ((length(time_sim) - burning) / N_snap) == 0
            CSV.write(name_save * "_nsave_" * repr(n_save) * ".csv", Tables.table(landscape), writeheader=false)
            n_save += 1

        end
    end


    return landscape, d

end



"Function to map the neighbors in the lattice"
function Spatial_grid(size_landscape)

    R"position2rowcol = function(position, size = 4){
      row = ceiling(position/size)
      column = (position-1)%%size + 1
      return(c(row, column))
    }"

    R"
    grid.distance = function(pos1, pos2, size = 4,border = T){
      pos1.rowcol = position2rowcol(pos1, size = size)
      pos2.rowcol = position2rowcol(pos2, size = size)
      
      dist.vector = abs(pos1.rowcol - pos2.rowcol)
      if (border == F & dist.vector[1] == size - 1){
        dist.vector[1] = 1
      }
      if (border == F & dist.vector[2] == size - 1){
        dist.vector[2] = 1
      }
      return(sum(dist.vector))
    }
    "
    R"disp_matrix_grid_4 = function(size = 4, border = F){
      npatches = size*size
      neigh = matrix(rep(0, npatches**2), nrow = npatches)
      for (k in 1:npatches){
        for (l in 1:npatches){
          if (grid.distance(k, l, size = size, border = border) == 1){
            neigh[k, l] = 1
            neigh[l, k] = 1
          }
        }
      }
      return(neigh)
    }"

    @rput size_landscape
    R"A=disp_matrix_grid_4(size_landscape)"
    @rget A
    return (A)
end


"Function for faster computing of CA using Gillespie Tau leeping method"
function Gillespie_tau_leeping(; landscape, param, time)

    d2 = zeros(time, 5) #Allocating
    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    emax = param["emax"]
    cg = param["cg"]
    alpha11 = param["alpha11"]
    alpha12 = param["alpha12"]
    alpha21 = param["alpha21"]
    alpha22 = param["alpha22"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]
    tau_leap = param["tau_leap"]


    rules_change = transpose([0 0 0 1 2 -1; 1 2 -1 0 0 0])

    #Global densities
    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    nb_cell = size(landscape)[1]

    #Allocating 
    Rate_landscape = zeros(nb_cell, nb_cell, 6)


    for t in 1:time


        @rput landscape
        R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

        @rget neigh_1
        @rget neigh_2

        Rate_landscape[:, :, 1] .= @. beta * (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(emax * (1 - S * (1 - e)) - (cg * (rho_1 + rho_2) + (alpha11 * (neigh_1 / z) + alpha21 * (neigh_2 / z))))
        Rate_landscape[:, :, 1] .= Rate_landscape[:, :, 1] .* (landscape .== 0)

        Rate_landscape[:, :, 2] .= @. beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * (neigh_2 / z) + alpha12 * (neigh_1 / z))))
        Rate_landscape[:, :, 2] .= Rate_landscape[:, :, 2] .* (landscape .== 0)

        # calculate regeneration, degradation & mortality rate
        Rate_landscape[:, :, 3] .= d .* (landscape .== 0)
        Rate_landscape[:, :, 4] .= m .* (landscape .== 1)
        Rate_landscape[:, :, 5] .= m .* (landscape .== 2)
        Rate_landscape[:, :, 6] .= @.(r + f * neigh_1 / z) .* (landscape .== -1)

        #calculate propensity

        propensity = [sum(Rate_landscape[:, :, k]) for k in 1:size(Rate_landscape)[3]]

        nb_events = map(x -> rand(Poisson(x)), propensity * tau_leap)

        for event in 1:length(nb_events) #for each type of events
            patches = findall(landscape .== rules_change[event, 1])

            if nb_events[event] != 0 && length(patches) > nb_events[event]
                landscape[wsample(patches, Rate_landscape[patches, event], nb_events[event])] .= rules_change[event, 2]
            end

        end


        Rate_landscape = zeros(nb_cell, nb_cell, 6)

        rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
        rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
        rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
        rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded


        @views d2[t, :] = [t rho_1 rho_2 rho_f rho_d]


    end


    return d2, landscape


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
#region : Step 1-- Clustering along competitive gradient



param = Get_classical_param();
push!(param, 0.05); #tau leap
size_landscape = 25
landscape = Get_initial_lattice(size_mat=size_landscape)
time = range(1, 1000, step=1)

d, state = Run_CA_2_species(time_sim=range(1, 5000, step=1), param=copy(param), landscape=copy(landscape))
plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")




#Making the loop along interspecific competition gradient
c_seq = collect(range(0.1, stop=0.5, length=50))
S_seq = [0, 0.25, 0.5]
param = Get_classical_param()
size_landscape = 25
ini = Get_initial_lattice(size_mat=size_landscape)
replicate = 1

d2 = zeros(length(S_seq) * length(c_seq) * replicate, 6) #Allocating
count = 1

for stress in S_seq
    param[11] = stress
    for c1 in c_seq
        for rep in 1:replicate
            param[9] = c1
            d, state = Run_CA_2_species(time_sim=range(1, 10000, step=1), param=copy(param), landscape=copy(ini))

            #compute neighbors 
            @rput state
            R"neigh_1= simecol::neighbors(x =state,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
            R"neigh_2= simecol::neighbors(x =state,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

            @rget neigh_1
            @rget neigh_2

            # clustering the the two species
            rho_1 = length(findall((state .== 1))) / length(state)

            if rho_1 > 0
                q12 = mean(neigh_2[findall((state .== 1))] / 4) #average number of species 2 around species 1 
            else
                q12 = 0
            end
            c12 = q12 / rho_1

            # clustering total vegetation
            rho_p = (length(findall((state .== 2))) + length(findall((state .== 1)))) / length(state)


            if length(findall((state .== 2))) > 0 & length(findall((state .== 1))) > 0
                qpp = mean(neigh_2[findall((state .== 1))] + neigh_1[findall((state .== 1))] + neigh_2[findall((state .== 2))]) / 4 #vegetation value
            else
                qpp = 0
            end
            cpp = qpp / rho_p

            @views d2[count, :] = [rep q12 c12 cpp c1 stress]

            count += 1
        end
    end
end


CSV.write("../Table/2_species/Clustering_species_interspe_compet_gradient.csv", Tables.table(d2), writeheader=false)




#endregion

#
#region : Step 2-- Clustering along stress gradient 

# My intuition is that interspecific clusters of species is observed when (i) stress_tol -> competitive is low but competitive -> stress_tol high
# Facilitation must be high, stress_tol intraspecific competition is high so that we insure that only competitive can recruit near adult facilitators
# Low intraspecific competition between species with exploitative strategies


#Making the loop along stress gradient
S_seq = collect(range(0, stop=1, length=20))
param = Get_classical_param()
param["alpha11"] = 0.2
param["alpha12"] = 0.01
param["alpha21"] = 0.2
param["alpha22"] = 0.01
size_landscape = 25
ini = Get_initial_lattice(size_mat=size_landscape)
replicate = 1
max_time = 5000

d2 = zeros(length(S_seq) * replicate, 12) #Allocating
count = 1

for stress in S_seq
    param["S"] = stress
    for rep in 1:replicate
        state, d = Run_CA_2_species(time_sim=range(1, max_time, step=1), param=copy(param), landscape=copy(ini))

        #compute neighbors 
        @rput state
        R"neigh_1= simecol::neighbors(x =state,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        R"neigh_2= simecol::neighbors(x =state,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

        @rget neigh_1
        @rget neigh_2

        # clustering of the two species : q12, q21
        rho_1 = length(findall((state .== 1))) / length(state)

        if rho_1 > 0
            q12 = mean(neigh_2[findall((state .== 1))] / 4) #average number of species 2 around species 1 
            c12 = q12 / rho_1

            q11 = mean(neigh_1[findall((state .== 1))] / 4)
            c11 = q11 / rho_1
        else
            q12 = 0
            c12 = 0
            q11 = 0
            c11 = 0
        end

        rho_2 = length(findall((state .== 2))) / length(state)

        if rho_2 > 0
            q21 = mean(neigh_1[findall((state .== 2))] / 4) #average number of species 1 around species 2 
            c21 = q21 / rho_2

            q22 = mean(neigh_2[findall((state .== 2))] / 4)
            c22 = q22 / rho_2
        else
            q21 = 0
            c21 = 0
            q22 = 0
            c22 = 0
        end


        # clustering total vegetation
        rho_p = (length(findall((state .== 2))) + length(findall((state .== 1)))) / length(state)


        if length(findall((state .== 2))) > 0 & length(findall((state .== 1))) > 0 #we can compute q++
            qpp = mean(neigh_2[findall((state .== 1))] + neigh_1[findall((state .== 1))] + neigh_2[findall((state .== 2))]) / 4 #vegetation value
            cpp = qpp / rho_p
        else
            qpp = 0
            cpp = 0
        end


        @views d2[count, :] = [rep q12 q21 q11 q22 c12 c21 c11 c22 qpp cpp stress]

        count += 1
    end
end


CSV.write("../Table/2_species/Clustering_species_stress_gradient.csv", Tables.table(d2), writeheader=false)

#endregion

#
#region : Step 3-- Gillespie vs normal

p = Get_classical_param();
size_landscape = 25
init = Get_initial_lattice(size_mat=size_landscape)
time = range(1, 5000, step=1)


@time for k in 1:10
    Gillespie_tau_leeping(; landscape=init, param=p, time)
end

@time land, dynamics = Gillespie_tau_leeping(; landscape=init, param=p, time)
plot(dynamics[:, 1], dynamics[:, 2], seriescolor=:lightgreen, label="stress_tol")
plot!(dynamics[:, 1], dynamics[:, 3], seriescolor=:blue, label="competitive")
plot!(dynamics[:, 1], dynamics[:, 4], seriescolor=:orange, label="fertile")
plot!(dynamics[:, 1], dynamics[:, 5], seriescolor=:grey, label="degraded")



@time for k in 1:10
    Run_CA_2_species(time_sim=time, param=p, landscape=init)
end

@time d, state = Run_CA_2_species(time_sim=time, param=p, landscape=init)
plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")


#endregion

#
#region : Step 4-- Patch size distribution : different competition scenario  

S_seq = collect(range(0, stop=1, length=6))
param = Get_classical_param()
param["alpha11"] = 0.3
param["alpha12"] = 0.01
param["alpha21"] = 0.2
param["alpha22"] = 0.01
param["S"] = 0.6
param["tau_leap"] = 0.1
size_landscape = 100
ini = Get_initial_lattice(size_mat=size_landscape)
max_time = 2000
rep = 3

for stress in S_seq
    for nrep in 1:3
        param["S"] = stress
        #state, d = Gillespie_tau_leeping(time=range(1, max_time, step=1), param=copy(param), landscape=copy(ini))
        state, d = Run_CA_2_species(time_sim=range(1, max_time, step=1), param=copy(param), landscape=copy(ini))
        CSV.write("../Table/2_species/Patch_size/Example_sim/Landscape_size_100_stress_" * repr(stress) * "_" * repr(nrep) * ".csv", Tables.table(state), writeheader=false)
    end
end



S_seq = collect(range(0, stop=1, length=6))
param = Get_classical_param()
param["S"] = 0.6
param["tau_leap"] = 0.1
size_landscape = 100
ini = Get_initial_lattice(size_mat=size_landscape)
max_time = 2000
rep = 3
for rela_c in [1, 2.5, 4]
    for scena in (1:2)

        if scena == 1 #low inter
            name_scena = "low_inter"
            param["cg"] = 0.01
            param["alpha12"] = 0.05
            param["alpha21"] = param["alpha12"] * rela_c
            param["alpha11"] = 0.2
            param["alpha22"] = 0.2
        else #inter=intra
            name_scena = "intra=inter"
            param["cg"] = 0.01
            param["alpha12"] = 0.1
            param["alpha21"] = param["alpha12"] * rela_c
            param["alpha11"] = 0.1
            param["alpha22"] = 0.1
        end

        for stress in S_seq
            for nrep in 1:3
                param["S"] = stress
                #state, d = Gillespie_tau_leeping(time=range(1, max_time, step=1), param=copy(param), landscape=copy(ini))
                state, d = Run_CA_2_species(time_sim=range(1, max_time, step=1), param=copy(param), landscape=copy(ini))
                CSV.write("../Table/2_species/Patch_size/Competition_regime/Landscape_size_100_stress_" * repr(stress) * "_" * name_scena * "_" * repr(nrep) * "_relat_compet_" * repr(rela_c) * ".csv", Tables.table(state), writeheader=false)
            end
        end
    end
end

#endregion

#
#region : Step 5-- Fitting power laws for patch size distribution 
# S varies for 3-4 competition strength

max_time = 10000
n_snapshot = 25
burning_phase = 3000
S_seq = collect(range(0, stop=0.8, length=30))
param = Get_classical_param()

size_landscape = 100
ini = Get_initial_lattice(size_mat=size_landscape)

name_scena = ["low_inter", "inter=intra"]

for scena in (1)
    for stress in S_seq
        param["S"] = stress
        for rela_c in [1, 2.5, 4]

            if scena == 1 #low inter
                name_scena = "low_inter"
                param["cg"] = 0.01
                param["alpha12"] = 0.05
                param["alpha21"] = param["alpha12"] * rela_c
                param["alpha11"] = 0.2
                param["alpha22"] = 0.2
            else #inter=intra
                name_scena = "intra=inter"
                param["cg"] = 0.01
                param["alpha12"] = 0.1
                param["alpha21"] = param["alpha12"] * rela_c
                param["alpha11"] = 0.1
                param["alpha22"] = 0.1
            end


            param["alpha21"] = param["alpha12"] * rela_c
            name_landscape = "../Table/2_species/Patch_size/Big_sim/Landscape_" * name_scena[scena] * "_S_" * repr(round(stress, digits=2)) * "_a21_" * repr(rela_c)
            state, d = Run_CA_2_species(time_sim=range(1, max_time, step=1), param=copy(param),
                landscape=copy(ini), save=true, name_save=name_landscape, burning=burning_phase, N_snap=n_snapshot)

        end
    end
end










#endregion

#
#region : Step 6-- Patagonian steppe structure


function Ca_2_species2(; landscape, param)


    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded


    # Neighbors :

    #using simcol package from R 
    @rput landscape
    R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

    @rget neigh_1
    @rget neigh_2
    @rget neigh_f
    @rget neigh_d

    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    emax = param["emax"]
    cg = param["cg"]
    alpha11 = param["alpha11"]
    alpha12 = param["alpha12"]
    alpha21 = param["alpha21"]
    alpha22 = param["alpha22"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]


    colonization1 = @. beta * ((1 - delta) * rho_1 + delta * neigh_1 / z) * @.(emax * (1 - S * (1 - e)) - (cg * (rho_1 + rho_2) + (alpha11 * (neigh_1 / z) + alpha21 * (neigh_2 / z))))
    colonization2 = @. beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * (neigh_2 / z) + alpha12 * (neigh_1 / z))))

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

function Run_CA_2_species2(; time_sim, param, landscape)

    d = Array{Float64}(undef, length(time_sim) + 1, 5) #Allocating

    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    d[1, :] = [1 rho_1 rho_2 rho_f rho_d] #the dataframe 

    @inbounds for k = Base.OneTo(length(time_sim))

        rho_1, rho_2, rho_f, rho_d, landscape = Ca_2_species2(landscape=landscape, param=param)
        @views d[k+1, :] = [k + 1 rho_1 rho_2 rho_f rho_d]

    end

    return landscape, d

end

#two examples for influence of dispersal range
param["alpha11"] = 0.3
param["alpha12"] = 0.01
param["alpha21"] = 0.2
param["alpha22"] = 0.01
param["cg"] = 0.01
param["S"] = 0.2
param["tau_leap"] = 0.1
size_landscape = 100
ini = Get_initial_lattice(size_mat=size_landscape)
max_time = 1000
rep = 3

#similar dispersal
state, d = Run_CA_2_species(time_sim=range(1, max_time, step=1), param=copy(param), landscape=copy(ini),
    save=false, burning=10, name_save="", N_snap=0)
Plot_landscape(state)
savefig("../Figures/2_species/CA/Patagonian_type_same_disp.png")


#different dispersal
state, d = Run_CA_2_species2(time_sim=range(1, max_time, step=1), param=copy(param), landscape=copy(ini))
Plot_landscape(state)
Plot_dynamics(d)
savefig("../Figures/2_species/CA/Patagonian_type_different_disp.png")



@time state, d = Gillespie_tau_leeping(landscape=copy(ini), param=copy(param), time=range(1, max_time, step=1))
Plot_dynamics(d)

S_seq = collect(range(0, stop=8, length=25))
a_seq = collect(range(0.1, stop=1.5, length=25))
param = Get_classical_param()

#endregion
