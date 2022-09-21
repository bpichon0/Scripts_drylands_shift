using StatsBase, RCall, Plots, Random, Agents, DifferentialEquations, BenchmarkTools, CSV, Tables, Distributions

#region : 0-- Functions, CA, Gillespie

function Get_params(; r, d, f, beta, m, e, emax, cintra, cinter1, cinter2, S, delta, z)
    return vec(Float64[r d f beta m e emax cintra cinter1 cinter2 S delta z])
end

function Get_classical_param()
    return Get_params(r=0.02, d=0.1, f=0.9, beta=0.8, m=0.05, e=0, emax=1,
        cintra=0.1, cinter1=0.1, cinter2=0.1, S=0, delta=0.1, z=4)
end

function Get_initial_lattice(; frac=[0.4, 0.4, 0.1, 0.1], size_mat=25)

    ini_vec = sample([1, 2, 0, -1], Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end

#params=Get_classical_param()

#sum(landscape[Base.findall(landscape.==1)])
#sum(landscape[Base.findall(landscape.==2)])



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

    r, d, f, beta, m, e, emax, cintra, cinter1, cinter2, S, delta, z = param

    colonization1 = @. beta * (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(emax * (1 - S * (1 - e)) - (cintra * rho_1 + cinter1 * rho_2))
    colonization2 = @. beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(emax * (1 - S) - (cintra * rho_2 + cinter2 * rho_1))

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




function Run_CA_2_species(; time_sim, param, landscape)

    d = Array{Float64}(undef, length(time_sim) + 1, 5) #Allocating

    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    d[1, :] = [1 rho_1 rho_2 rho_f rho_d] #the dataframe 

    @inbounds for k = Base.OneTo(length(time_sim))

        rho_1, rho_2, rho_f, rho_d, landscape = Ca_2_species(landscape=landscape, param=param)
        @views d[k+1, :] = [k + 1 rho_1 rho_2 rho_f rho_d]

    end

    return d, landscape

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

    d2 = zeros(length(time), 5) #Allocating
    r, d, f, beta, m, e, emax, cintra, cinter1, cinter2, S, delta, z, tau_leap = param


    rules_change = transpose([0 0 0 1 2 -1; 1 2 -1 0 0 0])

    #Global densities
    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded

    nb_cell = size(landscape)[1]

    #Allocating 
    Rate_landscape = zeros(nb_cell, nb_cell, 6)


    for t in time


        @rput landscape
        R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        R"neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
        R"neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

        @rget neigh_1
        @rget neigh_2
        @rget neigh_f
        @rget neigh_d

        Rate_landscape[:, :, 1] .= @. beta * (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(emax * (1 - S * (1 - e)) - (cintra * rho_1 + cinter1 * rho_2))
        Rate_landscape[:, :, 2] .= @. beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(emax * (1 - S) - (cintra * rho_2 + cinter2 * rho_1))

        # calculate regeneration, degradation & mortality rate
        Rate_landscape[:, :, 3] .= d
        Rate_landscape[:, :, 4] .= m
        Rate_landscape[:, :, 5] .= m
        Rate_landscape[:, :, 6] .= @.(r + f * neigh_1 / z)

        #calculate propensity

        propensity = [sum(Rate_landscape[:, :, k][findall(landscape .== rules_change[k, 1])]) for k in 1:size(Rate_landscape)[3]]

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


    return landscape, d2


end
#endregion


function Plot_dynamics(d)

    plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
    plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
    plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
    plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")


end

#region : Step 1-- Clustering coefficient



param = Get_classical_param();
push!(param, 0.05); #tau leap
size_landscape = 25
landscape = Get_initial_lattice(size_mat=size_landscape)
time = range(1, 1000, step=1)

d, state = Run_CA_2_species(time_sim=range(1, 1000, step=1), param=copy(param), landscape=copy(ini_land))
plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")

@time @inbounds for k in 1:10
    d, state = Run_CA_2_species(time_sim=range(1, 1000, step=1), param=param, landscape=copy(ini_land))
end



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




#Making the loop along stress gradient
c_seq = 0.3
S_seq = collect(range(0, stop=1, length=50))
param = Get_classical_param()
size_landscape = 25
ini = Get_initial_lattice(size_mat=size_landscape)
replicate = 5

d2 = zeros(length(S_seq) * length(c_seq) * replicate, 6) #Allocating
count = 1

for stress in S_seq
    param[11] = stress
    for c1 in c_seq
        for rep in 1:replicate
            param[9] = c1
            d, state = Run_CA_2_species(time_sim=range(1, 5000, step=1), param=copy(param), landscape=copy(ini))

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


CSV.write("../Table/2_species/Clustering_species_stress_gradient.csv", Tables.table(d2), writeheader=false)



#endregion

#region : Step 2-- Gillespie vs normal

p = Get_classical_param();
push!(p, 0.5); #tau leap
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
