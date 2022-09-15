using StatsBase, RCall, Plots, Random, Agents, DifferentialEquations, BenchmarkTools, CSV, Tables


function Get_params(; delta, beta, emax, cintra, cinter1, cinter2, m, d, r, f, e, S, z)
    return vec(Float64[delta beta emax cintra cinter1 cinter2 m d r f e S z])
end

function Get_classical_param()
    return Get_params(delta=0.1, beta=0.8, emax=1, cintra=0.1, cinter1=0.1, cinter2=0.1,
        m=0.05, d=0.1, r=0.02, f=0.9, e=0.1, S=0, z=4)
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

    delta, beta, emax, cintra, cinter1, cinter2, m, d, r, f, e, S, z = param

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
    landscape_update[findall((landscape .== -1) .& (rnum .> colonization1 .+ colonization2) .& (rnum .<= (colonization1 .+ colonization2 .+ degradation)))] .= 3

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







param = Get_classical_param()
size_landscape = 25
landscape = Get_initial_lattice(size_mat=size_landscape)

d, state = Run_CA_2_species(time_sim=range(1, 1000, step=1), param=copy(param), landscape=copy(ini_land))
plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")

@time @inbounds for k in 1:10
    d, state = Run_CA_2_species(time_sim=range(1, 1000, step=1), param=param, landscape=copy(ini_land))
end



#Making the loop along interspecific competition gradient
c_seq = collect(range(0.1, stop=0.5, length=30))
S_seq = [0, 0.25, 0.5]
param = Get_classical_param()
size_landscape = 25
ini = Get_initial_lattice(size_mat=size_landscape)
replicate = 5

d2 = Array{Float64}(undef, length(S_seq) * length(c_seq) * replicate, 6) #Allocating
count = 1

for stress in S_seq
    param[12] = stress
    for c1 in c_seq
        for rep in 1:replicate
            param[5] = c1
            d, state = Run_CA_2_species(time_sim=range(1, 2000, step=1), param=copy(param), landscape=copy(ini))

            #compute neighbors 
            @rput state
            R"neigh_1= simecol::neighbors(x =state,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
            R"neigh_2= simecol::neighbors(x =state,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

            @rget neigh_1
            @rget neigh_2

            # clustering the the two species
            rho_1 = length(findall((landscape .== 2))) / length(landscape)
            q12 = mean(neigh_2[findall((landscape .== 1))] / 4) #average number of species 2 around species 1 
            c12 = q12 / rho_1

            # clustering total vegetation
            rho_p = (length(findall((landscape .== 2))) + length(findall((landscape .== 1)))) / length(landscape)
            qpp = mean(neigh_2[findall((landscape .== 1))] / 4) + mean(neigh_1[findall((landscape .== 1))] / 4) + mean(neigh_2[findall((landscape .== 2))] / 4) #vegetation value
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

d2 = Array{Float64}(undef, length(S_seq) * length(c_seq) * replicate, 6) #Allocating
count = 1

for stress in S_seq
    param[12] = stress
    for c1 in c_seq
        for rep in 1:replicate
            param[5] = c1
            d, state = Run_CA_2_species(time_sim=range(1, 2000, step=1), param=copy(param), landscape=copy(ini))

            #compute neighbors 
            @rput state
            R"neigh_1= simecol::neighbors(x =state,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
            R"neigh_2= simecol::neighbors(x =state,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"

            @rget neigh_1
            @rget neigh_2

            # clustering the the two species
            rho_1 = length(findall((landscape .== 2))) / length(landscape)
            q12 = mean(neigh_2[findall((landscape .== 1))] / 4) #average number of species 2 around species 1 
            c12 = q12 / rho_1

            # clustering total vegetation
            rho_p = (length(findall((landscape .== 2))) + length(findall((landscape .== 1)))) / length(landscape)
            qpp = mean(neigh_2[findall((landscape .== 1))] / 4) + mean(neigh_1[findall((landscape .== 1))] / 4) + mean(neigh_2[findall((landscape .== 2))] / 4) #vegetation value
            cpp = qpp / rho_p

            @views d2[count, :] = [rep q12 c12 cpp c1 stress]

            count += 1
        end
    end
end


CSV.write("../Table/2_species/Clustering_species_stress_gradient.csv", Tables.table(d2), writeheader=false)
