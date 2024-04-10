using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings,
    Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames, ColorSchemes



## 2 species functions



function Get_classical_param_2species()

    r = 0.01
    d = 0.025
    f = 0.9
    beta = 1
    m = 0.15
    e = 0.1
    cg = 0
    alpha_e = 0.3
    cintra = 0.3
    S = 0
    delta = 0.1
    z = 4
    tau_leap = 0.5

    return Dict{String,Any}(
        "r" => r,
        "d" => d,
        "f" => f,
        "beta" => beta,
        "m" => m,
        "e" => e,
        "cg" => cg,
        "alpha_e" => alpha_e,
        "cintra" => cintra,
        "S" => S,
        "delta" => delta,
        "z" => z,
        "tau_leap" => tau_leap
    )
end

function Get_initial_lattice_2species(; frac=[0.4, 0.4, 0.1, 0.1], size_mat=25)
    ini_vec = sample([1, 2, 0, -1], Weights(frac), size_mat * size_mat)
    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
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
    alpha_e = param["alpha_e"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]

    colonization1 = @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - (cintra * rho_1 + alpha_e * (1 + exp(-1)) * rho_2))
    colonization2 = @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - (cintra * rho_2 + alpha_e * rho_1))
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
    cintra = param["cintra"]
    alpha_e = param["alpha_e"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]

    colonization1 = @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - (cintra * neigh_1 / z + alpha_e * (1 + exp(-1)) * neigh_2 / z))
    colonization2 = @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - (cintra * neigh_2 / z + alpha_e * neigh_1 / z))
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




function Ca_2_species_medium(; landscape, param)


    rho_1 = length(findall((landscape .== 1))) / length(landscape) #fraction stress_tol
    rho_2 = length(findall((landscape .== 2))) / length(landscape) #fraction competitive
    rho_f = length(findall((landscape .== 0))) / length(landscape) #fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_f # fraction degraded


    # Neighbors :

    #using simcol package from R 
    @rput landscape
    R"neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_1_comp= simecol::neighbors(x =landscape,state = 1, wdist =  matrix(c(rep(1,12),0,rep(1,12)),5,5),bounds = 1)"
    R"neigh_2_comp= simecol::neighbors(x =landscape,state = 2, wdist =  matrix(c(rep(1,12),0,rep(1,12)),5,5),bounds = 1)"


    @rget neigh_1
    @rget neigh_2
    @rget neigh_1_comp
    @rget neigh_2_comp

    r = param["r"]
    d = param["d"]
    f = param["f"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    cintra = param["cintra"]
    alpha_e = param["alpha_e"]
    S = param["S"]
    delta = param["delta"]
    z = param["z"]

    colonization1 = @. (delta * rho_1 + (1 - delta) * neigh_1 / z) * @.(beta * (1 - S * (1 - e)) - (cintra * neigh_1_comp / 24 + alpha_e * (1 + exp(-1)) * neigh_2_comp / 24))
    colonization2 = @. (delta * rho_2 + (1 - delta) * neigh_2 / z) * @.(beta * (1 - S) - (cintra * neigh_2_comp / 24 + alpha_e * neigh_1_comp / 24))
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
    elseif type_competition=="medium"#medium competition

        @inbounds for k in 1:time

            rho_1, rho_2, rho_f, rho_d, landscape = Ca_2_species_medium(landscape=copy(landscape), param=param)
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



function Plot_dynamics_2species(d)

    plot(d[:, 1], d[:, 2], seriescolor=:lightgreen, label="stress_tol")
    plot!(d[:, 1], d[:, 3], seriescolor=:blue, label="competitive")
    plot!(d[:, 1], d[:, 4], seriescolor=:orange, label="fertile")
    plot!(d[:, 1], d[:, 5], seriescolor=:grey, label="degraded")


end

function Plot_landscape_2species(landscape)
    colGRAD = cgrad([colorant"#696969", colorant"#D8CC7B", colorant"blue", colorant"#ACD87B"])
    heatmap(landscape, yflip=true, fill=true, c=colGRAD)
end






## N species functions
#Parameters functions 

function Get_classical_param_Nspecies(; N_species=4, alpha_e=0.1, cintra=0.2, scenario_trait="spaced", trade_off=1,type_kernel="both")

    r = 0.01
    d = 0.025
    f = 0.9
    beta = 1
    m = 0.15
    e = 0.1
    cg = 0
    S = 0
    delta = 0.1
    z = 4
    tau_leap = 0.05
    trait = Get_species_traits(Nsp=N_species, scena=scenario_trait)
    trade_off = trade_off

    Interaction_mat = Get_interaction_matrix(Nsp=N_species, alpha_e=alpha_e, cintra=cintra, trait_sp=trait, trade_off=1,type_kernel=type_kernel)

    return [N_species, r, d, f, beta, m, e, cg, S, delta, z, tau_leap, trait, Interaction_mat, trade_off]

end



function Get_classical_param_dict(; N_species=4, alpha_e=0.1, cintra=0.2, scenario_trait="spaced")

    r = 0.01
    d = 0.025
    f = 0.9
    beta = 1
    m = 0.15
    e = 0.1
    cg = 0
    S = 0
    delta = 0.1
    z = 4
    tau_leap = 0.05
    trait = Get_species_traits(Nsp=N_species, scena=scenario_trait)

    Interaction_mat = Get_interaction_matrix(Nsp=N_species, alpha_e=alpha_e, cintra=cintra, trait_sp=trait)

    return Dict{String,Any}(
        "Nsp" => N_species,
        "r" => r,
        "d" => d,
        "f" => f,
        "beta" => beta,
        "m" => m,
        "e" => e,
        "cg" => cg,
        "S" => S,
        "delta" => delta,
        "z" => z,
        "tau_leap" => tau_leap,
        "trait" => trait,
        "alpha" => Interaction_mat)
end


function Get_interaction_matrix(; Nsp, trait_sp, cintra=0.5, alpha_e, trade_off=1,type_kernel="both")

    Interaction_mat = zeros(Nsp, Nsp)
    for i in 1:Nsp
        Interaction_mat[i, i] = cintra
    end
    
    if type_kernel=="both"
      for i in 1:(Nsp-1)
          for j in (i+1):Nsp
              if i != j
                  Interaction_mat[i, j] = alpha_e * (1 + exp(-abs(trait_sp[i] - trait_sp[j])) * trait_sp[i]^trade_off) #j on i
                  Interaction_mat[j, i] = alpha_e * (1 + exp(-abs(trait_sp[i] - trait_sp[j])) * trait_sp[j]^trade_off) #i on j
              end
          end
      end
    else 
      for i in 1:(Nsp-1)
        for j in (i+1):Nsp
            if i != j
                Interaction_mat[i, j] = alpha_e * (1 + trait_sp[i]^trade_off) #j on i
                Interaction_mat[j, i] = alpha_e * (1 + trait_sp[j]^trade_off) #i on j
            end
        end
      end
    end

    return Interaction_mat
end





function Get_species_traits(; Nsp, scena)
    "Scena can either be equal to 'random' or 'spaced' or to the fraction of stress tolerant species in the community"
    if scena == "random"
        trait = rand(Uniform(0, 1), Nsp)
    elseif scena == "spaced"
        trait = collect(range(1, stop=0, length=Nsp))
    else #community composition is given
        if scena > 1 #in case of percentages
            scena /= 100
        end
        trait_stress_tol = rand(Uniform(0.7, 1), Int64(trunc(Nsp * scena)))
        trait_compet = rand(Uniform(0, 0.7), Nsp - length(trait_stress_tol))
        trait = [trait_stress_tol; trait_compet]
    end

    return trait
end


#Initial state functions 

function Get_initial_lattice_Nspecies(; param, branch="Degradation", size_mat=25, type_ini="equal")

    frac = Get_initial_state(Nsp=param["Nsp"], type=type_ini, branch=branch)
    species_vec = [k for k in 1:param["Nsp"]]
    push!(species_vec, 0)
    push!(species_vec, -1)
    ini_vec = sample(species_vec, Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end

function Get_initial_state(; Nsp, type, branch, PA=false)
    if branch == "Degradation"  #degradation branch

        if type == "equal"
            frac = [0.8 / Nsp for k in 1:Nsp]
        elseif type == "random"
            frac = rand(Nsp)
            frac = (frac / sum(frac)) * 0.8
        end
        push!(frac, 0.1)
        push!(frac, 0.1)

    else #restoration branch

        if type == "equal"
            frac = [0.01 / Nsp for k in 1:Nsp]
        elseif type == "random"
            frac = rand(Nsp)
            frac = (frac / sum(frac)) * 0.01
        end
        push!(frac, 0.50)
        push!(frac, 0.49)
    end
    if PA

        with_rho_ii = vcat(vec(frac), vec(frac .* frac))
        with_rho_i0 = vcat(vec(with_rho_ii), vec(frac[1:Nsp] .* frac[Nsp+1]))
        with_rho_0m = vcat(vec(with_rho_i0), frac[Nsp+1] * frac[Nsp+2])
        with_rho_mi = vcat(vec(with_rho_0m), vec(frac[1:Nsp] .* frac[Nsp+2]))

        ini_rho_ij = zeros(Nsp, Nsp)
        for i in 1:(Nsp-1)
            for j in (i+1):Nsp
                ini_rho_ij[i, j] = frac[i] * frac[j]
            end
        end

        frac = copy(vcat(vec(with_rho_mi), vec(transpose(ini_rho_ij)[tril!(trues(size(ini_rho_ij)), -1)])))
    end
    return frac
end

# MF functions 
function MF_N_species_dict(du, u, p, t)
    rho_0 = u[p["Nsp"]+1]
    rho_d = u[p["Nsp"]+2]

    for k in 1:p["Nsp"]
        du[k] = rho_0 * (u[k] * (p["beta"] * (1 - p["S"] * (1 - p["trait"][k] * p["e"])) - (sum([p["alpha"][k, i] * u[i] for i in 1:p["Nsp"]])))) - u[k] * p["m"]
    end
    du[p["Nsp"]+1] = -p["d"] * rho_0 + rho_d * (p["r"] + p["f"] * (sum([p["trait"][i] * u[i] for i in 1:p["Nsp"]])))

    for k in 1:p["Nsp"]
        du[p["Nsp"]+1] = du[p["Nsp"]+1] - rho_0 * (u[k] * (p["beta"] * (1 - p["S"] * (1 - p["trait"][k] * p["e"])) - (sum([p["alpha"][k, i] * u[i] for i in 1:p["Nsp"]])))) + u[k] * p["m"]
    end

    du[p["Nsp"]+2] = p["d"] * rho_0 - rho_d * (p["r"] + p["f"] * (sum([p["trait"][i] * u[i] for i in 1:p["Nsp"]])))

end


# MF functions 
function MF_N_species(du, u, p, t)
    Nsp = p[1]
    r = p[2]
    d = p[3]
    f = p[4]
    beta = p[5]
    m = p[6]
    e = p[7]
    S = p[9]
    trait = p[13]
    alpha = p[14]
    trade_off = p[16]

    rho_0 = u[Nsp+1]
    rho_d = u[Nsp+2]

    for k in 1:Nsp
        du[k] = rho_0 * (u[k] * (beta * (1 - S * (1 - trait[k]^trade_off * e)) - (sum([alpha[k, i] * u[i] for i in 1:Nsp])))) - u[k] * m
    end
    du[Nsp+1] = -d * rho_0 + rho_d * (r + f * (sum([trait[i]^trade_off * u[i] for i in 1:Nsp])))

    for k in 1:Nsp
        du[Nsp+1] = du[Nsp+1] - rho_0 * (u[k] * (beta * (1 - S * (1 - trait[k]^trade_off * e)) - (sum([alpha[k, i] * u[i] for i in 1:Nsp])))) + u[k] * m
    end

    du[Nsp+2] = d * rho_0 - rho_d * (r + f * (sum([trait[i]^trade_off * u[i] for i in 1:Nsp])))

end


function Reorder_dynamics(sol)
    sol_dyn = reduce(hcat, sol.u)'     #transforming solutions in clean matrix
    sol_dyn = cat(sol_dyn, sol.t, dims=2)
    return sol_dyn
end

function PA_N_species(du, u, p, t)

    Nsp = p[1]
    r = p[2]
    d = p[3]
    f = p[4]
    beta = p[5]
    m = p[6]
    e = p[7]
    S = p[9]
    delta = p[10]
    z = p[11]
    trait = p[13]
    alpha = p[14]


    rho_i = u[1:(Nsp)]
    rho_0 = u[Nsp+1]
    rho_m = u[Nsp+2]
    rho_ii = u[(Nsp+3):(2*Nsp+4)]
    rho_i0 = u[(2*Nsp+5):(3*Nsp+5)] #The last one being rho_i0
    rho_im = u[(3*Nsp+6):(4*Nsp+5)]
    rho_ij = u[(4*Nsp+6):Int(((Nsp^2 + 7 * Nsp) / 2) + 5)]


    #rho_i -> OK
    for i in 1:Nsp
        du[i] = rho_0 * (delta * rho_i[i] + (1 - delta) * (rho_i0[i] / rho_0)) *
                (beta * (1 - S * (1 - trait[i] * e)) -
                 (sum([rho_i[k] * alpha[i, k] for k in 1:(Nsp)]))) -  #recruitment of plant i
                rho_i[i] * m #mortality

    end

    #rho_0 -> OK
    du[Nsp+1] = rho_m * (r + f * (sum([trait[k] * (rho_im[k] / rho_m) for k in 1:Nsp]))) - #restoration
                d * rho_0 - #degradation
                sum([du[k] for k in 1:Nsp]) #sum of plant recruitment

    #rho_m -> OK
    du[Nsp+2] = d * rho_0 - #degradation
                rho_m * (r + f * (sum([trait[k] * (rho_im[k] / rho_m) for k in 1:Nsp]))) #restoration


    #rho_ii -> OK
    for i in (Nsp+3):(2*Nsp+2)
        du[i] = 2 * rho_i0[i-(Nsp+2)] * (delta * rho_i[i-(Nsp+2)] + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_i0[i-(Nsp+2)] / rho_0)) *
                (beta * (1 - S * (1 - trait[i-(Nsp+2)] * e)) -
                 (sum([rho_i[k] * alpha[i-(Nsp+2), k] for k in 1:(Nsp)]))) - #recruitment
                2 * m * rho_ii[i-(Nsp+2)] #mortality
    end

    #rho_00 -> OK
    du[2*Nsp+3] = 2 * sum([rho_i0[k] * m for k in 1:Nsp]) + #plant mortality
                  2 * rho_i0[Nsp+1] * (r + ((z - 1) / z) * f * (sum([trait[k] * (rho_im[k] / rho_m) for k in 1:Nsp]))) - #restoration from rho_0m
                  2 * rho_ii[Nsp+1] * d - #degradation
                  2 * sum([rho_ii[Nsp+1] * ((delta * rho_i[k] + (1 - delta) * ((z - 1) / z) * (rho_i0[k] / rho_0)) *
                                            (beta * (1 - S * (1 - trait[k] * e)) -
                                             (sum([rho_i[x] * alpha[k, x] for x in 1:(Nsp)])))) for k in 1:Nsp]) #recruitment, we sum over all species recruitment possibilities


    #rho_mm -> OK
    du[2*Nsp+4] = 2 * rho_i0[Nsp+1] * d -
                  2 * rho_ii[Nsp+2] * (r + ((z - 1) / z) * f * (sum([trait[k] * (rho_im[k] / rho_m) for k in 1:Nsp])))




    #rho_m0 -> need rho_mi so the value is fixed after
    du[3*Nsp+5] = 1


    #rho_mi -> OK
    for i in (3*Nsp+6):(4*Nsp+5)

        du[i] = rho_i0[i-(3*Nsp+5)] * d + #degradation
                rho_i0[Nsp+1] * (delta * rho_i[i-(3*Nsp+5)] + (1 - delta) * ((z - 1) / z) * (rho_i0[i-(3*Nsp+5)] / rho_0)) *
                (beta * (1 - S * (1 - trait[i-(3*Nsp+5)] * e)) -
                 (sum([rho_i[k] * alpha[(i-(3*Nsp+5)), k] for k in 1:(Nsp)]))) - #recruitment
                rho_im[i-(3*Nsp+5)] * m - #mortality
                rho_im[i-(3*Nsp+5)] * (r + (f * trait[(i-(3*Nsp+5))] / z) + ((z - 1) / z) * f * (sum([trait[k] * (rho_im[k] / rho_m) for k in 1:Nsp]))) #restoration
    end

    #now we fix rho_m0
    du[3*Nsp+5] = du[Nsp+2] - sum(du[(3*Nsp+6):(4*Nsp+5)]) - du[2*Nsp+4]#rho_m0 = rho_m - sum_j rho_mj (j species) - rho_mm





    #rho_ij species pairs --> OK
    #We complete the triangular sup of a matrix before turning that into a vector
    species_pairs = zeros(Nsp, Nsp)
    if (Nsp > 2)
        for i in 1:(Nsp-1)
            for j in (i+1):Nsp

                if i == 1
                    line = 0
                else
                    line = sum([Nsp - k for k in 1:(i-1)])
                end

                species_pairs[i, j] = rho_i0[i] * (delta * rho_i[j] + (1 - delta) * ((z - 1) / z) * (rho_i0[j] / rho_0)) *
                                      (beta * (1 - S * (1 - trait[j] * e)) -
                                       (sum([rho_i[k] * alpha[j, k] for k in 1:(Nsp)]))) + #recruitment j next to i
                                      rho_i0[j] * (delta * rho_i[i] + (1 - delta) * ((z - 1) / z) * (rho_i0[i] / rho_0)) *
                                      (beta * (1 - S * (1 - trait[i] * e)) -
                                       (sum([rho_i[k] * alpha[i, k] for k in 1:(Nsp)]))) - #recruitment j
                                      2 * m * rho_ij[line+abs(j - i)] #mortality. We need to get the index that coresponds to the proper rho_ij in the vector rho_ij

            end
        end

    else
        for i in 1:(Nsp-1)
            for j in (i+1):Nsp

                species_pairs[i, j] = rho_i0[i] * (delta * rho_i[j] + (1 - delta) * ((z - 1) / z) * (rho_i0[j] / rho_0)) *
                                      (beta * (1 - S * (1 - trait[j] * e)) -
                                       (sum([rho_i[k] * alpha[j, k] for k in 1:(Nsp)]))) + #recruitment j next to i
                                      rho_i0[j] * (delta * rho_i[i] + (1 - delta) * ((z - 1) / z) * (rho_i0[i] / rho_0)) *
                                      (beta * (1 - S * (1 - trait[i] * e)) -
                                       (sum([rho_i[k] * alpha[i, k] for k in 1:(Nsp)]))) - #recruitment i next to i
                                      2 * rho_ij[1] * m #mortality. We need to get the index that coresponds to the proper rho_ij in the vector rho_ij

            end
        end
    end

    species_pairs += transpose(species_pairs)

    du[(4*Nsp+6):Int(((Nsp^2 + 7 * Nsp) / 2) + 5)] = transpose(species_pairs)[tril!(trues(size(species_pairs)), -1)] #to get triangular sup matrix


    #rho_i0 = rho_i - rho_mi - rho_ii
    for i in (2*Nsp+5):(3*Nsp+4)
        du[i] = du[i-(2*Nsp+5)+1] - du[(i+Nsp+1)] - du[(i-(Nsp+3)+1)] - sum(species_pairs[i-(2*Nsp+5)+1, :])
    end

end


#CA functions
function Get_neighbors_matrix(; landscape, sp)
    @rput sp
    @rput landscape
    R"neighbors_mat = simecol::neighbors(x =landscape,state = sp, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    @rget neighbors_mat
end

function Gillespie_CA_N_species(; param, landscape, tmax, type_competition)


    r = param["r"]
    d = param["d"]
    f = param["f"]
    interaction_mat = param["alpha"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    Nsp = param["Nsp"]
    delta = param["delta"]
    S = param["S"]
    z = param["z"]
    tau_leap = param["tau_leap"]
    trait_sp = param["trait"]

    d2 = zeros(tmax, Nsp + 3) #Allocating

    rules_change = hcat(vcat([0 for x in 1:(Nsp)], [k for k in 1:Nsp], -1, 0), vcat([x for x in 1:Nsp], [0 for k in 1:(Nsp)], 0, -1))

    nb_cell = size(landscape)[1]

    rho = [length(findall((landscape .== k))) / length(landscape) for k in 1:Nsp] #vegetation

    #Allocating
    Rate_landscape = zeros(nb_cell, nb_cell, 2 * Nsp + 2) #death and recruitment rates, degradation + restoration

    if type_competition == "global"

        for dt in 1:tmax

            neigh_matrix = zeros(nb_cell, nb_cell, Nsp)
            for k in 1:Nsp
                neigh_matrix[:, :, k] = Get_neighbors_matrix(landscape=copy(landscape), sp=copy(k))
            end

            for sp in 1:Nsp
                Rate_landscape[:, :, sp] .= beta .* (delta * rho[sp] .+ (1 - delta) * neigh_matrix[:, :, sp] / z) .* #dispersal
                                            ((1 - S * (1 - trait_sp[sp] * e)) .- #recruitment
                                             ((sum([interaction_mat[sp, k] * #global competition
                                                    rho[k] for k in 1:Nsp])))) .* (landscape .== 0) #local competition

                Rate_landscape[:, :, Nsp+sp] .= m .* (landscape .== sp) #to have mortality only in species cells

            end
            # calculate regeneration, degradation & mortality rate
            Rate_landscape[:, :, 2*Nsp+1] .= (r .+ f .* sum([trait_sp[k] .* neigh_matrix[:, :, k] for k in 1:Nsp]) ./ z) .* (landscape .== -1)
            Rate_landscape[:, :, 2*Nsp+2] .= d .* (landscape .== 0)




            propensity = [sum(Rate_landscape[:, :, k]) for k in 1:size(Rate_landscape)[3]]

            nb_events = map(x -> rand(Poisson(x)), propensity * tau_leap)

            for event in 1:length(nb_events) #for each type of events : colonization, death etc...
                patches = findall(landscape .== rules_change[event, 1])

                if nb_events[event] != 0 && length(patches) > nb_events[event]
                    landscape[wsample(patches, Rate_landscape[patches, event], nb_events[event])] .= rules_change[event, 2]
                end

            end #end event loop

            rho = [length(findall((landscape .== k))) / length(landscape) for k in 1:Nsp] #vegetation
            push!(rho, length(findall((landscape .== 0))) / length(landscape)) #empty sites
            push!(rho, length(findall((landscape .== -1))) / length(landscape)) #degraded sites
            push!(rho, dt)

            @views d2[dt, :] = rho

            Rate_landscape = zeros(nb_cell, nb_cell, 2 * Nsp + 2) #death and recruitment rates, degradation + restoration

        end #end time loop

    else

        for dt in 1:tmax

            neigh_matrix = zeros(nb_cell, nb_cell, Nsp)
            for k in 1:Nsp
                neigh_matrix[:, :, k] = Get_neighbors_matrix(landscape=copy(landscape), sp=copy(k))
            end

            for sp in 1:Nsp
                Rate_landscape[:, :, sp] .= beta .* (delta * rho[sp] .+ (1 - delta) * neigh_matrix[:, :, sp] / z) .* #dispersal
                                            ((1 - S * (1 - trait_sp[sp] * e)) .- #recruitment
                                             ((sum([interaction_mat[sp, k] *
                                                    neigh_matrix[:, :, k] for k in 1:Nsp]) / z))) .* (landscape .== 0)

                Rate_landscape[:, :, Nsp+sp] .= m .* (landscape .== sp) #to have mortality only in species cells

            end
            # calculate regeneration, degradation & mortality rate
            Rate_landscape[:, :, 2*Nsp+1] .= (r .+ f .* sum([trait_sp[k] .* neigh_matrix[:, :, k] for k in 1:Nsp]) ./ z) .* (landscape .== -1)
            Rate_landscape[:, :, 2*Nsp+2] .= d .* (landscape .== 0)




            propensity = [sum(Rate_landscape[:, :, k]) for k in 1:size(Rate_landscape)[3]]

            nb_events = map(x -> rand(Poisson(x)), propensity * tau_leap)

            for event in 1:length(nb_events) #for each type of events : colonization, death etc...
                patches = findall(landscape .== rules_change[event, 1])

                if nb_events[event] != 0 && length(patches) > nb_events[event]
                    landscape[wsample(patches, Rate_landscape[patches, event], nb_events[event])] .= rules_change[event, 2]
                end

            end #end event loop

            rho = [length(findall((landscape .== k))) / length(landscape) for k in 1:Nsp] #vegetation
            push!(rho, length(findall((landscape .== 0))) / length(landscape)) #empty sites
            push!(rho, length(findall((landscape .== -1))) / length(landscape)) #degraded sites
            push!(rho, dt)

            @views d2[dt, :] = rho

            Rate_landscape = zeros(nb_cell, nb_cell, 2 * Nsp + 2) #death and recruitment rates, degradation + restoration

        end #end time loop










    end

    return d2, landscape

end

function Plot_dynamics_Nspecies(; d, Nsp, name_x_axis="time", PA=false, black=false)
    if black
        plot(xlabel=name_x_axis, ylabel="Densities", legend=false, bg="black")
    else
        plot(xlabel=name_x_axis, ylabel="Densities", legend=false)
    end

    colors = palette([colorant"#077D10", colorant"#2A9026", colorant"#4DA33D", colorant"#71B754",
            colorant"#94CB6B", colorant"#A8D881", colorant"#9ED894", colorant"#93D8A7", colorant"#89D8B9",
            colorant"#7ED8CC", colorant"#6FC1D6", colorant"#5E9FDC", colorant"#4C7DE3", colorant"#3B5BE9", colorant"#2A39EF"], Nsp)

    names = ["Sp" * repr(x) for x in 1:Nsp]
    push!(names, "Fertile")
    push!(names, "Degraded")

    if PA
        d = d[:, vec(hcat(transpose([k for k in 1:Nsp]), transpose([Nsp + 1, Nsp + 2, size(d, 2)])))]
    end
    if black
        for k in 1:(size(d, 2)-1)
            if k != (size(d, 2) - 1) && k <= Nsp
                plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=5, thickness_scaling=1, color=colors[k])
            elseif k == Nsp + 1
                plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=5, thickness_scaling=1, color="#F5D592")
            else
                display(plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=5, thickness_scaling=1, color="black", legend=false))
            end
        end
    else

        for k in 1:(size(d, 2)-1)
            if k != (size(d, 2) - 1) && k <= Nsp
                plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=1.5, thickness_scaling=1, color=colors[k])
            elseif k == Nsp + 1
                plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=1.5, thickness_scaling=1, color="#F5D592")
            else
                display(plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=1.5, thickness_scaling=1, color="black", legend=:bottomright))
            end
        end

    end
end

function Plot_landscape_Nspecies(landscape)
    max_land = maximum(landscape)
    landscape = Matrix{Float64}(landscape)
    landscape[findall(landscape .== -1)] .= NaN
    landscape[findall(landscape .== 0)] .= NaN

    #generated from R package
    colGRAD = cgrad([colorant"#077D10", colorant"#2A9026", colorant"#4DA33D", colorant"#71B754",
        colorant"#94CB6B", colorant"#A8D881", colorant"#9ED894", colorant"#93D8A7", colorant"#89D8B9",
        colorant"#7ED8CC", colorant"#6FC1D6", colorant"#5E9FDC", colorant"#4C7DE3", colorant"#3B5BE9", colorant"#2A39EF"])

    heatmap(landscape, yflip=true, fill=true, c=colGRAD, clim=(1, max_land))
end

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
