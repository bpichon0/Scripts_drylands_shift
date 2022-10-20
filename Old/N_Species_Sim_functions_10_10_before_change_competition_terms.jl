using StatsBase, RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings,
    BenchmarkTools, Images, Tables, CSV, LinearAlgebra, Distributions, DataFrames


#Parameters functions 

function Get_classical_param(; N_species=4, type_interaction="nested", relative_competition=4, scenario_trait="random")

    r = 0.05
    d = 0.1
    f = 0.9
    beta = 0.8
    m = 0.1
    e = 0.1
    emax = 1.2
    cg = 0.1
    S = 0
    delta = 0.1
    z = 4
    tau_leap = 0.05
    trait = Get_species_traits(Nsp=N_species, scena=scenario_trait)

    Interaction_mat = Get_interaction_matrix(Nsp=N_species, type=type_interaction, rela_comp=relative_competition, trait_sp=trait)

    return Dict{String,Any}(
        "Nsp" => N_species,
        "r" => r,
        "d" => d,
        "f" => f,
        "beta" => beta,
        "m" => m,
        "e" => e,
        "emax" => emax,
        "cg" => cg,
        "S" => S,
        "delta" => delta,
        "z" => z,
        "tau_leap" => tau_leap,
        "trait" => trait,
        "alpha" => Interaction_mat
    )
end

function Get_interaction_matrix(; Nsp, type, trait_sp, rela_comp) #A discuter comment définir les interactions

    if (type == "random")
        Interaction_mat = rand(Nsp * Nsp) / (Nsp + 3)
        Interaction_mat = reshape(Interaction_mat, Nsp, Nsp)
    elseif (type == "equal")
        Interaction_mat = zeros(Nsp, Nsp) .+ 0.1

    elseif (type == "nested")
        Interaction_mat = zeros(Nsp, Nsp) .+ 0.1 #for intraspecific competition
        for i in 1:Nsp
            for j in 1:Nsp
                if (trait_sp[i] > trait_sp[j]) #i.e., less stress tolerant plant
                    Interaction_mat[j, i] += Interaction_mat[j, i] * abs(trait_sp[i] - trait_sp[j]) * rela_comp
                    #Interaction_mat[j, i] += Interaction_mat[j, i] * exp(-abs(trait_sp[i] - trait_sp[j])) * rela_comp
                end
            end
        end

    elseif type == "low_inter"
        Interaction_mat = zeros(Nsp, Nsp) .+ 0.05 #for interspecific competition
        for i in 1:size(Interaction_mat)[1]
            Interaction_mat[i, i] = 0.2
        end
        for i in 1:Nsp
            for j in 1:Nsp
                if (trait_sp[i] > trait_sp[j]) # i is more stress tolerant but less competitive than j
                    Interaction_mat[j, i] += Interaction_mat[j, i] * abs(trait_sp[i] - trait_sp[j]) * rela_comp
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

function Get_initial_lattice(; param, branch="Degradation", size_mat=25, type_ini="equal")

    frac = Get_initial_state(param=param, type=type_ini, branch=branch)
    species_vec = [k for k in 1:param["Nsp"]]
    push!(species_vec, 0)
    push!(species_vec, -1)
    ini_vec = sample(species_vec, Weights(frac), size_mat * size_mat)

    return reshape(ini_vec, size_mat, size_mat) #reshape by columns
end

function Get_initial_state(; param, type, branch)
    if branch == "Degradation"  #degradation branch

        if type == "equal"
            frac = [0.8 / param["Nsp"] for k in 1:param["Nsp"]]
        elseif type == "random"
            frac = rand(param["Nsp"])
            frac = (frac / sum(frac)) * 0.8
        end
        push!(frac, 0.1)
        push!(frac, 0.1)

    else #restoration branch

        if type == "equal"
            frac = [0.01 / param["Nsp"] for k in 1:param["Nsp"]]
        elseif type == "random"
            frac = rand(param["Nsp"])
            frac = (frac / sum(frac)) * 0.01
        end
        push!(frac, 0.50)
        push!(frac, 0.49)
    end
    return frac
end

# MF & PA functions 
function MF_N_species(du, u, p, t)
    rho_0 = u[p["Nsp"]+1]
    rho_d = u[p["Nsp"]+2]

    for k in 1:p["Nsp"]
        du[k] = rho_0 * (p["beta"] * u[k] * (p["emax"] * (1 - p["S"] * (1 - p["trait"][k] * p["e"])) - (p["cg"] * sum(u[1:p["Nsp"]]) + sum([p["alpha"][i, k] * u[i] for i in 1:p["Nsp"]])))) - u[k] * p["m"]
    end
    du[p["Nsp"]+1] = -p["d"] * rho_0 + rho_d * (p["r"] + p["f"] * (sum([p["trait"][i] * u[i] for i in 1:p["Nsp"]])))

    for k in 1:p["Nsp"]
        du[p["Nsp"]+1] = du[p["Nsp"]+1] - rho_0 * (p["beta"] * u[k] * (p["emax"] * (1 - p["S"] * (1 - p["trait"][k] * p["e"])) - (p["cg"] * sum(u[1:p["Nsp"]]) + sum([p["alpha"][i, k] * u[i] for i in 1:p["Nsp"]])))) + u[k] * p["m"]
    end

    du[p["Nsp"]+2] = p["d"] * rho_0 - rho_d * (p["r"] + p["f"] * (sum([p["trait"][i] * u[i] for i in 1:p["Nsp"]])))

end

function Reorder_dynamics(sol)
    sol_dyn = reduce(hcat, sol.u)'     #transforming solutions in clean matrix
    sol_dyn = cat(sol_dyn, sol.t, dims=2)
    return sol_dyn
end

function PA_N_species(du, u, p, t)

    f = p["f"]
    e = p["e"]
    tau_leap = p["tau_leap"]
    Nsp = p["Nsp"]
    emax = p["emax"]
    cg = p["cg"]
    alpha = p["alpha"]
    r = p["r"]
    delta = p["delta"]
    beta = p["beta"]
    S = p["S"]
    m = p["m"]
    trait_sp = p["trait"]
    z = p["z"]
    d = p["d"]

    # u : 

    #rho_1
    du[1] = (1 - rho_1 - rho_2 - rho_m) * beta * (delta * rho_1 + (1 - delta) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S * e) - (cg * (rho_1 + rho_2) + (alpha11 * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m)) + alpha21 * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m))))) - rho_1 * m

    #rho_2
    du[2] = (1 - rho_1 - rho_2 - rho_m) * beta * (delta * rho_2 + (1 - delta) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m)) + alpha12 * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m))))) - rho_2 * m

    #rho_m
    du[3] = (1 - rho_1 - rho_2 - rho_m) * d - rho_m * (r + f * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m)))

    #rho_12
    du[4] = ((rho_1 - rho_11 - rho_12 - rho_1m)) * beta * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m)) + (alpha12 / z) + alpha12 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m))))) +
            ((rho_2 - rho_22 - rho_12 - rho_2m)) * beta * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S * e) - (cg * (rho_1 + rho_2) + (alpha11 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m)) + (alpha21 / z) + alpha21 * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m))))) -
            2 * rho_12 * m

    #rho_1m
    du[5] = ((rho_1 - rho_11 - rho_12 - rho_1m)) * d + (rho_m - rho_mm - rho_1m - rho_2m) * beta * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S * e) - (cg * (rho_1 + rho_2) + (alpha11 * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m)) + alpha21 * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m))))) -
            rho_1m * m - rho_1m * (r + f * ((1 / z) + ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m))))

    #rho_2m
    du[6] = ((rho_2 - rho_22 - rho_12 - rho_2m)) * d + (rho_m - rho_mm - rho_1m - rho_2m) * beta * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1 - rho_1 - rho_2 - rho_m)) * (emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m)) + alpha12 * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m)))))) -
            rho_2m * m - rho_2m * (r + f * (((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m))))

    #rho_11
    du[7] = 2 * ((rho_1 - rho_11 - rho_12 - rho_1m)) * beta * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S * e) - (cg * (rho_1 + rho_2) + (alpha11 * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m)) * ((z - 1) / z) + (alpha11 / z) + alpha21 * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m))))) -
            2 * rho_11 * m

    #rho_22
    du[8] = 2 * ((rho_2 - rho_22 - rho_12 - rho_2m)) * beta * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1 - rho_1 - rho_2 - rho_m))) * (emax * (1 - S) - (cg * (rho_1 + rho_2) + (alpha22 * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1 - rho_1 - rho_2 - rho_m)) * ((z - 1) / z) + (alpha22 / z) + alpha12 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1 - rho_1 - rho_2 - rho_m))))) -
            2 * rho_22 * m

    #rho_mm
    du[9] = 2 * (rho_m - rho_mm - rho_1m - rho_2m) * d - 2 * rho_mm * (r + f * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1 - rho_1 - rho_2 - rho_m)))


end


#CA functions
function Get_neighbors_matrix(; landscape, sp)
    @rput sp
    @rput landscape
    R"neighbors_mat = simecol::neighbors(x =landscape,state = sp, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    @rget neighbors_mat
end

function Gillespie_CA_N_species(; param, landscape, tmax)


    r = param["r"]
    d = param["d"]
    f = param["f"]
    interaction_mat = param["alpha"]
    beta = param["beta"]
    m = param["m"]
    e = param["e"]
    emax = param["emax"]
    cg = param["cg"]
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


    for dt in 1:tmax

        neigh_matrix = zeros(nb_cell, nb_cell, Nsp)
        for k in 1:Nsp
            neigh_matrix[:, :, k] = Get_neighbors_matrix(landscape=copy(landscape), sp=copy(k))
        end

        for sp in 1:Nsp
            Rate_landscape[:, :, sp] .= beta .* (delta * rho[sp] .+ (1 - delta) * neigh_matrix[:, :, sp] / z) .* #dispersal
                                        (emax * (1 - S * (1 - trait_sp[sp] * e)) .- #recruitment
                                         (cg * sum(rho[1:Nsp]) .+ (sum([interaction_mat[k, sp] * #global competition
                                                                        neigh_matrix[:, :, k] for k in 1:Nsp]) / z))) .* (landscape .== 0) #local competition

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



    return d2, landscape

end

function Plot_dynamics(; d, Nsp, name_x_axis)
    plot(xlabel=name_x_axis, ylabel="Densities", legend=false)

    names = ["Sp" * repr(x) for x in 1:Nsp]
    push!(names, "Fertile")
    push!(names, "Degraded")

    for k in 1:(size(d, 2)-1)
        if k != (size(d, 2) - 1)
            plot!(d[:, size(d, 2)], d[:, k], label=names[k], lw=1.5, thickness_scaling=1)
        else
            display(plot!(d[:, size(d, 2)], d[:, k], label=names[k], legend=:topleft,
                lw=1.5, thickness_scaling=1))
        end
    end
end

function Plot_landscape(landscape)
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


function Randomization_landscape(; landscape, niter=10000)

    for i in 1:niter

        cell_1 = sample(1:size(landscape)[1], 2)
        cell_2 = sample(1:size(landscape)[1], 2)

        save = landscape[cell_1[1], cell_1[2]]
        landscape[cell_1[1], cell_1[2]] = landscape[cell_2[1], cell_2[2]]
        landscape[cell_2[1], cell_2[2]] = save

    end

    return landscape

end


function Get_species_pairs(; landscape, Nsp=15)

    landscape[findall(landscape .== -1)] .= 0
    neighbor_array = zeros(size(landscape)[1], size(landscape)[1], Nsp + 1)

    #get neighbors of each cell
    for sp in 1:(Nsp+1)
        neighbor_array[:, :, sp] = Get_neighbors_matrix(landscape=copy(landscape), sp=copy(sp - 1)) #to account for fertile stes
    end

    #for each species get the mean number of times species i is in its surrounding
    d = zeros(Nsp + 1, Nsp + 1)

    for sp in 0:Nsp
        mean_neigh_cells = mean(neighbor_array[findall(landscape .== sp), :], dims=1) #mean nuber of neighbors per species !
        d[sp+1, :] = mean_neigh_cells
    end

    return d
end



function Compute_z_score_pairs(; final_state, Nsp, N_null)

    null_matrices = zeros(size(final_state)[1], size(final_state)[1], N_null)

    #generating randomized landscapes
    for i in 1:N_null
        null_matrices[:, :, i] = Randomization_landscape(landscape=copy(final_state), niter=100000)
    end


    observed_neigh = Get_species_pairs(landscape=copy(final_state), Nsp=Nsp)
    Plot_landscape(observed_neigh)


    null_neigh = zeros(size(observed_neigh)[1], size(observed_neigh)[1], N_null)
    for i in 1:size(null_matrices)[3]
        null_neigh[:, :, i] = Get_species_pairs(landscape=copy(null_matrices[:, :, i]), Nsp=Nsp)
    end
    mean_neigh_null = mean(null_neigh, dims=3)
    sd_neigh_null = std(null_neigh, dims=3)
    z_score_neigh = (observed_neigh .- mean_neigh_null) ./ (sd_neigh_null)

    return z_score_neigh[:, :, 1]

end
