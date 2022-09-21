using StatsBase ,RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings, Images
using Tables, CSV, LinearAlgebra, Distributions


#region :  1-- Parameter functions

"""
"Function that create a vector of species traits"
"""

function Create_traits_species(nb_sp)
  return rand(Uniform(0,1),nb_sp)
end




"""
    Interaction_matrix()
    
    "Function that build the interaction matrix between N plant species"
"""
function Interaction_matrix(nb_sp,sp_traits)

  Adj = Matrix{Float64}(undef,nb_sp,nb_sp)
  Adj[diagind(Adj)] .= 1

  for i in 1:size(Adj)[1]
    for j in (i+1):size(Adj)[2]
      if (sp_traits[i]>sp_traits[j]) # i is less competitive compared to j 
        Adj[j,i] = rand(Uniform(0.5,2))
        Adj[i,j] = rand(Uniform(0,Adj[j,i]))
        
      else  # i is more competitive
        Adj[i,j] = rand(Uniform(0.5,2))
        Adj[j,i] = rand(Uniform(0,Adj[i,j]))
      end
    end
  end

  return Adj

end




"""
    Get_params()

    "Function that return a dictionary of parameter of the model. S is the number of species."
"""
function Get_params(;cg,cl,fs,nb_sp,z,S)
  
  # Dispersal : fraction of global dispersal
  δ0=0.2
  d=3/4 # minimal dispersal fraction when surounded by shrub = 5% of global dispersal

  # Seed production
  β=0.8

  # Mortality
  m0=0.01
  mmax= 1

  # Maximal germination rate
  ϵmax = 1

  # Facilitation & competition intensities
  fs = fs
  cg = cg
  cl = cl

  ## Aridity niche
  #Aopt = Array{Float64}(undef,S,1)  #we assume A ∈ [0,1] 1 being the harsher conditions
  #Aopt[1:N_shrub]    .=0.8
  #Aopt[N_shrub+1:S]  .=0.3

  #Tol = Array{Float64}(undef,S,1) 
  #Tol[1:N_shrub,1]     .=1.5
  #Tol[N_shrub+1:S,1]   .=1.5


  #Species traits
  traits=Create_traits_species(nb_sp)

  #interaction matrix
  Adj = Interaction_matrix(nb_sp,traits)

  # neighboring 
  z=z

  #stess & scaling constant of stress
  e=0.5
  S=0.5

  return Dict{String,Any}(

    "nb_sp"       =>    nb_sp, 
    "d"           =>    d, 
    "δ0"          =>    δ0, 
    "β"           =>    β, 
    "m0"          =>    m0, 
    "mmax"        =>    mmax, 
    "ϵmax"        =>    ϵmax, 
    "cg"          =>    cg, 
    "cl"          =>    cl, 
    "fs"          =>    fs, 
#    "Aopt"       =>    Aopt, 
#    "Tol"        =>    Tol, 
    "Adj"         =>    Adj, 
    "z"           =>    z, 
    "S"           =>    S, # Will be changed in loops : stress intensity
    "e"           =>    e, 
    "dt"          =>    1,
    "psi"         =>    traits,
  )
end



"""
    Get_initial_lattice()

    Function that build the initial lattice state of the model
"""
function Get_initial_lattice(;param,size_mat=25)
  ini_vec=sample(0:param["nb_sp"],size_mat*size_mat)

  return reshape(ini_vec,size_mat,size_mat) #reshape by columns
end


"""
    Get_facilitation_landscape()

"""

function Get_all_rates_landscape(param,landscape)
 

  rho=[length(findall((landscape .== k)))   / length(landscape) for k in 1:param["nb_sp"]] #vegetation
  pushfirst!(rho,length(findall((landscape .== 0)))   / length(landscape) ) #empty sites

  nb_cell=size(landscape)[1]
 
  #Allocating 
  Stress_landscape=Array{Float64}(undef,nb_cell,nb_cell)
  Dispersal_landscape=Array{Float64}(undef,nb_cell,nb_cell)
  Competition_landscape=Array{Float64}(undef,nb_cell,nb_cell)
  Recruitment_landscape=Array{Float64}(undef,nb_cell,nb_cell)

  spatial_grid=Spatial_grid(nb_cell) #landscape of interaction

  for  i in 1:nb_cell
    for j in 1:nb_cell

      cell=Int(spatial_grid[i,j]) #trait of plant or empty cell

      neighbors_traits = landscape[findall((spatial_grid[i,:] .== 1))] #traits of neighboring cells 
      deleteat!(neighbors_traits, findall(x->x==0,neighbors_traits)) #remove empty cells

      if cell==0 #cell is empty
        
        # 1) local stress 
        Stress_landscape[i,j] = param["S"] * (1+ param["fs"] * sum([  param["psi"][k] * (1/param["z"] ) for k in neighbors_traits  ])) #stress value to the cell
        
        # 2) dispersal
        Dispersal_landscape[i,j] = 0  #there is no vegetation so we put 0 by default
        
        # 3) competition 
        Competition_landscape[i,j] = 0 #same

        # 4) Recruitment : we compute for each sp, the probability that she colonizes this cell

        for sp in 1:param["nb_sp"]
          s
          
          Recruitment_landscape[i,j] = param["ϵmax"] * (1-Stress_landscape[i,j]*(1-param["e"]*param["psi"][cell]))    
        
        
        
        end

      else

        Stress_landscape[i,j]= param["S"] #if there is facilitation, it does not change the stress value of the cell 
        delta=param["δ0"] #* (1-param["d"] * sum([  param["psi"][k] * (1/param["z"] ) for k in neighbors_traits ])) #fraction of local dispersal

        Dispersal_landscape[i,j] =  param["β"] *( delta*rho[Int(cell+1)] + (1-delta)) #cell+1 as the first element is empty cells

        #param["β"] *  Dispersal_landscape * rho_g + (1 - param["δg"]) * neigh_g / param["z"])
        Competition_landscape[i,j] = param["cg"]*(length(findall((landscape .!= 0))) / length(landscape)) + # global competition 
                                     param["cl"] * sum(param["Adj"][k,cell] * (1 / param["z"]) for k in neighbors_traits ) # local competition

        Recruitment_landscape[i,j] = 0


      end
    end
  end

  return Stress_landscape,Dispersal_landscape,Competition_landscape,Recruitment_landscape,rho

end








#endregion
#region : 2-- Dynamics functions



"""
    CA_shift(;landscape, param)

    Main function, that do one step of the CA model of interacting shrubs & grasses
"""
function CA_network_shift(;landscape, param)
  


    
    Stress_landscape,Dispersal_landscape,Competition_landscape,Recruitment_landscape,rho=Get_all_rates_landscape(landscape,param)


    # colonization 
    colonization=  @. param["β"] *  (Dispersal_landscape * rho_g + (1 - param["δg"]) * neigh_g / param["z"]) *(niche_sp - Cg  )*param["dt"]
      
    # mortality 
    death_shrub = param["ms"] *param["dt"]
    death_grass = param["mg"] *param["dt"]

    ## regeneration & degradation
    #regeneration = @.(param["r"] + param["fs"] * neigh_s / param["z"])*param["dt"]
    #degradation = param["d"]*param["dt"] 
      
    # Apply rules
    rnum = reshape(rand(length(landscape)),Int64(sqrt(length(landscape))),Int64(sqrt(length(landscape))))# one random number between 0 and 1 for each cell
    landscape_update=copy(landscape)
    
    # New vegetation (shrub & grass)
    landscape_update[findall((landscape .==0) .& (rnum .<= colonization_grass))] .= 2
    landscape_update[findall((landscape .==0) .& (rnum .> colonization_shrub) .& (rnum .<= colonization_shrub .+ colonization_grass))] .= 1
    
    # New fertile
    landscape_update[findall((landscape .== 1) .& (rnum .<= death_shrub))] .= 0
    landscape_update[findall((landscape .== 2) .& (rnum .<= death_grass))] .= 0
      
    #Do cumsum for N species

    rho_s = length(findall((landscape .== 1)))   / length(landscape) #fraction vegetation
    rho_g = length(findall((landscape .== 2)))   / length(landscape) #fraction vegetation
    rho_0 = 1 - rho_g - rho_s
      
    return  rho_s , rho_g, rho_0 , landscape_update
    
end


"""
    Run_CA_shift(time_sim,param,landscape)

    Function that runs and save output of CA in a dataframe (matrix) here
"""
function Run_CA_shift(time_sim,param,landscape)
  
  d=Array{Float64}(undef,length(time_sim)+1,4) #Allocating
  
  rho_s = length(findall((landscape .== 1)))   / length(landscape) #fraction vegetation
  rho_g = length(findall((landscape .== 2)))   / length(landscape) #fraction vegetation
  rho_0 = 1 - rho_g - rho_s

  d[1,:]=[1 rho_s rho_g rho_0] #the dataframe 
  merge!(param,Dict("dt"=>time_sim[2]-time_sim[1])) #time step

  @inbounds for k = Base.OneTo(length(time_sim))  
    
    rho_s, rho_g ,rho_0,landscape=CA_shift(landscape=landscape, param=param)
    @views d[k+1,:] = [k+1 rho_s rho_g rho_0]

  end
  return d , landscape
end


function PA_shift(du,u,param,t)
  
  # m = - and  p = + 
  ρ_S , ρ_SS , ρ_G , ρ_GG , ρ_m , ρ_SG , ρ_Sm , ρ_Gm , ρ_mm  = u 

  #conservation equations
  ρ_0  = copy(  1 - ρ_S  - ρ_G  - ρ_m)
  ρ_G0 = copy(ρ_G - ρ_Gm - ρ_SG - ρ_GG)
  ρ_S0 = copy(ρ_S - ρ_Sm - ρ_SG - ρ_SS)
  ρ_0m = copy(ρ_m - ρ_mm - ρ_Sm - ρ_Gm)
  
  niche_sp=Compute_species_niche(param)
    
  # competition coefficients
  Cs =  param["cg"] * (ρ_G + ρ_S) + param["cl"] * (param["Adj"][1,1] * (ρ_SS/ρ_S) * ((param["z"] - 1) / param["z"]) + param["Adj"][2,1] * (ρ_SG/ρ_S) * ((param["z"]-1) / param["z"]))
  Cg =  param["cg"] * (ρ_G + ρ_S) + param["cl"] * (param["Adj"][1,2] * (ρ_GG/ρ_G) * ((param["z"] - 1) / param["z"]) + param["Adj"][2,2] * (ρ_SG/ρ_G) * ((param["z"]-1) / param["z"]))

  # colonization 
  colonization_shrub =  (param["βs"] * (param["δs"] * ρ_S + (1 - param["δs"]) * (ρ_S0/ρ_0) * ((param["z"] - 1) / param["z"])) * (niche_sp[1] - Cs))[1]
  colonization_grass =  (param["βg"] * (param["δg"] * ρ_G + (1 - param["δg"]) * (ρ_G0/ρ_0) * ((param["z"] - 1) / param["z"])) * (niche_sp[2] - Cg))[1]



  #dρ_S 
  du[1] =  colonization_shrub * ρ_0 - param["ms"] * ρ_S

  #dρ_G 
  du[2] =  colonization_grass * ρ_0 - param["mg"] * ρ_G

  #dρ_- 
  du[3] =  param["d"] * ρ_0 - (param["r"] + param["fs"] * (ρ_Sm / ρ_m)) * ρ_m 

  #dp_SG
  du[4] =  colonization_grass * ρ_G0 + colonization_shrub * ρ_S0 - (param["ms"] +    param["mg"]) * ρ_SG

  #dp_SS
  du[5] =  2 * (colonization_shrub + ((1 - param["δs"]) / param["z"])) * ρ_S0       - 2 * param["ms"] * ρ_SS

  #dp_GG
  du[6] =  2 * (colonization_grass + ((1 - param["δg"]) / param["z"])) * ρ_G0       - 2 * param["mg"] * ρ_GG

  #dp_mm
  du[7] =  2 * param["d"] * ρ_0m  - 2 * (param["r"] * param["fs"] * ((param["z"] - 1) / param["z"]) * (ρ_Sm / ρ_m)) * ρ_mm 

  #dp_Sm
  du[8] =    param["d"] * ρ_S0  + colonization_shrub * ρ_0m - param["ms"] * ρ_Sm - (param["r"] + param["fs"] * ((param["z"] - 1) / param["z"]) * (ρ_Sm / ρ_m) + param["fs"] / param["z"]) * ρ_Sm

  #dp_Gm
  du[9] =    param["d"] * ρ_G0  + colonization_grass * ρ_0m - param["mg"] * ρ_Gm - (param["r"] + param["fs"] * ((param["z"] - 1) / param["z"]) * (ρ_Sm / ρ_m)) * ρ_Gm


end

function Reorder_dynamics(sol)
  sol_dyn =reduce(hcat,sol.u)'     #transforming solutions in clean matrix
  sol_dyn=cat(sol_dyn,sol.t,dims=2)
  sol_dyn=sol_dyn[:,[10, 1, 2, 3 ,4,5,6,7,8,9]]
  
  return sol_dyn
end

#endregion
#region 3-- Plot functions



function Plot_landscape(landscape)
  @rput landscape
  R"p=image(landscape,xaxt='n',yaxt='n')"
  return 
end


function Plot_dynamics(d,type)
  color_plot=["forestgreen","lightgreen","#EC9B6A","black"]
  label_plot=["S","G","F","D"]
  plot(xlabel=type,ylabel="Densities",legend=false)
  
  for k in 2:size(d,2)
    if k!=size(d,2)
      plot!(d[:,1],d[:,k],label=label_plot[k-1],color=color_plot[k-1],lw=2)
    else 
      display(plot!(d[:,1],d[:,size(d,2)],label=label_plot[size(d,2)-1],color=color_plot[k-1],legend=:topleft,lw=2))
    end
  end
  

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
  return(A)
end


#endregion
