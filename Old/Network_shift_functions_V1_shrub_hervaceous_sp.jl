using StatsBase ,RCall, Plots, StatsPlots, Random, DifferentialEquations, LaTeXStrings, Images
using Tables, CSV, LinearAlgebra, Distributions


#region :  1-- Parameter functions

"""
"Function that create a vector of species traits"
"""

function Create_traits_species(S)
  return rand(Uniform(0,1),S)
end




"""
    Interaction_matrix(S,N_shrub;equal=false)
    
    "Function that build the interaction matrix between N plant species"
"""
function Interaction_matrix(S,equal=false)

  psi=Create_traits_species(S)

  Adj = Matrix{Float64}(undef,S,S)
  Adj[diagind(Adj)] .= 1

  for i in 1:size(Adj)[1]
    for j in (i+1):size(Adj)[2]
      if (psi[i]>psi[j]) # i is less competitive compared to j 
        Adj[j,i] = rand(Uniform(0.5,2))
        Adj[i,j] = rand(Uniform(0,Adj[j,i]))
        
      else  # i is more competitive
        Adj[i,j] = rand(Uniform(0.5,2))
        Adj[j,i] = rand(Uniform(0,Adj[i,j]))
      end
    end
  end

  return Adj,psi

end




"""
    Get_params(cg,cl,fs,S,z,frac_shrub)

    "Function that return a dictionary of parameter of the model. S is the number of species."
"""
function Get_params(;cg,cl,fs,S,z,frac_shrub)
  
  # Richness : number of species in each functional group
  #N_shrub = Int64(round(S*frac_shrub))
  #N_herb  = S-N_shrub

  # Dispersal : fraction of global dispersal
  δs = .5
  δg = .5

  # Seed production
  #βs = Array{Float64}(undef,N_shrub,1) .+0.8
  #βg = Array{Float64}(undef,N_herb,1)  .+0.8
  βg,βs=0.8,0.8

  # Mortality
  #ms = Array{Float64}(undef,N_shrub,1) .+ 0.02  # 20  years
  #mg = Array{Float64}(undef,N_herb,1)  .+ 0.05   # 5 years
  mg,ms=0.05,0.02

  # Germination
  ϵs = Array{Float64}(undef,N_shrub,1) .+ 1
  ϵg = Array{Float64}(undef,N_herb,1)  .+ 1 

  # Aridity niche
  Aopt = Array{Float64}(undef,S,1)  #we assume A ∈ [0,1] 1 being the harsher conditions
  Aopt[1:N_shrub]    .=0.8
  Aopt[N_shrub+1:S]  .=0.3

  Tol = Array{Float64}(undef,S,1) 
  Tol[1:N_shrub,1]     .=1.5
  Tol[N_shrub+1:S,1]   .=1.5

  #interaction matrix
  Adj = Interaction_matrix(S,N_shrub)

  #degradation & restoration
  #d = 0.1
  #r = 0.0001  

  return Dict{String,Any}(

    "N_shrub"     =>    N_shrub, 
    "N_herb"      =>    N_herb, 
    "S"           =>    S, 
    "δs"          =>    δs, 
    "δg"          =>    δg, 
    "βs"          =>    βs, 
    "βg"          =>    βg, 
    "ms"          =>    ms, 
    "mg"          =>    mg, 
    "cg"          =>    cg, 
    "cl"          =>    cl, 
    "fs"          =>    fs, 
    #"r"           =>    r, 
    #"d"           =>    d, 
    "ϵs"          =>    ϵs, 
    "ϵg"          =>    ϵg, 
    "Aopt"        =>    Aopt, 
    "Tol"         =>    Tol, 
    "Adj"         =>    Adj, 
    "z"           =>    z, 
    "A"           =>    0.5, # Will be changed in loops 
    "dt"           =>    1    
  )
end


"""
    Set_equal_param(param)

    "Function that set all parameters equal exept for niche optimum"
"""
function Set_equal_param(param)
  param["δs"]=param["δg"]  
  param["βs"]=param["βg"]
  param["ms"]=param["mg"]
  param["ϵs"]=param["ϵg"]
  param["Aopt"]=[0.5, 0.5]
  param["Adj"] = Interaction_matrix(param["S"],param["N_shrub"],equal=true)

  return param
end







"""
    Compute_species_niche(param)

    "Function that computes the species niche and more specifically, their maximal recruitment rate"
"""
function Compute_species_niche(param,landscape_aridity,landscape)

  S,N_shrub=param["S"],param["N_shrub"]
  landscape_aridity
  niche_sp = Array{Float64}(undef,size(landscape_aridity)[1],size(landscape_aridity)[2])

  for i in 1:size(landscape_aridity)[1]
    for j in 1:size(landscape_aridity)[2]
        
      if landscape[i,j]==1 #ie a shrub
        niche_sp[i,j] = param["ϵs"][1] * exp(-((landscape_aridity[i,j]-param["Aopt"][1])/param["Tol"][1])^2)
      else
        niche_sp[i,j] = param["ϵg"][1] * exp(-((landscape_aridity[i,j]-param["Aopt"][2])/param["Tol"][2])^2)
      end
    
    end
  end
  return niche_sp
end






"""
    Get_initial_lattice(;frac=[.8 ,.1, .1],size_mat=25)

    Function that build the initial lattice state of the model
"""
function Get_initial_lattice(;frac=[.8 ,.1, .1],size_mat=25)
  frac=[frac[3], frac[2], frac[1]/2, frac[1]/2]
  ini_vec=sample([-1, 0, 1, 2],Weights(frac),size_mat*size_mat)

  return reshape(ini_vec,size_mat,size_mat) #reshape by columns
end

#endregion
#region : 2-- Dynamics functions



"""
    CA_shift(;landscape, param)

    Main function, that do one step of the CA model of interacting shrubs & grasses
"""
function CA_shift(;landscape, param)
  
    # Variables : 0 = empty, 1 = shrub, 2 = grass   
    rho_s = length(findall((landscape .== 1)))   / length(landscape) #fraction vegetation
    rho_g = length(findall((landscape .== 2)))   / length(landscape) #fraction vegetation
    rho_0 = 1 - rho_g - rho_s
    
    # Neighbors :

    # using simcol package from R 
    @rput landscape
    R"neigh_s= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_g= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_0= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    
    @rget neigh_s
    @rget neigh_g
    @rget neigh_0

    #local environmental condition
    Aridity_landscape=@.(param["A"]*(1-param["fs"]*neigh_s*(param["A"]-0)))

    niche_sp=Compute_species_niche(param,Aridity_landscape,landscape)
    

    # competition coefficients
    Cs = @.( param["cg"]*(rho_g+rho_s) + param["cl"]*(param["Adj"][1,1]*neigh_s/param["z"]+param["Adj"][2,1]*neigh_g/param["z"]))
    Cg = @.( param["cg"]*(rho_g+rho_s) + param["cl"]*(param["Adj"][1,2]*neigh_s/param["z"]+param["Adj"][2,2]*neigh_g/param["z"]))

    #grass dispersal
    delta_g = @.(param["δg"]/ (1+4*(neigh_s/param["z"])))

    # colonization 
    colonization_shrub = @. param["βs"] * (param["δs"] * rho_s + (1 - param["δs"]) * neigh_s / param["z"]) *(niche_sp - Cs  )*param["dt"]
    colonization_grass = @. param["βg"] * (param["δg"] * rho_g + (1 - param["δg"]) * neigh_g / param["z"]) *(niche_sp - Cg  )*param["dt"]
      
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



function Plot_species_niche(param)

  S,N_shrub=param["S"],param["N_shrub"]
  Tol,Aopt=param["Tol"],param["Aopt"]

  A=0:0.001:1.5
  plot(xlabel=L"Aridity (A)",ylabel=L"ϵ \ .\ \exp(-(\frac{A-A_{opt}}{T})^2)",legend=:topleft,fill=(0,:0.5))

  for i in 1:N_shrub
    plot!(A,vec(@. param["ϵs"][i]*exp(-((A-Aopt[i])/Tol[i])^2)),fillrange=0,fillalpha=0.25,label="Shrub")
  end
  for i in 1:(S-N_shrub)
    if i==(S-N_shrub)
      display(plot!(A,vec(@. param["ϵg"][i]*exp(-((A-Aopt[N_shrub+i])/Tol[N_shrub+i])^2)),fillrange=0,fillalpha=0.25,label="Grass"))
    else
      plot!(A,vec(@. param["ϵg"][i]*exp(-((A-Aopt[N_shrub+i])/Tol[N_shrub+i])^2)),fillrange=0,fillalpha=0.25,label="Grass")
    end
  end
end


function Plot_landscape(landscape)
  @rput landscape
  R"p=image(landscape,col=c('0'='#EC9B6A','1'='forestgreen','-1'='black','2'='lightgreen'),xaxt='n',yaxt='n')"
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



#endregion