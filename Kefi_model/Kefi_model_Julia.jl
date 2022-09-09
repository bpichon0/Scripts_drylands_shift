using StatsBase ,RCall, Plots, Random, Agents, DifferentialEquations,BenchmarkTools


function Get_params(;δ,b,c,m,d,r,f,z)
    return vec(Float64[δ b c m d r f z])
end

function Get_classical_param()
    return Get_params( δ = 0.1,b=0.57,c=0.2,m=0.2,d=0.1,r=0,f=0.9,z=4)
end
  
function Get_initial_lattice(;frac=[.8 ,.1, .1],size_mat=25)

    ini_vec=sample([1, 2, 3],Weights(frac),size_mat*size_mat)

    return reshape(ini_vec,size_mat,size_mat) #reshape by columns
end

#params=Get_classical_param()

#sum(landscape[Base.findall(landscape.==1)])
#sum(landscape[Base.findall(landscape.==2)])



function CA_Kefi(;landscape, param)
  
    # Variables : 1 = vegetation, 2 = fertile, 3 = degraded   
    rho_v = length(findall((landscape .== 1)))   / length(landscape) #fraction vegetation
    rho_f = length(findall((landscape .== 2)))    / length(landscape) #fraction fertile
    rho_d = 1-rho_v-rho_f # fraction degraded
    
    # Neighbors :

    #using simcol package from R 
    @rput landscape
    R"neigh_v= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_f= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    R"neigh_d= simecol::neighbors(x =landscape,state = 3, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)"
    
    @rget neigh_v
    @rget neigh_f
    @rget neigh_d

    δ,b,c,m,d,r,f,z,dt= param  

    colonization = @.(δ * rho_v + (1 - δ) * neigh_v / z) * @.(b - c * neigh_v )*dt
      
    # calculate regeneration, degradation & mortality rate
    death = m *dt
    regeneration = @.(r + f * neigh_v / z)*dt
    degradation = d*dt 
      
    # Apply rules
    rnum = reshape(rand(length(landscape)),Int64(sqrt(length(landscape))),Int64(sqrt(length(landscape))))# one random number between 0 and 1 for each cell
    landscape_update=copy(landscape)
      
    ## New vegetation
    landscape_update[findall((landscape .==2) .& (rnum .<= colonization))] .= 1
      
    ## New fertile
    landscape_update[findall((landscape .== 1) .& (rnum .<= death))] .= 2
    landscape_update[findall((landscape .== 3) .& (rnum .<= regeneration))] .= 2
      
    ## New degraded 
    landscape_update[findall((landscape .== 2) .& (rnum .> colonization) .& (rnum .<= (colonization .+ degradation)))] .= 3
  
    rho_v = length(findall((landscape_update .== 1))) / length(landscape_update)
    rho_f = length(findall((landscape_update .== 2))) / length(landscape_update)
    rho_d = length(findall((landscape_update .== 3))) / length(landscape_update)
      
    return  rho_v , rho_f ,  rho_d , landscape_update
    
end

function Run_CA_kefi(;time_sim,param,landscape)
  
  d=Array{Float64}(undef,length(time_sim)+1,4) #Allocating
  
  Rho_v = length(findall((landscape .== 1))) / length(landscape)
  Rho_f = length(findall((landscape .== 2))) / length(landscape)
  Rho_d = length(findall((landscape .== 3))) / length(landscape)

  d[1,:]=[1 Rho_v Rho_f Rho_d] #the dataframe 
  push!(param,time_sim[2]-time_sim[1]) #time step

  @inbounds for k = Base.OneTo(length(time_sim))  
    
    Rho_v,Rho_f,Rho_d,landscape=CA_Kefi(landscape=landscape, param=param)
    @views d[k+1,:] = [k+1 Rho_v Rho_f Rho_d]

  end
  
  return d , landscape
  
end




param=Get_classical_param()
fraction_cover=[.8,.1,.1]
size_landscape=25
ini_land=Get_initial_lattice( frac=fraction_cover,size_mat=size_landscape)

d,state=Run_CA_kefi(time_sim=range(1,1000,step=1),param=copy(param),landscape=copy(ini_land))
plot(d[:,1],d[:,2], seriescolor = :green,label="vegetation")
plot!(d[:,1],d[:,3], seriescolor = :grey,label="fertile")
plot!(d[:,1],d[:,4], seriescolor = :black,label="degraded")

@time @inbounds for k in 1:10 d,state=Run_CA_kefi(time_sim=range(1,1000,step=1),param=param,landscape=copy(ini_land)) end

# about 3-4 times faster compared to R simulations. Surely can be more optimized, I don't know yet julia well


param=Get_params( δ = 0.1,b=0.57,c=0.3,m=0.05,d=0.2,r=0.01,f=0.3,z=4)
fraction_cover=[.3,.4,.3];size_landscape=25
ini_land=Get_initial_lattice( frac=fraction_cover,size_mat=size_landscape)


#Simulation using CA instead of PA or MF models
Random.seed!(8)
d2=zeros(1,4)
for b in 0:0.01:0.8
  param[2]=b
  d,state=Run_CA_kefi(time_sim=range(1,1000,step=1),param=copy(param),landscape=copy(ini_land))
  d2=cat(d2,[b mean(d[700:1000,2]) mean(d[700:1000,3]) mean(d[700:1000,4])],dims=1)
end

plot(d2[:,1],d2[:,2], seriescolor = :green,label="vegetation density")



# PA 

function PA_Kefi(du,u,p,t)
  δ,b,c,m,d,r,f,z= p

  ρ_pp, ρ_pm, ρ_mm, ρ_p, ρ_m = u 

  #dp++
  du[1] = 2* (ρ_p-ρ_pm-ρ_pp)  * (δ *ρ_p  +  (1-δ)/z + ((z-1)/z)*(1-δ)*((ρ_p-ρ_pm-ρ_pp)/(1-ρ_p-ρ_m)))*
          (b-c*ρ_p) - 2 *ρ_pp*m

  #dp+-
  du[2] =  d*(ρ_p-ρ_pm-ρ_pp)+ (ρ_m-ρ_mm-ρ_pm) *(δ *ρ_p+((z-1)/z)*(1-δ)*((ρ_p-ρ_pm-ρ_pp)/(1-ρ_p-ρ_m)))*
  (b-c*ρ_p) - ρ_pm* (r+ f/z + ((z-1)/z)*f*((ρ_pm)/(ρ_m))  +m)

  #dp--
  du[3] =   2*d*(ρ_m-ρ_mm-ρ_pm) -2 *ρ_mm * (r+((z-1)/(z))*f*(ρ_pm/ρ_m))

  #dp+ 
  du[4] =  (δ *ρ_p+(1-δ)*((ρ_p-ρ_pm-ρ_pp)/(1-ρ_p-ρ_m))) * (b-c*ρ_p)*(1-ρ_p-ρ_m)-m*ρ_p

  #dp-
  du[5] =  d*(1-ρ_p-ρ_m)-(r+f*(ρ_pm/ρ_m))*ρ_m

end


u0 = [0.1 ; 0.1 ; 0.1 ; 0.8; 0.1]
tspan = (0.0,1000.0)
t=0:.1:1000
param=Get_classical_param()

prob = ODEProblem(PA_Kefi,u0,tspan,param)
sol = solve(prob,Tsit5(),saveat=t)


#Comparizon CA & PA

fraction_cover=[.3,.4,.3];size_landscape=25
ini_land=Get_initial_lattice( frac=fraction_cover,size_mat=size_landscape)
param=Get_classical_param()

d,state=Run_CA_kefi(time_sim=range(1,1000,step=1),param=copy(param),landscape=copy(ini_land))

u0 = [0.1 ; 0.1 ; 0.1 ; 0.8; 0.1]
tspan = (0.0,1000.0)
t=0:.1:1000
prob = ODEProblem(PA_Kefi,u0,tspan,param)
sol = solve(prob,Tsit5(),saveat=t)

#transforming solutions in clean matrix
sol_dyn =reduce(hcat,sol.u)'
sol_dyn=cat(sol_dyn,sol.t,dims=2)

##ploting
plot(d[:,1],d[:,2], seriescolor = :green,label="vegetation")
plot!(d[:,1],d[:,3], seriescolor = :grey,label="fertile")
plot!(d[:,1],d[:,4], seriescolor = :black,label="degraded")
plot!(sol_dyn[:,6],sol_dyn[:,4], seriescolor = :green,label="PA_V",linestyle=:dash)
plot!(sol_dyn[:,6],sol_dyn[:,5], seriescolor = :black,label="PA_D",linestyle=:dash)
plot!(sol_dyn[:,6],1 .- sol_dyn[:,4] .- sol_dyn[:,5], seriescolor = :grey,label="PA_F",linestyle=:dash)














#testing with dictionaries
u0 = [0.1 ; 0.1 ; 0.1 ; 0.8; 0.1]
tspan = (0.0,1000.0)
t=0:.1:1000


param=Dict{String,Any}(
    "δ"     =>    0.1, "b"      =>    0.57, "m"           =>    0.2, 
    "d"          =>    0.1, "c"          =>    0.2, 
    "r"          =>    0.001,  "z"          =>    4,"f"          =>    0.9)


function PA_Kefi(du,u,p,t)
      δ,b,c,m,d,r,f,z= param["δ"],param["b"],param["c"],param["m"],param["d"],param["r"],param["f"],param["z"]
    
      ρ_pp, ρ_pm, ρ_mm, ρ_p, ρ_m = u 
    
      #dp++
      du[1] = 2* (ρ_p-ρ_pm-ρ_pp)  * (δ *ρ_p  +  (1-δ)/z + ((z-1)/z)*(1-δ)*((ρ_p-ρ_pm-ρ_pp)/(1-ρ_p-ρ_m)))*
              (b-c*ρ_p) - 2 *ρ_pp*m
    
      #dp+-
      du[2] =  d*(ρ_p-ρ_pm-ρ_pp)+ (ρ_m-ρ_mm-ρ_pm) *(δ *ρ_p+((z-1)/z)*(1-δ)*((ρ_p-ρ_pm-ρ_pp)/(1-ρ_p-ρ_m)))*
      (b-c*ρ_p) - ρ_pm* (r+ f/z + ((z-1)/z)*f*((ρ_pm)/(ρ_m))  +m)
    
      #dp--
      du[3] =   2*d*(ρ_m-ρ_mm-ρ_pm) -2 *ρ_mm * (r+((z-1)/(z))*f*(ρ_pm/ρ_m))
    
      #dp+ 
      du[4] =  (δ *ρ_p+(1-δ)*((ρ_p-ρ_pm-ρ_pp)/(1-ρ_p-ρ_m))) * (b-c*ρ_p)*(1-ρ_p-ρ_m)-m*ρ_p
    
      #dp-
      du[5] =  d*(1-ρ_p-ρ_m)-(r+f*(ρ_pm/ρ_m))*ρ_m
    
end
prob = ODEProblem(PA_Kefi,u0,tspan,param)
sol = solve(prob,Tsit5(),saveat=t)
sol_dyn =reduce(hcat,sol.u)'
sol_dyn=cat(sol_dyn,sol.t,dims=2)

##ploting
plot(sol_dyn[:,6],sol_dyn[:,4], seriescolor = :green,label="PA_V",linestyle=:dash)
plot!(sol_dyn[:,6],sol_dyn[:,5], seriescolor = :black,label="PA_D",linestyle=:dash)
plot!(sol_dyn[:,6],1 .- sol_dyn[:,4] .- sol_dyn[:,5], seriescolor = :grey,label="PA_F",linestyle=:dash)

