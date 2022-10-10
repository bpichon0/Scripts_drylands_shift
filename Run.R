rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#

## 1) Exploration along competitive ability ----

tspan = c(0, 6000)
t = seq(0, 6000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


for (scena in c("global_comp","local_comp")){
  
  
  
  #setting intraspecific competition to .5 (maximal competition strength explored)
  param=Get_PA_parameters(type_comp = gsub("_comp","",scena))
  param["cintra"]=.5
  param["r"]=.05
  param["d"]=.1
  
  N_rep = 100
  S_seq = seq(0, 1, length.out = N_rep)
  alpha_seq = seq(0, .5, length.out = N_rep)
  
  
  
  
  #first species alone for the fundamental niche of species 
  # 1) Competitive species alone
  
  state=Get_PA_initial_state(c(0,.8,.1))
  julia_assign("state", state)
  d2 = tibble()
  
  for (S in S_seq) {
    
    param["S"] = S
    
    if (scena=="global_comp"){
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(PA_two_species_global_comp, state, tspan, p)")
    }else{
      param["cg"]=0
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(PA_two_species_local_comp, state, tspan, p)")
      
    }
    sol = de$solve(prob, de$Tsit5(), saveat = t)
    d = as.data.frame(t(sapply(sol$u, identity)))
    colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
    
    d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S))
  }
  d2[d2 < 10^-3] = 0
  
  S_critic2=d2$S[min(which(d2$rho_2==0))]
  
  
  # 2) Stress-tolerant species alone
  
  state=Get_PA_initial_state(c(0.8,0,.1))
  julia_assign("state", state)
  
  d2 = tibble()
  S_seq = seq(0, 1, length.out = 100)
  
  for (S in S_seq) {
    
    param["S"] = S
    julia_assign("p", param)
    
    if (scena=="global_comp"){
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(PA_two_species_global_comp, state, tspan, p)")
    }else{
      param["cg"]=0
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(PA_two_species_local_comp, state, tspan, p)")
      
    }
    sol = de$solve(prob, de$Tsit5(), saveat = t)
    d = as.data.frame(t(sapply(sol$u, identity)))
    colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
    
    d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S))
  }
  d2[d2 < 10^-3] = 0
  
  S_critic1=d2$S[min(which(d2$rho_1==0))]
  
  
  #main loop with coexisting species
  d2 = tibble()
  
  state=Get_PA_initial_state()
  julia_assign("state", state)
  
  for (S in S_seq) {
    
    for (alpha0 in alpha_seq) {
      
      param["S"] = S
      param["alpha_0"] = alpha0
      julia_assign("p", param)
      
      if (scena=="global_comp"){
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(PA_two_species_global_comp, state, tspan, p)")
      }else{
        param["cg"]=0
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(PA_two_species_local_comp, state, tspan, p)")
        
      }
      
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
      
      d2 = rbind(d2, as_tibble(t(get_mean_densities(d))) %>% add_column(S = S, alpha_0 = alpha0))
    }
  }
  d2[d2 < 10^-4] = 0
  colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
  d2$rho_plus = d2$rho_1 + d2$rho_2
  
  # postprocessing
  post_processing_2species_PA(d2, alpha_seq, S_seq,file_name = scena,S_critic1,S_critic2)
}
