rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#










tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
Nsim=10000
S_seq = seq(0, 1, length.out = Nsim)
c_seq = seq(0,.3, length.out = Nsim)
epsilon=10^(-9)
disp_seq=c(.1,.9)
scale_comp_seq=c("Local","Global")

d_RNE=d_Net_effect=tibble()


for (scale_c in scale_comp_seq){
  
  for (disp in disp_seq){
    
    for (S in S_seq) {
      
      for (comp in c_seq) {
        
        # 1) RII 
        
        param=Get_PA_parameters()
        param["S"] = S
        param["alpha_0"] = comp
        param["delta"]=disp
        
        
        
        
        #only competitive species
        state =Get_PA_initial_state(c(0,.8,.1,.1))
        julia_assign("state", state)
        julia_assign("p", param)
        
        if (scale_c == "Local"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }else{  
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d_RNE = rbind(d_RNE, d[nrow(d),c(1,2) ] %>% add_column(S = S,
                                                               alpha_0 = comp,
                                                               Sp="Competitive",
                                                               Scale_comp=scale_c,
                                                               Dispersal=disp))
        
        
        
        #only stress tolerant one
        state =Get_PA_initial_state(c(0.8,0,.1,.1))
        
        julia_assign("state", state)
        
        if (scale_c == "Local"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }else{  
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d_RNE = rbind(d_RNE, d[nrow(d),c(1,2) ] %>% add_column(S = S,
                                                               alpha_0 = comp,
                                                               Sp="Stress-tolerant",
                                                               Scale_comp=scale_c,
                                                               Dispersal=disp))
        
        
        
        
        #both species
        
        state =Get_PA_initial_state(c(0.4,.4,.1,.1))
        julia_assign("state", state)
        
        if (scale_c == "Local"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }else{  
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d_RNE = rbind(d_RNE, d[nrow(d),c(1,2) ] %>% add_column(S = S, 
                                                               alpha_0 = comp,
                                                               Sp="Both",
                                                               Scale_comp=scale_c,
                                                               Dispersal=disp))
        
        
        
        
        # 2) Net-effects 
        
        param=Get_PA_parameters()  
        param=c(param[1:3],"beta1"=1,"beta2"=1,param[5:12])
        
        for (sp in 1:2) { #for both species
          
          param["S"] = S
          param["alpha_0"] = comp # update parameters
          julia_assign("p", param)
          
          
          # first without press perturbation
          
          if (scale_c == "Local"){ 
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F_press, state, tspan, p)")
            
          }else{  
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
            
          }
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
          colnames(mean_densities)="Densities"
          
          d_Net_effect = rbind(d_Net_effect, mean_densities %>% add_column(S = S, 
                                                                           alpha_0 = comp, 
                                                                           Type = "control",
                                                                           Species=c(1,2)[sp],
                                                                           Scale_comp=scale_c,
                                                                           Dispersal=disp)) 
          
          
          # with press perturbation on growth rate proxy
          param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
          julia_assign("p", param)
          
          if (scale_c == "Local"){ 
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F_press, state, tspan, p)")
            
          }else{  
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
            
          }
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
          colnames(mean_densities)="Densities"
          
          d_Net_effect = rbind(d_Net_effect, mean_densities %>% add_column(S = S, 
                                                                           alpha_0 = comp, 
                                                                           Type = "press",
                                                                           Species=c(1,2)[sp],
                                                                           Scale_comp=scale_c,
                                                                           Dispersal=disp))
          
          param[paste0("beta",sp)] = 1
        }
      }
    }
  }
}

d_RNE$rho_plus = d_RNE$rho_1 + d_RNE$rho_2


d_compe=filter(d_RNE,Sp=="Competitive")
d_stesstol=filter(d_RNE,Sp=="Stress-tolerant")
d_both=filter(d_RNE,Sp=="Both")

d_RNE=tibble(S=d_both$S,alpha_0=d_both$alpha_0, 
             RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1), #actually its RII
             RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2),
             NintA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
             NintA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)))
d_RNE[,3:6][is.na(d_RNE[,3:6])] = NA



d_Net_effect[d_Net_effect < epsilon] = 0
colnames(d_Net_effect) = c("Eq", "S", "alpha_0", "Type","Species","Branches")

net_effect =sapply(seq(1, nrow(d_Net_effect) , by = 2),function(x){
  return((d_Net_effect$Eq[x+1] - d_Net_effect$Eq[x]) / epsilon)
})

d_net=d_Net_effect%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect

#for latex format in facet
appender <- function(string) {
  TeX(paste("$\\alpha_e = $", string))  
}


d_all=list(d_net,d_RNE)
save(file="../Table/2_species/PA/Comparing_net_effects.RData",d_all)
