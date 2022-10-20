
tspan = c(0, 5000)
t = seq(0, 5000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()

N_sim=50
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.3,length.out=N_sim)[100]
name_scena=c("global_C_local_F","global_C_global_F")[1]
delta_seq=c(.1,.9)[1]
branches=c("Degradation","Restoration")

for (branch in branches){
  d2=tibble()
  
  if (branch =="Degradation"){
    state =Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
  }else {
    state =Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49))[-1])
  }
  
  for (disp in delta_seq) {
    
    for (scena_ID in 1){ #for each scenario of species pairs
      
      for (ccomp in c_seq){
        
        for (S in S_seq) { #varying dispersal scale
          
          julia_assign("state", state)
          param=Get_PA_parameters()
          param["cintra"]=.3
          param["S"] = S
          param["alpha_0"] = .3
          param["delta"]=disp
          param["e"]=.1
          param["r"]=0.01
          param["d"]=0.025
          param["beta"] = 1
          param["m"]=.15
          julia_assign("p", param)
          
          if (scena_ID==1){ 
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
            
          }else{  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
            
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = ccomp,
                                                     Scena=name_scena[scena_ID],
                                                     Delta=disp,Branch=branch))
          
        }
      }
    }
  }
  
  if (branch == "Degradation"){
    plot(d2$S,d2$rho_1,type="l",col="green")
    lines(d2$S,d2$rho_2,col="red")
    
  }else{
    
    lines(d2$S,d2$rho_1,col="forestgreen")
    lines(d2$S,d2$rho_2,col="blue")
  }
  
  
  
  
}

