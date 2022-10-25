rm(list = ls())
source("./2_Species_analysis_functions.R")






tspan = c(0, 8000)
t = seq(0, 8000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.3,length.out=3)
name_scena=c("global_C_local_F","global_C_global_F")
delta_seq=c(.1,.9)
branches=c("Degradation","Restoration")

for (branch in branches){
  
  if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
    state =Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
    S_seq=seq(0,1, length.out = N_sim)
  }else {
    state =Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
    S_seq=rev(seq(0,1, length.out = N_sim))
  }
  
  for (disp in delta_seq) {
    
    for (scena_ID in 1:2){ #for each scenario of species pairs
      
      for (ccomp in c_seq){
        
        for (S in S_seq) { #varying dispersal scale
          
          julia_assign("state", state)
          param=Get_PA_parameters()
          param["cintra"]=.3
          param["S"] = S
          param["alpha_0"] = ccomp
          param["delta"]=disp
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
}

d2[d2 < 10^-4] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2
d2=d2[,c(1:2,10:15)]
colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","Scena","Delta","Branches","Rho_plus")


d2t=transform(d2,
              Delta = factor(Delta, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
              Scena=factor(Scena,levels=c("global_C_global_F","global_C_local_F"),
                           labels=c("Global facilitation","Local facilitation")),
              alpha_0 = factor(alpha_0, levels=c(0,.15,.3), labels=c("alpha[0] : 0", "alpha[0] : 0.15",
                                                                   "alpha[0] : 0.3")))




p=ggplot(d2t%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
  geom_line(aes(x=Stress,y=value,linetype=Branches,color=variable))+
  facet_grid(Scena+Delta~alpha_0,labeller=labeller(Delta=label_parsed,alpha_0=label_parsed))+
  the_theme+theme(strip.text.x = element_text(size=12))+
  scale_color_manual(values=as.character(color_rho[c(4,2)]))+
  labs(x="Stress (S)",y="Densities",linetype="",color="")

ggsave("../Figures/2_species/Bifu_bistability.pdf",width = 8,height = 7)

