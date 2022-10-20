rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()


tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

d2=tibble()
n_random_ini=100
S_seq = seq(0,1, length.out = 100)
c_seq=c(.3)
name_scena=c("global_C_local_F","global_C_global_F")
delta_seq=c(.1,.9)

for (disp in delta_seq) {
  
  for (scena_ID in 1:2){ #for each scenario of species pairs
    
    for (S in S_seq) { #varying dispersal scale
      
      for (n_ini in 1:n_random_ini){ # for each combination, we draw random communities
        
        state =Get_PA_initial_state(Get_MF_initial_state(type="random")[-4])
        julia_assign("state", state)
        param=Get_PA_parameters()
        param["cintra"]=.3
        param["S"] = S
        param["alpha_0"] = .3
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
        colnames(d) = c("Stress_tolerant", "Competitive", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d2 = rbind(d2, d[nrow(d), c(1,2)] %>% add_column(S = S, alpha_0 = .3,N_ini=n_ini,Scena=name_scena[scena_ID],Delta=disp))
      }
    }
  }
}


d2[d2 < 10^-4] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2
colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","N_ini","Scena","Delta")
set.seed(123)
v_vec=runif(2)
d2$CSI=d2[,1]*v_vec[1]+d2[,2]*v_vec[2]
write.table(d2,"../Table/2_species/Random_ini_com_PA.csv",sep=";")


d2t=transform(d2,
              Delta = factor(Delta, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
              Scena=factor(Scena,levels=c("global_C_global_F","global_C_local_F"),
                           labels=c("Global facilitation","Local facilitation")))




p1=ggplot(d2t%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
  geom_point(aes(x=Stress,y=value,color=variable),size=.3)+
  facet_grid(Scena~Delta,labeller=labeller(Delta=label_parsed))+
  the_theme+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=as.character(color_rho)[c(4,2)],labels=c("Stress-tolerant","Competitive"))
p2=ggplot(d2t%>%melt(., measure.vars=c("CSI")))+
  geom_point(aes(x=Stress,y=value,color=variable),size=.3)+
  facet_grid(Scena~Delta,labeller=labeller(Delta=label_parsed))+
  the_theme+
  labs(x="Stress (S)",y="CSI",color="")+
  scale_color_manual(values="black")

p_tot=ggarrange(p1+xlab(""),p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),nrow=2,heights = c(1,.8))

ggsave("../Figures/2_species/PA/CSI_2_species.pdf",p_tot,width = 5,height = 10)

