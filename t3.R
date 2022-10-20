
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()






tspan = c(0, 7000) #to avoid long transient
t = seq(0, 7000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 100
N_rep2 = 10
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(0, .3, length.out = N_rep2)
f_seq=seq(0,1,length.out=N_rep2)
delta_seq=c(.1, .9)
h_seq=c(1,1.5)[1]


name_scena=c("global_C_local_F","global_C_global_F")

d_niche=d_RNE=d_all_dyn=tibble() #initializing the tibble

for (scena_ID in 1:2){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (h in h_seq){ #varying competitive advantage strength
      
      for (f in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) {
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          
          param["cintra"]=.3
          param["f"]=f
          param["h"]=h
          param["delta"]=disp
          param["alpha_0"]=alpha0
          
          
          
          # 1) Competitive species alone
          
          state=Get_PA_initial_state(c(0,.8,.1))
          julia_assign("state", state)
          d2 = d3 = tibble()
          
          for (S in S_seq) {
            
            param["S"] = S
            
            if (scena_ID==1){ #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
              
            }else if (scena_ID==2){  #global C, global F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Competitive"))
            d3 = rbind(d3, d[nrow(d), ] %>% add_column(S = S,Sp="Competitive",alpha_0=alpha0))
            
          }
          d2[d2 < 10^-4] = 0
          d3[d3 < 10^-4] = 0
          
          S_critic2_alone=abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) #range of values where competitive species is
          
          
          
          # 2) Stress-tolerant species alone
          
          state=Get_PA_initial_state(c(0.8,0,.1))
          julia_assign("state", state)
          
          d2 = tibble()
          S_seq = seq(0, 1, length.out = 100)
          
          for (S in S_seq) {
            
            param["S"] = S
            julia_assign("p", param)
            
            if (scena_ID==1){ #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
              
            }else if (scena_ID==2){  #global C, global F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }
            
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Stress_tolerant"))
            d3 = rbind(d3, d[nrow(d), ] %>% add_column(S = S,Sp="Stress_tolerant",alpha_0=alpha0))
            
          }
          
          d2[d2 < 10^-4] = 0
          d3[d3 < 10^-4] = 0
          
          S_critic1_alone=abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
          
          
          
          #3) Coexisting species
          
          state=Get_PA_initial_state()
          julia_assign("state", state)
          
          #varying the global interspecific competition
          
          d2 = tibble()
          
          for (S in S_seq) { #varying the stress 
            
            param["S"] = S
            param["alpha_0"] = alpha0
            julia_assign("p", param)
            
            if (scena_ID==1){ #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
              
            }else if (scena_ID==2){  #global C, global F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, Sp="Both"))
            d3 = rbind(d3, d[nrow(d), ] %>% add_column(S = S,Sp="Both",alpha_0=alpha0))
            
          }
          d2[d2 < 10^-4] = 0
          d3[d3 < 10^-4] = 0
          colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
          d2$rho_plus = d2$rho_1 + d2$rho_2
          
          S_critic1_both = abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
          S_critic2_both = abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) 
          
          
          
          
          #Putting niche in the tibble
          
          d_niche=rbind(d_niche,tibble(
            Facilitation = f, Local_disp = disp, H = h, Scena = name_scena[scena_ID], alpha_0 =alpha0,
            Delta_niche_1 = 100 * (S_critic1_both - S_critic1_alone) / S_critic1_alone,  #making it a percentage of initial niche
            Delta_niche_2 = 100 * (S_critic2_both - S_critic2_alone) / S_critic2_alone)) #making it a percentage of initial niche
          
          
          
          d3$rho_plus = d3$rho_1 + d3$rho_2
          
          
          d_compe=filter(d3,Sp=="Competitive")
          d_stesstol=filter(d3,Sp=="Stress_tolerant")
          d_both=filter(d3,Sp=="Both")
          
          d_RNE=rbind(d_RNE,tibble(S=d_both$S,alpha_0=d_both$alpha_0, Facilitation=f,H=h,Local_disp=disp,Scena=name_scena[scena_ID], 
                                   RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1),
                                   RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2)))
          
          
          d_all_dyn=rbind(d3,d_all_dyn)
          
          
        } #end competition loop
        
      } #end facilitation loop
      
    } #end h loop
    
  } #end dispersal loop
  
} #end scenario loop
write.table(d_niche,"../Table/2_species/Niche_expansion_PA.csv",sep=";")

d_RNE[,6:7][is.na(d_RNE[,6:7])] = NA

write.table(d_RNE,"../Table/2_species/RNE_PA.csv",sep=";")
write.table(d_all_dyn,"../Table/2_species/All_dyn_PA_RNE.csv",sep=";")









### 2) Analyse + plot graphics ----

#First, niche of species
d_niche=read.table("../Table/2_species/Niche_expansion_PA.csv",sep=";")
alpha_seq=seq(0,.3,length.out=10)
d_niche$alpha_0=alpha_seq

for (scale_facil in c("global_C_local_F","global_C_global_F")){
  
  p=ggplot(d_niche%>%
             melt(.,measure.vars=c("Delta_niche_2"))%>%
             filter(., Scena==scale_facil))+
    geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
    the_theme+
    facet_grid(.~Local_disp+H,labeller = label_both)+
    scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
    labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="Niche expansion (%)")
  
  ggsave(paste0("../Figures/2_species/PA/Niche_expansion_PA_facilitation_",ifelse(scale_facil=="global_C_local_F","local","global"),".pdf"),p,width = 12,height = 4)
}


#Second, RII analysis
d_RNE=read.table("../Table/2_species/RNE_PA.csv",sep=";")

# We compute the fraction of S values where RII is positive 
# And the slope along the stress gradient

d_RNE_analysis=tibble()


for (disp in unique(d_RNE$Local_disp)){
  
  for (scena in unique(d_RNE$Scena)){
    
    for (hcomp in unique(d_RNE$H)){
      
      for (a0 in unique(d_RNE$alpha_0)){
        
        for (facil in unique(d_RNE$Facilitation)){
          
          d_fil=filter(d_RNE,Facilitation==facil,Scena==scena,H==hcomp,alpha_0==a0,Local_disp==disp)
          slope=lm(d_fil$RNE_competitive~d_fil$S)$coeff[2]
          frac_net_facilitation=length(d_fil$RNE_competitive[d_fil$RNE_competitive>0])/length(d_fil$RNE_competitive[d_fil$RNE_competitive!=0])
          d_RNE_analysis=rbind(d_RNE_analysis,
                               tibble(alpha_0=a0, Facilitation=facil,H=hcomp,
                                      Local_disp=disp,Scena=scena,
                                      Slope=slope,Frac_facil=frac_net_facilitation))
          
        }
      }
    }
  }
}


d_RNE_analysis=transform(d_RNE_analysis,
                         Local_disp = factor(Local_disp, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
                         Scena=factor(Scena,levels=c("global_C_global_F","global_C_local_F"),
                                      labels=c("Global facilitation","Local facilitation")))

p=ggplot(d_RNE_analysis%>%filter(., H==1))+
  geom_tile(aes(x=Facilitation,y=alpha_0,fill=Frac_facil))+
  the_theme+facet_grid(Scena~Local_disp,labeller=labeller(Local_disp=label_parsed))+
  scale_fill_gradientn(colours=colorRampPalette(c("#000000","#737774","#FFFFFF"))(10))+
  labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="Fraction of net positive interaction along the stress gradient ")

ggsave("../Figures/2_species/PA/Fraction_positive_RNE.pdf",width = 7,height = 6)

