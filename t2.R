
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()


tspan = c(0, 7000) #to avoid long transient
t = seq(0, 7000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 50
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(0, .3, length.out = N_rep)
f_seq=c(0,.3,.9)
delta_seq=c(.1, .9)

id_for_RNE = seq(1:N_rep)[c(1, length(seq(1:N_rep))/2, length(seq(1:N_rep)))]

name_scena=c("local_C_local_F","global_C_global_F",
             "local_C_global_F","global_C_local_F")

for (scena_ID in c(1:4)[c(2,4)]){ #for each scenario of species pairs
  
  for (f in f_seq){ #varying facilitation strength
    
    
    for (disp in delta_seq){
      
      
      name_fig=paste0("_d_",disp,"_f_",f,"_h_",1,"_scena_",name_scena[scena_ID])
      
      
      #setting intraspecific competition to .2 (maximal competition strength explored)
      param=Get_PA_parameters()
      
      param["cintra"]=.3
      param["f"]=f
      param["delta"]=disp
      
      
      
      
      
      # 1) Competitive species alone
      
      state=Get_PA_initial_state(c(0,.8,.1))
      julia_assign("state", state)
      d2 = tibble()
      
      for (S in S_seq) {
        
        param["S"] = S
        
        if (scena_ID==1){ #local C, local F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }else if (scena_ID==2){  #global C, global F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
          
        }else if (scena_ID==3){  #local C, global F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
          
        }else if (scena_ID==4){ #global C, local F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S))
      }
      d2[d2 < 10^-4] = 0
      
      S_critic2=d2$S[min(which(d2$rho_2==0))]
      
      
      
      
      # 2) Stress-tolerant species alone
      
      state=Get_PA_initial_state(c(0.8,0,.1))
      julia_assign("state", state)
      
      d2 = tibble()
      S_seq = seq(0, 1, length.out = 100)
      
      for (S in S_seq) {
        
        param["S"] = S
        julia_assign("p", param)
        
        if (scena_ID==1){ #local C, local F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }else if (scena_ID==2){  #global C, global F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
          
        }else if (scena_ID==3){  #local C, global F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
          
        }else if (scena_ID==4){ #global C, local F
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S))
      }
      
      d2[d2 < 10^-4] = 0
      
      S_critic1=d2$S[min(which(d2$rho_1==0))]
      
      
      #writing tipping points species alone
      write.table(tibble(Sp=c(1,2),Stress_crit=c(S_critic1,S_critic2)),
                  paste0("../Table/2_species/Sim_PA_scales/Tipping_sp",name_fig,".csv"),sep=";")
      
      
      
      
      
      
      
      #3) Coexisting species
      
      d2 = tibble()
      
      state=Get_PA_initial_state()
      julia_assign("state", state)
      
      for (S in S_seq) {
        
        for (alpha0 in alpha_seq) {
          
          param["S"] = S
          param["alpha_0"] = alpha0
          julia_assign("p", param)
          
          if (scena_ID==1){ #local C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==3){  #local C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==4){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
            
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, alpha_0 = alpha0))
        }
      }
      d2[d2 < 10^-4] = 0
      colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
      d2$rho_plus = d2$rho_1 + d2$rho_2
      
      write.table(d2,paste0("../Table/2_species/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")
      
      
      
      
      # 4) Doing the RII analysis 
      
      
      #Doing the RII analysis. We replicate analyse for clarity of the code, knowing that it will take more computation time
      
      d_RNE = tibble() 
      
      for (S in S_seq) {
        
        param["S"] = S
        
        for (comp in alpha_seq[id_for_RNE]) {
          
          param["alpha_0"] = comp
          
          #only competitive species
          state =Get_PA_initial_state(c(0,.8,.1))
          julia_assign("state", state)
          julia_assign("p", param)
          
          
          
          if (scena_ID==1){ #local C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==3){  #local C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==4){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
            
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          
          d_RNE = rbind(d_RNE, d[nrow(d), ] %>% add_column(S = S,alpha_0 = comp,Sp="Competitive"))
          d_RNE[d_RNE < 10^-4] = 0
          colnames(d_RNE) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm","S","alpha_0","Sp")
          
          #only stress tolerant one
          state =Get_PA_initial_state(c(0.8,0,.1))
          julia_assign("state", state)
          
          if (scena_ID==1){ #local C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==3){  #local C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==4){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
            
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d_RNE = rbind(d_RNE, d[nrow(d), ] %>% add_column(S = S,alpha_0 = comp,Sp="Stress-tolerant"))
          d_RNE[d_RNE < 10^-4] = 0
          
          #both species
          state =Get_PA_initial_state()
          julia_assign("state", state)
          
          if (scena_ID==1){ #local C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==3){  #local C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==4){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
            
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d_RNE = rbind(d_RNE, d[nrow(d), ] %>% add_column(S = S, alpha_0 = comp,Sp="Both"))
        }
      }
      d_RNE$rho_plus = d_RNE$rho_1 + d_RNE$rho_2
      
      write.table(d_RNE,paste0("../Table/2_species/Sim_PA_scales/2_species_RNE",name_fig,".csv"),sep=";")
      
      
    } #end dispersal loop
    
    
  } #end facilitation loop
  
} #end scenario loop
















N_rep = 50
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(0, .3, length.out = N_rep)
h_seq=c(1,1.5)[1]
f_seq=c(0,.3,.9)

id_for_RNE = seq(1:N_rep)[c(1, length(seq(1:N_rep))/2, length(seq(1:N_rep)))]

name_scena=c("local_C_local_F","global_C_global_F",
             "local_C_global_F","global_C_local_F")

dir.create("../Figures/2_species/PA/Sim_PA_scales",showWarnings = F)
dir.create("../Figures/2_species/PA/Clustering",showWarnings = F)

for (scena_ID in c(1:4)[c(2,4)]){ #for each scenario of species pairs
  
  for (f in f_seq){ #varying facilitation strength
    
    for (h in h_seq){ #varying competitive advantage strength
      
      for (disp in c(.1, .9)){
        
        name_fig=paste0("_d_",disp,"_f_",f,"_h_",h,"_scena_",name_scena[scena_ID])
        
        #loading data
        d_state=read.table(paste0("../Table/2_species/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")
        Scritic=read.table(paste0("../Table/2_species/Sim_PA_scales/Tipping_sp",name_fig,".csv"),sep=";")
        d_RNE=read.table(paste0("../Table/2_species/Sim_PA_scales/2_species_RNE",name_fig,".csv"),sep=";")
        
        
        
        #1) First classical figure with state diagram, RII and shifts
        
        d_state$state = sapply(1:nrow(d_state), function(x) {
          if (is.na(d_state$rho_1[x]) & is.na(d_state$rho_2[x])) {
            return("NA")
          }else{
            
            
            if (d_state$rho_1[x] > 0 & d_state$rho_2[x] > 0) {
              return("Coexistence")
            }
            if (d_state$rho_1[x] > 0 & d_state$rho_2[x] == 0) {
              return("Stress-tolerant")
            }
            if (d_state$rho_1[x] == 0 & d_state$rho_2[x] > 0) {
              return("Competitive")
            }
            if (d_state$rho_1[x] == 0 & d_state$rho_2[x] == 0) {
              return("Desert")
            }
          }
          
        })
        a_values_bifu =  alpha_seq[id_for_RNE]
        
        
        color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress-tolerant" = "#7BD8D3")
        
        # state at equilibrium
        p1 = ggplot(d_state) +
          geom_tile(aes(x = S, y = alpha_0 , fill = state)) +
          the_theme +
          annotate("text", x = rep(1.05, 3), y = a_values_bifu+max(a_values_bifu)/35 , label = c("e, f", "c, d", "a, b"), color = "black",size=5) +
          scale_fill_manual(values = color_rho) +
          theme(legend.position = "bottom") +
          labs(x = "Stress (S)", y = TeX(r'(Competition strength \ $\alpha_e)'), fill = "") +
          geom_hline(yintercept = a_values_bifu , lwd = .1, color = "gray40") +
          theme(legend.text = element_text(size = 12))
        
        
        
        #bifurcation diagrams
        S_critic1=Scritic[1,2]
        S_critic2=Scritic[2,2]
        
        d_bifu = filter(d_state, round(alpha_0, 4) %in% round(a_values_bifu, 4))
        for (a_bifu in 1:3) {
          assign(
            paste0("p2_", a_bifu),
            ggplot(d_bifu %>% filter(., round(alpha_0, 4) == round(a_values_bifu[a_bifu], 4)) %>% melt(., measure.vars = c("rho_1", "rho_2")) %>%
                     mutate(., variable = recode_factor(variable, "rho_1" = "Stress-tolerant", "rho_2" = "Competitive"))) +
              geom_point(aes(x = S, y = value, color = variable), size = .4) +
              geom_segment(
                x = S_critic2, y = .3+.1,xend = S_critic2, yend = .3,
                lineend = "round",linejoin = "round", size = .3, 
                arrow = arrow(length = unit(0.1, "inches")),
                colour = color_rho[2]) + 
              geom_segment(
                x = S_critic1, y = .3+.1,xend = S_critic1, yend = .3,
                lineend = "round",linejoin = "round", size = .3, 
                arrow = arrow(length = unit(0.1, "inches")),
                colour = color_rho[4])+
              labs(x = "Stress (S)", y = "Density", color = "") +
              the_theme +
              scale_color_manual(values = color_rho[c(2, 4)]) +
              theme(legend.position = "none")
          )
        }
        
        #RII
        
        d_compe=filter(d_RNE,Sp=="Competitive")
        d_stesstol=filter(d_RNE,Sp=="Stress-tolerant")
        d_both=filter(d_RNE,Sp=="Both")
        
        d_RNE=tibble(S=d_both$S,alpha_0=d_both$alpha_0, 
                     RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1), #acutally its RII
                     RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2),
                     NintA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
                     NintA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)))
        d_RNE[,3:6][is.na(d_RNE[,3:6])] = NA
        # d_RNE$RNE_stress_tol[246] = NA #to avoidnumeric problems
        
        
        
        for (a0 in 1:3) {
          assign(
            paste0("p_", a0),
            ggplot(d_RNE %>% filter(., round(alpha_0,4) == round(a_values_bifu[a0],4)) %>%melt(., measure.vars=c("NintA_st","NintA_comp")) %>%
                     mutate(., variable = recode_factor(variable, "NintA_st" = "Stress-tolerant", "NintA_comp" = "Competitive"))) +
              geom_line(aes(x=S,y=value,color=variable),lwd=.75)+the_theme+
              geom_hline(yintercept = 0,linetype=9)+
              geom_segment(
                x = S_critic2, y = -1+.3,xend = S_critic2, yend = -1,
                lineend = "round",linejoin = "round", size = .3, 
                arrow = arrow(length = unit(0.1, "inches")),
                colour = color_rho[2]) + 
              geom_segment(
                x = S_critic1, y = -1+.3,xend = S_critic1, yend = -1,
                lineend = "round",linejoin = "round", size = .3, 
                arrow = arrow(length = unit(0.1, "inches")),
                colour = color_rho[4]) + 
              theme(legend.text = element_text(size = 12))+
              labs(x="Stress (S)",y="RII",color="")+ylim(-1,1)+
              scale_color_manual(values=c(as.character(color_rho)[c(4,2)]))+
              theme(legend.text = element_text(hjust = 5))
          )
        }
        
        pright_1=ggarrange(p2_3+xlab(""),p_3+xlab("")+theme(legend.position = "none"),ncol=2,align = 'v',labels = letters[c(1,2)])
        pright_2=ggarrange(p2_2+xlab(""),p_2+xlab("")+theme(legend.position = "none"),ncol=2,align = 'v',labels = letters[c(3,4)])
        pright_3=ggarrange(p2_1+theme(legend.text = element_text(size=12))+  guides(color = guide_legend(override.aes = list(size = 3))),
                           p_1+theme(legend.text = element_text(size=12))+  guides(color = guide_legend(override.aes = list(size = 3))),
                           ncol=2,align = 'h',labels = letters[c(5,6)],common.legend = T,legend = "bottom")
        pright_tot=ggarrange(pright_1,pright_2,pright_3,nrow = 3,heights = c(1,1,1.2),align = "hv")
        
        
        
        p_tot = ggarrange(p1,pright_tot, ncol = 2, widths = c(2, 2))
        
        ggsave(paste0("../Figures/2_species/PA/Sim_PA_scales/Fig",name_fig,".pdf"),width = 12,height = 6)
        
        
        
      }
    }
  }
}