# Other) Recruitment rate & 4 grid matrices----
rm(list = ls())
library(tidyverse)
d = expand_grid(emax = 1, e = .1, S = seq(0, 1, length.out = 200), psi = c(0, .5, 1))

p = ggplot(d) +
    geom_line(aes(x = S, y = emax * (1 - S * (1 - e * psi)), color = as.factor(psi), group = psi)) +
    scale_color_manual(values = c("blue", "green", "red")) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
        x = "S",
        y = TeX("$(1-S (1-e \\psi_i))$"), color = TeX("$\\psi_i$")
    )



ggsave("../Figures/Recruitment_rate.pdf", p, width = 6, height = 4)



# Step 1) Two species MF model ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#



## 1) State diagram with competitive ability and stress for different type of competition ----



# Intra + interspecific competition


# preparing the 3 types of simulations


tspan = c(0, 5000)
t = seq(0, 5000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
Nsim=100
S_seq = seq(0, 1, length.out = Nsim)
c_seq = seq(0,.3, length.out = Nsim)
c_rne=c_seq[c(1,Nsim/2,Nsim)]
#Doing the simulation

for (f in c(0,.3,.9)){
  
  for (h in c(1,1.5)){
    
    name_fig=paste0("_f_",f,"_h_",h)
    
    # 1) Competitive species alone
    
    state =c(0,.8,.1,.1)
    julia_assign("state", state)
    param=Get_MF_parameters()    
    
    param["f"]=f
    param["h"]=h
    
    d2 = tibble()
    
    for (S in S_seq) {
      
      param["S"] = S
      param["cintra"] =  .3
      julia_assign("p", param)
      
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S))
    }
    d2[d2 < 10^-4] = 0
    colnames(d2) = c("rho_1", "rho_2", "rho_d", "rho_0", "S")
    
    S_critic2=d2$S[min(which(d2$rho_2==0))]
      
    
    # 2) Stress-tolerant species alone
    
    state =c(0.8,0,.1,.1)
    julia_assign("state", state)
    
    d2 = tibble()
    for (S in S_seq) {
      
      param["S"] = S
      param["cintra"] =  .3
      julia_assign("p", param)
      
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S))
    }
    d2[d2 < 10^-4] = 0
    colnames(d2) = c("rho_1", "rho_2", "rho_d", "rho_0", "S")
    
    S_critic1=d2$S[min(which(d2$rho_1==0))]
    

    
    # 3) Simulation with both species
    state =Get_MF_initial_state()
    julia_assign("state", state)
    d2 = tibble()
    
    for (scena in c( "main")){
      
      
      for (S in S_seq) {
        for (comp in c_seq) {
          
          param["S"] = S
          param["alpha_0"] = comp
          julia_assign("p", param)
          
          prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = comp))
        }
      }
      d2[d2 < 10^-4] = 0
      colnames(d2) = c("rho_1", "rho_2", "rho_d", "rho_0", "S", "alpha_0")
      d2$rho_plus = d2$rho_1 + d2$rho_2
      
    }
    
    
    
    
    #compute RII = biomass species with-without/(with+without) neighbors
    
    
    d_RNE=tibble()
    
    
    S_seq = seq(0, 1, length.out = 100)
    
    
    for (S in S_seq) {
      param["S"] = S
      
      for (comp in c_rne) {
        
        param["alpha_0"] = comp
        
        #only competitive species
        state =c(0,.8,.1,.1)
        julia_assign("state", state)
        julia_assign("p", param)
        
        prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
        
        d_RNE = rbind(d_RNE, d[nrow(d), ] %>% add_column(S = S,alpha_0 = comp,Sp="Competitive"))
        d_RNE[d_RNE < 10^-4] = 0
        colnames(d_RNE) = c("rho_1", "rho_2", "rho_d", "rho_0", "S","alpha_0","Sp")
        
        #only stress tolerant one
        state =c(0.8,0,.1,.1)
        julia_assign("state", state)
        
        prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
        
        d_RNE = rbind(d_RNE, d[nrow(d), ] %>% add_column(S = S,alpha_0 = comp,Sp="Stress-tolerant"))
        d_RNE[d_RNE < 10^-4] = 0
        
        #both species
        state =c(.4,.4,.1,.1)
        julia_assign("state", state)
        prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
        
        d_RNE = rbind(d_RNE, d[nrow(d), ] %>% add_column(S = S, alpha_0 = comp,Sp="Both"))
      }
    }
    d_RNE$rho_plus = d_RNE$rho_1 + d_RNE$rho_2
    
    
    d_compe=filter(d_RNE,Sp=="Competitive")
    d_stesstol=filter(d_RNE,Sp=="Stress-tolerant")
    d_both=filter(d_RNE,Sp=="Both")
    
    d_RNE=tibble(S=d_both$S,alpha_0=d_both$alpha_0, 
                 RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1), #acutally its RII
                 RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2),
                 NintA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
                 NintA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)))
    d_RNE[,3:6][is.na(d_RNE[,3:6])] = NA

    
    
    #Plotting 
    
    
    d2$state = sapply(1:nrow(d2), function(x) {
      if (d2[x, 1] > 0 & d2[x, 2] > 0) {
        return("Coexistence")
      }
      if (d2[x, 1] > 0 & d2[x, 2] == 0) {
        return("Stress_tolerant")
      }
      if (d2[x, 1] == 0 & d2[x, 2] > 0) {
        return("Competitive")
      }
      if (d2[x, 1] == 0 & d2[x, 2] == 0) {
        return("Desert")
      }
    })
    

    
    color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
    
    # state at equilibrium
    p1 = ggplot(d2) +
      geom_tile(aes(x = S, y = alpha_0, fill = state)) +
      theme_classic() +
      scale_fill_manual(values = color_rho) +
      annotate("text", x = rep(1.05, 3), y = c_rne+max(c_rne)/20 , label = c("C", "B", "A"), color = "black") +
      theme(legend.position = "bottom") +
      labs(x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'), fill = "") +
      geom_hline(yintercept = c_rne , lwd = .1, color = "gray40") +
      theme(legend.text = element_text(size = 11))
    
    
    # density of global vegetation
    # density_col = colorRampPalette(c("red", "white", "blue"))
    # p2 = ggplot(d2) +
    #   geom_tile(aes(x = S, y = alpha_0 , fill = rho_plus)) +
    #   theme_classic() +
    #   scale_fill_gradientn(colours = density_col(100)) +
    #   theme(legend.position = "bottom") +
    #   labs(
    #     x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'),
    #     fill = TeX(r'(Global vegetation \ $\rho_1 + \rho_2$)')
    #   )
    
    
    # some bifurcation diagrams
    
    d_bifu = filter(d2, round(alpha_0, 4) %in% round(c_rne, 4))
    for (c_bifu in 1:3) {
      assign(
        paste0("p2_", c_bifu),
        ggplot(d_bifu %>% filter(., round(alpha_0, 4) == round(c_rne[c_bifu], 4)) %>% melt(., measure.vars = c("rho_1", "rho_2")) %>%
                 mutate(., variable = recode_factor(variable, "rho_1" = "Stress_tolerant", "rho_2" = "Competitive"))) +
          geom_point(aes(x = S, y = value, color = variable), size = .5) +
          geom_segment(
            x = S_critic2, y = .3+.075,xend = S_critic2, yend = .3,
            lineend = "round",linejoin = "round", size = .3, 
            arrow = arrow(length = unit(0.1, "inches")),
            colour = color_rho[2]) + 
          geom_segment(
            x = S_critic1, y = .3+.075,xend = S_critic1, yend = .3,
            lineend = "round",linejoin = "round", size = .3, 
            arrow = arrow(length = unit(0.1, "inches")),
            colour = color_rho[4]) + 
          labs(x = "Stress (S)", y = "Density", color = "") +
          the_theme +
          scale_color_manual(values = color_rho[c(2, 4)]) +
          theme(legend.text = element_text(size = 12))
      )
    }


    for (a0 in 1:3) {
      assign(
        paste0("p3_", a0),
        ggplot(d_RNE %>% filter(., round(alpha_0,4) == round(c_rne[a0],4)) %>%melt(., measure.vars=c("NintA_comp","NintA_st")) %>%
                 mutate(., variable = recode_factor(variable, "NintA_st" = "Stress-tolerant", "NintA_comp" = "Competitive"))) +
          geom_line(aes(x=S,y=value,color=variable))+the_theme+
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
          labs(x="Stress (S)",y="RII",color="")+ylim(-1,2)+
          scale_color_manual(values=c(as.character(color_rho)[c(4,2)]))
      )
    }

    
    pright_1=ggarrange(p2_3+theme(legend.position = "none")+xlab(""),p3_3+xlab("")+theme(legend.position = "none"),ncol=2,align = 'v',labels = letters[c(1,2)])
    pright_2=ggarrange(p2_2+theme(legend.position = "none")+xlab(""),p3_2+xlab("")+theme(legend.position = "none"),ncol=2,align = 'v',labels = letters[c(3,4)])
    pright_3=ggarrange(p2_1+theme(legend.text = element_text(size=12))+  guides(color = guide_legend(override.aes = list(size = 3))),
                       p3_1+theme(legend.text = element_text(size=12))+  guides(color = guide_legend(override.aes = list(size = 3))),
                       ncol=2,align = 'h',labels = letters[c(5,6)],common.legend = T,legend = "bottom")
    pright_tot=ggarrange(pright_1,pright_2,pright_3,nrow = 3,heights = c(1,1,1.2),align = "hv")
    
    
    
    p_tot = ggarrange(p1,pright_tot, ncol = 2, widths = c(2, 2))
    
    ggsave(paste0("../Figures/2_species/MF/2_species_MF",name_fig,".pdf"),width = 12,height = 6)
    
    
  }
}




# d_diversity=tibble()
# d2=read.table(paste0("../Table/2_species/2_species_main.csv"),sep=";")
# 
# for (k in unique(d2$alpha_0)[c(1, round(length(unique(d2$alpha_0)) / 2), length(unique(d2$alpha_0)))]){ #for each alpha_0
#   
#   d2_f=filter(d2,alpha_0==k)
#   for (i in 1:nrow(d2_f)){ #for each simulation along the stress gradient
#     
#     densi=d2_f[i,1:2]
#     if (any(densi==0)){
#       d_diversity=rbind(d_diversity,tibble(FD_0 =0, FD_1=0,  FD_2 =0,  D_0=0,   D_1=0,   D_2=0,)%>%add_column(., Stress=d2_f$S[i],alpha_0=k))
#       
#     } else{
#       d_diversity=rbind(d_diversity,Get_diversity_community(trait =c(1,0) ,densities =densi )%>%add_column(., Stress=d2_f$S[i],alpha_0=k))
#       
#     }
#     
#   }
# }
# 
# p=ggplot(d_diversity%>%melt(., id.vars=c("Stress","alpha_0")))+
#   geom_line(aes(x=Stress,y=value,color=variable))+
#   facet_grid(alpha_0~variable)+the_theme
# ggsave("../Figures/2_species/MF/Test_trait_diversity.pdf",p,width = 10,height = 4)


## 2) Hysteresis size as a function of competition coefficient ----


Run_dynamics_hysteresis(plot = T, N_seq_c = 3, N_seq_S = 300)



## 3) Net effects between species ----

# we want to calculate the net effect between both species. We evaluate that by increasing slightly the recruitment rate of 1 at eq

state = Get_MF_initial_state()
tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)



length_seq = 300;epsilon=10^(-5)
S_seq = seq(0, 1, length.out = length_seq)

for (func_MF in c("normal","SGH")[1]){

  d2 = d3 = tibble()
  param=Get_MF_parameters()  
  param=c(param[1:3],"beta1"=.8,"beta2"=.8,param[5:10])
  
  c_seq = seq(0, .3, length.out = length_seq)  #we constrain the competition coefficient explored to allow species to coexist.
  param["cintra"]=.3
  
  C_for_analyse = c_seq[c(1,round((3/4)*length(c_seq)),length(c_seq))]
  
  for (comp in C_for_analyse) {
    for (S in S_seq) {
      for (sp in 1:2) { #for both species
        
        param["S"] = S
        param["alpha_0"] = comp # update parameters
        julia_assign("p", param)
        
        
        # first without press perturbation
        
        
        if (func_MF=="normal"){
          prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
        } else {
          prob = julia_eval("ODEProblem(MF_two_species_SGH_press, state, tspan, p)")
        }
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
        colnames(mean_densities)="Densities"
        d2 = rbind(d2, mean_densities %>% add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp])) 
        d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, alpha_0 = comp))
        
        
        
        # with press perturbation on growth rate
        
        param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
        julia_assign("p", param)
        
        
        if (func_MF=="normal"){
          prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
        } else {
          prob = julia_eval("ODEProblem(MF_two_species_SGH_press, state, tspan, p)")
        }
        
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
        colnames(mean_densities)="Densities"
        d2 = rbind(d2, mean_densities %>% add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp]))
        
        param[paste0("beta",sp)] = .8
      }
    }
  }
  
  d2[d2 < epsilon] = 0
  colnames(d2) = c("eq", "S", "alpha_0", "Type","Species")
  
  
  net_effect =sapply(1:length(seq(1, nrow(d2) , by = 2)),function(x){
    i=seq(1, nrow(d2) , by = 2)[x]
    return((d2$eq[i+1] - d2$eq[i]) / epsilon)
  })
  
  d_net=d2%>%
    filter(., Type=="control")%>%
    select(.,-Type)
  d_net$value=net_effect
  
  #for latex format in facet
  appender <- function(string) {
    TeX(paste("$\\alpha_e = $", string))  
  }
  
  p=ggplot(d_net %>%
             mutate(., alpha_0=as.factor(round(alpha_0,2)))) +
    geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
    the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
    geom_hline(yintercept = 0,linetype=9)+
    facet_wrap(.~alpha_0)+
    labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial b_j}$)'))+
    theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank())
  
  
  
  colnames(d3)=c("rho_1","rho_2","rho_d","rho_0","S","alpha_0")
  
  p2=ggplot(d3%>%select(., -rho_d,-rho_0)%>%
              melt(., id.vars=c("S","alpha_0"))%>%
              mutate(., alpha_0=as.factor(round(alpha_0,2))))+
    
    geom_line(aes(x=S,y=value,color=as.factor(variable),alpha=alpha_0),lwd=1,alpha=.5)+
    the_theme+scale_color_manual(values=c("blue","green"),labels=c("rho_1"= "Stess_tol","rho_2"="Competitive"))+
    facet_wrap(.~alpha_0,labeller = as_labeller(appender, default = label_parsed))+
    labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y="Densities")
  
  
  p_tot=ggarrange(p2,p,common.legend = T,legend = "bottom",nrow=2,align = "v")
  ggsave(paste0("../Figures/2_species/MF/Net_effects_","main","_","test",".pdf"),p_tot,width = 8,height = 5)
  
}


## 4) Random initial communities ----

tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


n_random_ini=20
d2 = tibble()
S_seq = seq(0,1, length.out = 100)
c_seq=c(0,.15,.3)

for (S in S_seq) {
  for (comp in c_seq) {
    for (n_ini in 1:n_random_ini){ # for each combination, we draw random communities
      
      state =Get_MF_initial_state(type="random")
      julia_assign("state", state)
      param=Get_MF_parameters()
      param["S"] = S
      param["alpha_0"] = comp
      param["cintra"]=.3
      julia_assign("p", param)
      
      
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = comp,n_ini=n_ini))
    }
  }
}
d2[d2 < 10^-4] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2
colnames(d2) = c("Stress_tolerant", "Competitive", "rho_d", "rho_0", "Stress", "alpha_0","N_ini","Community")
set.seed(123)
v_vec=runif(2)
d2$CSI=d2[,1]*v_vec[1]+d2[,2]*v_vec[2]
write.table(d2,"../Table/2_species/Random_ini_com_MF.csv",sep=";")

d2=read.table("../Table/2_species/Random_ini_com_MF.csv",sep=";")

d2t=transform(d2,
              alpha_0 = factor(alpha_0, levels=c(0,.15,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.15","alpha[e] : 0.3")))


p1=ggplot(d2t%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
  geom_point(aes(x=Stress,y=value,color=variable),size=.3)+
  facet_wrap(.~alpha_0,labeller=labeller(alpha_0=label_parsed))+
  the_theme+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=as.character(color_rho)[c(4,2)],labels=c("Stress-tolerant","Competitive"))
p2=ggplot(d2t%>%melt(., measure.vars=c("CSI")))+
  geom_point(aes(x=Stress,y=value,color=variable),size=.3)+
  facet_grid(.~alpha_0,labeller=labeller(alpha_0=label_parsed))+
  the_theme+
  labs(x="Stress (S)",y="CSI",color="")+
  scale_color_manual(values="black")

p_tot=ggarrange(p1+xlab(""),p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),nrow=2,heights = c(1,.8))

ggsave("../Figures/2_species/MF/CSI_2_species.pdf",p_tot,width = 6,height = 6)




# Step 2) Pair approximation (PA) ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#

## 1) Exploration along competitive ability ----
### a) Do the simulation ----
dir.create("../Table/2_species/Sim_PA_scales",showWarnings = F)
# As we can get oscillations, function to get the mean densities

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



### b) Analyse + plot graphics ----

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
              labs(x="Stress (S)",y="RII",color="")+ylim(-1,2)+
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

### c) Clean fig for the interaction between dispersal, facilitation scales ----

N_rep = 50
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(0, .3, length.out = N_rep)
h_seq=c(1,1.5)[1]
f_seq=c(.9)

id_for_RNE = seq(1:N_rep)[c(1, length(seq(1:N_rep))/2, length(seq(1:N_rep)))]

name_scena=c("local_C_local_F","global_C_global_F",
             "local_C_global_F","global_C_local_F")

for (hcomp in h_seq){
  
  d_state=tibble()
  for (scena_ID in c(2,4)){ #for each scenario of species pairs
    
    for (f in f_seq){ #varying facilitation strength
      
      for (disp in c(.1, .9)){
        
        name_fig=paste0("_d_",disp,"_f_",f,"_h_",hcomp,"_scena_",name_scena[scena_ID])
        
        #loading data
        d_state=rbind(d_state,read.table(paste0("../Table/2_species/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")%>%
                        add_column(., Delta=disp,scena=name_scena[scena_ID]))
        
      }
    }
  }
  
  d_state=filter(d_state,alpha_0 == .3 )
  
  
  d_state=transform(d_state,
                  Delta = factor(Delta, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
                  scena=factor(scena,levels=c("global_C_global_F","global_C_local_F"),
                               labels=c("Global facilitation","Local facilitation")))
  
  

  p=ggplot(d_state %>% melt(.,measure.vars=c("rho_1","rho_2")))+
    geom_point(aes(x=S,y=value,color=variable),size=.7)+
    facet_grid(scena~Delta,labeller=labeller(Delta=label_parsed))+
    the_theme+
    scale_color_manual(values=rev(as.character(color_rho)[c(2,4)]),labels=c("rho_1"="Stress-tolerant","rho_2"="Competitive"))+
    labs(x="Stress (S)",y="Densities",color="")+  
    guides(color = guide_legend(override.aes = list(size = 3)))+
    theme(strip.text.x = element_text(size=12))
  
  ggsave(paste0("../Figures/2_species/PA/Mecanism_scale_facilitation_dispersal_h_",hcomp,".pdf"),width = 8,height = 5)
}



## 2) CSI and random initial communities ----

tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


n_random_ini=50
d2 = tibble()
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
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = .3,n_ini=n_ini,Scena=name_scena[scena_ID],Delta=disp))
      }
    }
  }
}


d2[d2 < 10^-4] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2
d2=d2[,c(1:2,10:14)]
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




## 3) State diagram multi-stability ----


tspan = c(0, 300)
t = seq(0, 300, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=5
c_seq=seq(0,.3,length.out=N_sim)
name_scena=c("global_C_local_F","global_C_global_F")
delta_seq=c(.1,.9)
branches=c("Degradation","Restoration")

for (branch in branches){

  if (branch =="Degradation"){
    state =Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
    S_seq=seq(0,1, length.out = N_sim)
  }else {
    state =Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
    S_seq=rev(seq(0,1, length.out = N_sim))
  }
  
  for (disp in delta_seq) {
    
    for (scena_ID in 1:2){ 
      
      for (ccomp in c_seq){
        
        for (S in S_seq) { 
          
          julia_assign("state", state)
          param=Get_PA_parameters()
          param["cintra"]=ccomp
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
d2=d2[,-10][,c(1:2,10:15)]
colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","Scena","Delta","Branches","Rho_plus")
write.table(d2,"../Table/2_species/Multistability_mapping.csv",sep=";")



#Analysis

d2t=transform(d2,
              Delta = factor(Delta, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
              Scena=factor(Scena,levels=c("global_C_global_F","global_C_local_F"),
                           labels=c("Global facilitation","Local facilitation")))


d2t$state = sapply(1:nrow(d2t), function(x) {
  if (is.na(d2t$rho_1[x]) & is.na(d2t$rho_2[x])) {
    return("NA")
  }else{
    
    
    if (d2t$rho_1[x] > 0 & d2t$rho_2[x] > 0) {
      return("Coexistence")
    }
    if (d2t$rho_1[x] > 0 & d2t$rho_2[x] == 0) {
      return("Stress-tolerant")
    }
    if (d2t$rho_1[x] == 0 & d2t$rho_2[x] > 0) {
      return("Competitive")
    }
    if (d2t$rho_1[x] == 0 & d2t$rho_2[x] == 0) {
      return("Desert")
    }
  }
  
})









## 4) Scale and strength of facilitation, dispersal and competition on niche expansion ----
### 1) Do the simulation ----



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
  d_niche=transform(d_niche,
                           Local_disp = factor(Local_disp, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")))
  p=ggplot(d_niche%>%
             melt(.,measure.vars=c("Delta_niche_2"))%>%
             filter(., Scena==scale_facil))+
    geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
    the_theme+
    facet_wrap(.~Local_disp,labeller=labeller(Local_disp=label_parsed))+
    scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
    labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="Niche expansion (%)")
  
  ggsave(paste0("../Figures/2_species/PA/Niche_expansion_PA_facilitation_",
                ifelse(scale_facil=="global_C_local_F","local","global"),".pdf"),p,width = 7,height = 4)
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


## 5) Clustering between species ----
### a) Simulation ----


tspan = c(0, 7000) #to avoid long transient
t = seq(0, 7000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 10
S_seq = c(0,.1,.4,.7)
alpha_seq = seq(0, .3, length.out = N_rep)
f_seq=seq(0,1,length.out=N_rep)
delta_seq=c(.1, .9)
h_seq=c(1,1.5)[1]


name_scena=c("local_C_local_F","global_C_global_F","local_C_global_F","global_C_local_F")

d_clustering=tibble() #initializing the tibble

for (scena_ID in 1:4){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (h in h_seq){ #varying competitive advantage strength
      
      for (f in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) { #varying competition
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          
          param["cintra"]=.3
          param["f"]=f
          param["h"]=h
          param["delta"]=disp
          param["alpha_0"]=alpha0
          
          
          #3) Coexisting species
          
          state=Get_PA_initial_state()
          julia_assign("state", state)
          
          #varying the global interspecific competition
          
          d2 = tibble()
          
          for (S in S_seq) { #varying the stress 
            
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
            
            d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, alpha_0=alpha0))

          }
          d2[d2 < 10^-4] = 0
          colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
          d2$rho_plus = d2$rho_1 + d2$rho_2
          
                    
          
          d_clustering=rbind(d_clustering,tibble(
            q12 = d2$rho_12 / d2$rho_2, q21 = d2$rho_12 / d2$rho_1, 
            q11 = d2$rho_11/d2$rho_1,q22 = d2$rho_22/d2$rho_2,
            c12 = d2$rho_12/(d2$rho_1*d2$rho_2),
            c11 = d2$rho_11/(d2$rho_1*d2$rho_1),
            c22 = d2$rho_22/(d2$rho_2*d2$rho_2),
            S   = d2$S,alpha_0 = d2$alpha_0,h=h,f=f,delta=disp,Scena=scena_ID
          ))

          
        } #end competition loop
        
      } #end facilitation loop
      
    } #end h loop
    
  } #end dispersal loop
  
} #end scenario loop

write.table(d_clustering,"../Table/2_species/Clustering_PA.csv",sep=";")



### b) Analysis ----
d_clustering = read.table("../Table/2_species/Clustering_PA.csv",sep=";")
fig_col=colorRampPalette(c("yellow","orange","red"))

for (i in 5:7){
  selected_col=d_clustering[,i]
  selected_col[which(selected_col==0)]=NA
  d_clustering[,i]=selected_col
}

for (hcomp in unique(d_clustering$h)){
  
  for (del in unique(d_clustering$delta)){
    
    
    d_fil=filter(d_clustering,h==as.numeric(hcomp),delta==as.numeric(del))%>%
      mutate(., Scena=recode_factor(Scena,"1"="LC_LF","2"="GC_GF","3"="LC_GF","4"="GC_LF"))
    
    
    
    for (i in c("c12","c11",'c22')){
      
      p=ggplot(d_fil%>%melt(measure.vars=i))+
        geom_tile(aes(x=f,y=alpha_0,fill=value))+
        facet_grid(Scena~S,scales = "free")+
        the_theme+labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill=i)+
        scale_fill_gradientn(colours=fig_col(10))
      
      
      ggsave(paste0("../Figures/2_species/PA/Clustering/Clustering_",i,"_h_",hcomp,"_dispersal_",del,".pdf"),width = 8,height = 6)
      
    }
    
    
  }
}


### c) Clustering with high intraspecific competition, high dispersal and low interspecific competition ----

tspan = c(0, 7000) #to avoid long transient
t = seq(0, 7000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 10
S_seq = c(0,.1,.4,.7)
alpha_seq = seq(0, .3, length.out = N_rep)
f_seq=seq(0,1,length.out=N_rep)
delta_seq=c(.1, .9)
h_seq=c(1,1.5)[1]
c_intra_seq=c(.3,.5,.7)


name_scena=c("local_C_local_F","global_C_global_F","local_C_global_F","global_C_local_F")

d_clustering=tibble() #initializing the tibble

for (scena_ID in c(1:4)[c(2,4)]){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (cii in c_intra_seq){ #varying cintra loop
      print(cii)
        
      for (facil in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) { #varying competition
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          
          param["cintra"]=cii
          param["f"]=facil
          param["h"]=1
          param["delta"]=disp
          param["alpha_0"]=alpha0
          
          
          #3) Coexisting species
          
          state=Get_PA_initial_state()
          julia_assign("state", state)
          
          #varying the global interspecific competition
          
          d2 = tibble()
          
          for (S in S_seq) { #varying the stress 
            
            param["S"] = S
            param["alpha_0"] = alpha0
            julia_assign("p", param)
            
            if (scena_ID==2){ #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }else { #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
              
            }
            

            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, alpha_0=alpha0))
            
          }
          d2[d2 < 10^-4] = 0
          colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
          d2$rho_plus = d2$rho_1 + d2$rho_2
          
          
          
          d_clustering=rbind(d_clustering,tibble(
            rho_1 = d2$rho_1, rho2 = d2$rho_2, rho12=d2$rho_12,
            q12 = d2$rho_12 / d2$rho_2, q21 = d2$rho_12 / d2$rho_1, 
            q11 = d2$rho_11/d2$rho_1,q22 = d2$rho_22/d2$rho_2,
            c12 = d2$rho_12/(d2$rho_1*d2$rho_2),
            c11 = d2$rho_11/(d2$rho_1*d2$rho_1),
            c22 = d2$rho_22/(d2$rho_2*d2$rho_2),
            S   = d2$S,alpha_0 = d2$alpha_0,cintra=cii,f=facil,delta=disp,Scena=scena_ID
          ))
          
          
        } #end competition loop
        
      } #end facilitation loop
      
    }
      
  } #end dispersal loop
  
} #end scenario loop

write.table(d_clustering,"../Table/2_species/Clustering_PA_intra_comp.csv",sep=";")


### d) Analysis ----
d_clustering = read.table("../Table/2_species/Clustering_PA_intra_comp.csv",sep=";")
fig_col=colorRampPalette(c("yellow","orange","red"))

for (i in 8:10){
  selected_col=d_clustering[,i]
  selected_col[which(selected_col==0)]=NA
  d_clustering[,i]=selected_col
}

for (cii in unique(d_clustering$cintra)){
  
  for (del in unique(d_clustering$delta)){
    
    
    d_fil=filter(d_clustering,cintra==as.numeric(cii),delta==as.numeric(del))%>%
      mutate(., Scena=recode_factor(Scena,"1"="LC_LF","2"="GC_GF","3"="LC_GF","4"="GC_LF"))
    
  
    for (i in c("c12","c11",'c22',"rho_1","rho2")[c(1)]){
      
      p=ggplot(d_fil%>%melt(measure.vars=i))+
        geom_tile(aes(x=f,y=alpha_0,fill=value))+
        facet_grid(Scena~S,scales = "free")+
        the_theme+labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill=i)+
        scale_fill_gradientn(colours=fig_col(10))
      
      
      ggsave(paste0("../Figures/2_species/PA/Clustering/Clustering_",i,"_cintra_",cii,"_dispersal_",del,".pdf"),width = 8,height = 6)
      
    }
    
    
  }
}


## 6) Net effects from partial derivative ----


state = Get_PA_initial_state()
tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)



length_seq = 100;epsilon=10^(-5)
S_seq = seq(0, 1, length.out = length_seq)
name_scena=c("global_C_local_F","global_C_global_F")
disp_seq=c(.1,.9)
c_seq = seq(0, .3, length.out = length_seq)


d_niche=d_RNE=d_all_dyn=tibble() 
d2 = d3 = tibble()
param=Get_PA_parameters()  
param=c(param[1:3],"beta1"=.8,"beta2"=.8,param[5:12])
param["cintra"]=.3
C_for_analyse = seq(0, .3, length.out = length_seq)[c(1,length_seq/2,length_seq)]

for (scena_ID in 1:2){ #for each scenario of species pairs
  
  for (comp in C_for_analyse) {
    
    for (disp in disp_seq){
      
      for (S in S_seq) {
        
        for (sp in 1:2) { #for both species
          
          param["S"] = S
          param["alpha_0"] = comp # update parameters
          param["delta"] = disp
  
          
          # first without press perturbation
          
          if (scena_ID==1){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F_press, state, tspan, p)")
            
          }
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          d2 = rbind(d2, as_tibble(d[nrow(d),c(1,2)[-sp]]) %>% add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],
                     Scena=name_scena[scena_ID],Disp=disp)) 
          d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, alpha_0 = comp))
          
          
          
          # with press perturbation on growth rate
          
          param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
          julia_assign("p", param)
          
          
          if (scena_ID==1){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F_press, state, tspan, p)")
            
          }
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          d2 = rbind(d2, as_tibble(d[nrow(d),c(1,2)[-sp]]) %>% add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],
                                                               Scena=name_scena[scena_ID],Disp=disp))
          
          param[paste0("beta",sp)] = 1
        }
      }
    }
  }
}

d2[d2 < epsilon] = 0
colnames(d2) = c("Eq", "S", "alpha_0", "Type","Species","Scena","Disp")


write.table(d2,"../Table/2_species/Net_effects_partial_deriv_PA.csv",sep=";")


#Analysis
net_effect =sapply(1:length(seq(1, nrow(d2) , by = 2)),function(x){
  i=seq(1, nrow(d2) , by = 2)[x]
  return((d2$Eq[i+1] - d2$Eq[i]) / epsilon)
})


d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect

d_net=transform(d_net,
                Disp = factor(Disp, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
                Scena=factor(Scena,levels=c("global_C_global_F","global_C_local_F"),
                             labels=c("Global facilitation","Local facilitation")))

for (id_plot in 1:3){
  assign(paste0("p_",id_plot),
         ggplot(d_net %>%filter(round(alpha_0,4)==round(C_for_analyse[id_plot]))) +
           geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
           the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
           geom_hline(yintercept = 0,linetype=9)+
           facet_grid(Scena~Disp,labeller=labeller(Disp=label_parsed))+
           labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial b_j}$)'))+
           theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank())
  )
  
  
  
}
p=


colnames(d3)=c("rho_1","rho_2","rho_d","rho_0","S","alpha_0")

p2=ggplot(d3%>%select(., -rho_d,-rho_0)%>%
            melt(., id.vars=c("S","alpha_0"))%>%
            mutate(., alpha_0=as.factor(round(alpha_0,2))))+
  
  geom_line(aes(x=S,y=value,color=as.factor(variable),alpha=alpha_0),lwd=1,alpha=.5)+
  the_theme+scale_color_manual(values=c("blue","green"),labels=c("rho_1"= "Stess_tol","rho_2"="Competitive"))+
  facet_wrap(.~alpha_0,labeller = as_labeller(appender, default = label_parsed))+
  labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y="Densities")


p_tot=ggarrange(p2,p,common.legend = T,legend = "bottom",nrow=2,align = "v")
ggsave(paste0("../Figures/2_species/MF/Net_effects_","main","_",func_MF,".pdf"),p_tot,width = 8,height = 5)


# Step 3) CA between both species : example and clustering ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
#

## 0) Comparing Gillespie simulations with classical IBM  ----
rm(list = ls())
source("./2_Species_analysis_functions.R")


t_seq = seq(1, 2000, 1)
param = c(
  r = 0.02, d = 0.1, f = .9, beta = 0.8, m = 0.05, e = 0, emax = 1,
  cintra = 0.1, cinter1 = .1, cinter2 = .1, S = 0, delta = .1, z = 4, leap = 1
)
Lattice_ini = Get_initial_lattice()

plot_dynamics(Run_CA_2_species(t_seq, param, ini = Lattice_ini)$d)
plot_dynamics(Gillespie_tau_leaping_R(Lattice_ini, t_seq, param)$state)


t_compare = tibble()
for (k in 1:10) {
  t1 = system.time(Run_CA_2_species(t_seq, param, ini = Lattice_ini))[1]
  t2 = system.time(Gillespie_tau_leaping_R(Lattice_ini, t_seq, param))[1]
  t_compare = rbind(t_compare, tibble(Comp_time = c(t1, t2), Method = c("Classic", "Gillespie")))
}
t_compare$Comp_time = as.numeric(t_compare$Comp_time)
t_compare %>%
  group_by(Method) %>%
  summarise(m_t = mean(Comp_time), s_t = sd(Comp_time))



d = tibble()
for (tau in c(1, 0.1, .01, .001)) {
  param["leap"] = tau
  d = rbind(d, (Gillespie_tau_leaping_R(Lattice_ini, t_seq, param)$state) %>% add_column(., tau_leap = tau))
}
color_rho = c("fertile" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")
d$time = rep(t_seq, 4)
p = ggplot(d %>% melt(., id.vars = c("time", "tau_leap")) %>%
             mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive", "rho_0" = "fertile", "rho_d" = "desert"))) +
  geom_line(aes(x = time, y = value, color = variable), lwd = 1) +
  theme_classic() +
  scale_color_manual(values = color_rho) +
  labs(x = "Time", y = "Densities", color = "") +
  theme(legend.text = element_text(size = 11), legend.position = "bottom") +
  facet_wrap(. ~ tau_leap, labeller = label_both)

ggsave("../Figures/2_species/Tau_values_dynamics_Gillespie.pdf", width = 7, height = 4)

## 1) An example along interspecific competition gradient ----
t_seq = seq(1, 3000, 1)
param = Get_CA_parameters()
Lattice_ini = Get_initial_lattice(size = 100)
test = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
plot_dynamics(test$state)
test2 = Gillespie_tau_leaping_R(time = t_seq, param = param, landscape = Lattice_ini)
plot_dynamics(test2$state)

d = tibble()
for (s in c(0, .25, .5, .75, 1)) {
  param["S"] = s
  Lattice_ini = Get_initial_lattice()
  CA = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
  d = rbind(d, as_tibble(melt(CA$landscape)) %>% add_column(., S = s))
}

color_CA = c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
p2 = ggplot(d) +
  geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
  theme_transparent() +
  scale_fill_manual(values = color_CA, labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
  facet_grid(. ~ S, labeller = label_both) +
  theme(panel.border = element_blank()) +
  theme(legend.position = "bottom") +
  labs(fill = "")


p_tot = ggarrange(p1, p2, ncol = 2, widths = c(1, 4))

ggsave("../Figures/2_species/CA/CA_example.pdf", p_tot, width = 16, height = 4, device = "pdf")



# We may be able to have coexistence at S=0 with no difference in competitive ability by changing delta
param = c(r = 0.02, d = 0.1, f = .9, beta = 0.8, m = 0.05, e = .1, emax = 1, cintra = 0.1, cinter1 = .1, cinter2 = .1, S = 0, delta = .1, z = 4, tau_leap = 1)

d = tibble()
for (delta in c(0.1, .25, .5, .75, 1)) {
  param["delta"] = delta
  Lattice_ini = Get_initial_lattice()
  CA = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
  d = rbind(d, as_tibble(melt(CA$landscape)) %>% add_column(., delta = delta))
}
p = ggplot(d) +
  geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
  theme_transparent() +
  scale_fill_manual(values = color_CA, labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
  facet_grid(. ~ delta, labeller = label_both) +
  theme(panel.border = element_blank()) +
  theme(legend.position = "bottom") +
  labs(fill = "")

ggsave("../Figures/2_species/CA/Varying_delta_CA.pdf", p, width = 13, height = 4, device = "pdf")


## 2) Clustering metrics along inter-specific competition gradient ----

# we use the Julia Code 2_species_clustering.jl
param = c(r = 0.02, d = 0.1, f = .9, beta = 0.8, m = 0.05, e = .1, emax = 1, cintra = 0.1, cinter1 = .3, cinter2 = .1, S = 0, delta = .1, dt = 1, z = 4)

d = read.table("../Table/2_species/Clustering_species_interspe_compet_gradient.csv", sep = ",")
d[is.na(d)] = 0 # for sake of simplicity
colnames(d) = c("rep", "q12", "c12", "cpp", "Relative_compet", "S")

# q12
p1 = d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    q12_m = median(q12), q12_q1 = quantile(q12, 0.25), q12_q3 = quantile(q12, 0.75), .groups = "keep"
  ) %>%
  ggplot(.) +
  geom_line(aes(x = Relative_compet / param["cinter2"], y = q12_m), color = "#F1B943", lwd = 1) +
  geom_ribbon(aes(x = Relative_compet / param["cinter2"], ymin = q12_q1, ymax = q12_q3), fill = "#F1B943", alpha = .5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ S, labeller = label_both) +
  labs(x = "", y = TeX("$\\q_{1|2}$"))

# c12
p2 = d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    c12_m = median(c12), c12_q1 = quantile(c12, 0.25), c12_q3 = quantile(c12, 0.75), .groups = "keep"
  ) %>%
  ggplot(.) +
  geom_line(aes(x = Relative_compet / param["cinter2"], y = c12_m), color = "#119680") +
  geom_ribbon(aes(x = Relative_compet / param["cinter2"], ymin = c12_q1, ymax = c12_q3), fill = "#119680", alpha = .5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ S) +
  labs(x = "", y = TeX("$\\c_{12}$")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


# cpp
p3 = d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    cpp_m = median(cpp), cpp_q1 = quantile(cpp, 0.25), cpp_q3 = quantile(cpp, 0.75), .groups = "keep"
  ) %>%
  ggplot(.) +
  geom_line(aes(x = Relative_compet / param["cinter2"], y = cpp_m), color = "#E63535") +
  geom_ribbon(aes(x = Relative_compet / param["cinter2"], ymin = cpp_q1, ymax = cpp_q3), fill = "#E63535", alpha = .5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ S) +
  labs(x = "Relative competition ability", y = TeX("$\\c_{++}$")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )



p_tot = ggarrange(p1, p2, p3, nrow = 3)
ggsave("../Figures/2_species/CA/Clustering_figure_interspe.pdf", p_tot, width = 8, height = 6)




## 3) Same but along stress gradient  ----

d = read.table("../Table/2_species/Clustering_species_stress_gradient.csv", sep = ",")
colnames(d) = c("rep", "q12",  "q21",  "q11" , "q22" , "c12",  "c21",  "c11",  "c22",  "qpp" , "cpp" , "stress")


d %>%
  melt(., measure.vars=c("q12",'q21',"c21"))%>%
  ggplot(.) +
  geom_line(aes(x = stress, y = value,color=variable),lwd=1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Stress (S)", y = "")


# # cpp
# p4 = d %>%
#     group_by(stress) %>%
#     summarise(
#         cpp_m = median(cpp), cpp_q1 = quantile(cpp, 0.25), cpp_q3 = quantile(cpp, 0.75), .groups = "keep"
#     ) %>%
#     ggplot(.) +
#     geom_line(aes(x = stress, y = cpp_m), color = "#E63535") +
#     geom_ribbon(aes(x = stress, ymin = cpp_q1, ymax = cpp_q3), fill = "#E63535", alpha = .5) +
#     theme_classic() +
#     theme(legend.position = "bottom") +
#     labs(x = "Relative competition ability", y = TeX("$\\c_{++}$"))


p_tot = ggarrange(p1, p2, p3, nrow = 3)
ggsave("../Figures/2_species/CA/Clustering_figure_stress_gradient.pdf", width = 5, height = 5)




## 4) Testing CA on specific parameters values ----

t_seq = seq(1, 5000, 1)
param = Get_CA_parameters()
param["alpha11"] = 0.4
param["alpha12"] = 0.01
param["alpha21"] = 0.1
param["alpha22"] = 0.01
param["cg"] = 0.01
param["S"] = 0.5
param["tau_leap"] = 1
Lattice_ini = Get_initial_lattice()

test = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
plot_dynamics(test$state)

test2 = Gillespie_tau_leaping_R(time = t_seq, param = param, landscape = Lattice_ini)

plot_dynamics(test2$state)
Plot_landscape(test2$landscape)

#
## 5) Patch size distribution ----

### a) Example  ----

mapping(100,100)
S_seq=seq(0,.6,length.out=4)
N_rep = 3 

d_freq=d_size=tibble()
display_landscape=tibble()
for (stress in S_seq){
  for (rep in 1:N_rep){
    print(stress)
    if (stress == 0){
      stress="0.0"
    }
    if (stress == 1){
      stress="1.0"
    }
    
    landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Example_sim/Landscape_size_100_stress_",stress,"_",rep,".csv"),sep=","))
    
    if (rep==1){
      display_landscape = rbind(display_landscape, as_tibble(melt(landscape)) %>% add_column(., S = stress))
      
    }
    
    Freq_patches=Get_frequency_number_patches(landscape)
    d_freq=rbind(d_freq,Freq_patches$Patches_frequency%>%add_column(., S=stress,Nrep=rep))
    d_size=rbind(d_size,Freq_patches$Patches_size%>%add_column(., S=stress,Nrep=rep))
  }
}



p1=ggplot(d_freq)+geom_point(aes(x=log(Size),log(Frequency),color=Species))+
  the_theme+labs(x="Patch size (k)",y="Frequency patch > k")+
  facet_wrap(.~S)+scale_color_manual(values=c("gray40",as.character(color_rho[c(4,2)])))

p2=ggplot(display_landscape) +
  geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
  theme_transparent() +
  scale_fill_manual(values =  c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
                    , labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
  theme(panel.border = element_blank()) +
  theme(legend.position = "bottom") +
  labs(fill = "")+
  facet_wrap(.~S)

p_tot=ggarrange(p1,p2,nrow=2,labels=LETTERS[1:2])
ggsave("../Figures/2_species/CA/Example_PSD.pdf",p_tot,width = 7,height = 10)



### b) Analysing the simulations for different competition regimes ----


# the simulations were made on Julia (Step 1 : Patch size distribution : different competition scenario)

mapping(100,100)
S_seq=seq(0,.6,length.out=4)
N_rep = 3 

for (scena in c("global","local")[2]){
  for (alpha0 in c("0.0", "0.25", "0.5")){
  
    d_freq=d_size=d_max_size=tibble()
    display_landscape=tibble()
  
    for (stress in S_seq){
      
      for (rep in 1:N_rep){
        
        print(stress)
        
        if (stress == 0){
          stress="0.0"
        }
        if (stress == 1){
          stress="1.0"
        }
        
        landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Competition_regime/Landscape_size_100_stress_",stress,"_",scena,"_",rep,"_relat_compet_",alpha0,".csv"),sep=","))
        if (rep==1 & length(unique(as.numeric(landscape)))>2){ #i.e. we have vegetation
          display_landscape = rbind(display_landscape, as_tibble(melt(landscape)) %>% add_column(., S = stress))
        }

        Freq_patches=Get_frequency_number_patches(landscape)
        d_freq=rbind(d_freq,Freq_patches$Patches_frequency%>%add_column(., S=stress,Nrep=rep))
        d_size=rbind(d_size,Freq_patches$Patches_size%>%add_column(., S=stress,Nrep=rep))
        d_max_size=rbind(d_max_size,d_size%>%group_by(Species,S)%>%summarise(max_patch=max(patch_size),.groups = "keep")) #max across replicates
      }
    }

    p1=ggplot(d_freq)+geom_point(aes(x=log(Size),log(Frequency),color=Species))+
      the_theme+labs(x="Patch size (k)",y="Frequency patch > k")+
      facet_wrap(.~S,labeller = label_both)+scale_color_manual(values=c("gray40",as.character(color_rho[c(4,2)])))

    p2=ggplot(display_landscape) +
      geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
      theme_transparent() +
      scale_fill_manual(values =  c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
                        , labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
      theme(panel.border = element_blank()) +
      theme(legend.position = "bottom") +
      labs(fill = "")+
      facet_grid(.~S,labeller = label_both)

    p_tot=ggarrange(p1,p2,nrow=2,labels=LETTERS[1:2],heights = c(1,.5))
    ggsave(paste0("../Figures/2_species/CA/Patch_size_distribution_",scena,"_Intersp_comp_",alpha0,".pdf"),p_tot,width = 7,height = 8)
  }
}



#



### c) Fitting Power law ----
mapping(100,100)
n_save = 25
alpha0 = c("0.0", "0.25", "0.5")

#doing the loop on simulations
for (scena in c("global","local")[-2]){
  
  if (scena=="global"){
    S_seq = c(seq(0, 0.8, length.out=30),seq(0.8,.9,length.out=4)[-1])
  } else{
    S_seq = seq(0, 0.8, length.out=30)
  }
  
  for (stress in S_seq){
    d_PL=tibble()

    stress=round(stress, 2)

    if (stress == 0){
      stress="0.0"
    }
    if (stress == 1){
      stress="1.0"
    }

    for (alpha_0 in alpha0){
      print(alpha_0)

      for (n in 1:n_save){


        #get landscape
        landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Big_sim/Landscape_",scena,"_S_",stress,
                                    "_alpha0_", alpha_0,"_nsave_", n , ".csv"),sep=","))



        #psd = patch size distribution
        Freq_patches=Get_frequency_number_patches(landscape)
        d_psd=Freq_patches$Patches_frequency
        colnames(d_psd) =  c("n","p","size","Species")

        #perform analysis for each species or for global vegetation
        for (sp in c("+","1","2")){

          best_model= NA
          psd_sp=filter(d_psd,Species==sp)
          # if (n==21 & alpha_0=="0.0" & stress==0.61 & scena==1) best_model=2

          # 1: check if desert for the focal species

          if(sum(landscape %in% c(1,2))/length(landscape) <0.01& sp=="+" |
             sum(landscape %in% c(1))/length(landscape) <0.01   & sp=="1"  |
             sum(landscape %in% c(2))/length(landscape) <0.01 & sp=="2" ){ best_model = 1 }

          # 2: check if vegetated for all vegetation
          if(is.na(best_model) & sum(landscape %in% c(2))/length(landscape) >= 0.70  & sp=="2" |
             is.na(best_model) & sum(landscape %in% c(1,2))/length(landscape) >= 0.70  & sp=="+" |
             is.na(best_model) & sum(landscape %in% c(1))/length(landscape) >= 0.70  & sp=="1")  {
            best_model = 5
          }

          # 3: fit power law models and compare via AIC via function

          if(is.na(best_model)) {


            p_spanning <- tail(psd_sp$p,1)

            result <- fitPL(psd_sp, p_spanning)

            best_model = result$best
          }

          class = c("DEG", "DOWN","PL", "UP", "COV")[best_model]

          if (class %in% c("DEG","COV")){alpha=NA}
          if (class %in% c("UP")){alpha=coefficients(result$TPLup)["alpha"]}
          if (class %in% c("DOWN")){alpha=coefficients(result$TPLdown)["alpha"]}
          if (class %in% c("PL")){alpha=coefficients(result$PL)["alpha"]}

          d_PL=rbind(d_PL,tibble(Class=class,Max_patch=max(psd_sp$size),Alpha_expo=alpha,N_rep=n,
                                 Species=sp,alpha_0=alpha_0,Stress=stress))

        } #end loop species

      } #end loop replicates

    } #end loop alpha0

    write.table(d_PL,paste0("../Table/2_species/PL_summary/PL_",scena,"_S_",stress,".csv"),sep=";")

  } # end loop stress

  
  
  #merging dataframes
  d_PL=tibble()
  for ( i in round(S_seq,2)){if (i==0) {i = "0.0"}
    d_PL=rbind(d_PL,read.table(paste0("../Table/2_species/PL_summary/PL_",scena,"_S_",i,".csv"),sep=";"))}
  
  d_PL$Max_patch[is.infinite(d_PL$Max_patch)]=0
  
  d_PL$Class[which(d_PL$Species=="2" & d_PL$a21==4 & d_PL$Class=="PL")]="DOWN"
  
  #for graphical purposes
  appender <- function(string) {
    TeX(paste("$\\alpha_e = $", string))}
  
  
  color_sp=c("+"="gray70","1"=as.character(color_rho[4]),"2"= as.character(color_rho[2]))
  color_class_PL =  c("COV"="#68B15E", "UP"= "#ECD57A",   "DOWN" ="#D22D2D","PL"="#E88D35",   "DEG"="#ADA9A9")
  
  #Max patch size colored by species
  d_PL_maxsize_sp=d_PL%>%
    group_by(., Species,alpha_0,Stress)%>%
    summarise(.groups = "keep",Max_patch=mean(Max_patch))
  
  
  p=ggplot(NULL)+
    geom_point(data=d_PL,aes(x=Stress,Max_patch,color=Species),size=.5,alpha=.3)+
    geom_line(data=d_PL_maxsize_sp,aes(x=Stress,Max_patch,color=Species),lwd=.9)+
    the_theme+labs(x="Stress (S)", y="Max patch size",color="")+
    facet_wrap(.~alpha_0,labeller = as_labeller(appender, default = label_parsed))+
    scale_color_manual(values=color_sp,
                       labels=c("+"="All vegetation","1"="Stress-tolerant","2"= "Competitive"))+
    scale_y_log10()
  
  ggsave(paste0("../Figures/2_species/CA/Max_patch_size_by_species_",scena,".pdf"),p,width = 7,height = 4)
  
  
  #Max patch size colored by PL class
  d_PL2=transform(d_PL,
                  alpha_0 = factor(alpha_0, levels=c(0,.25,.5), labels=c("alpha[e] : 0", "alpha[e] : 0.25",
                                                                         "alpha[e] : 0.5")),
                  Species=factor(Species,levels=c("+","1","2"),labels=c("Total vegetation","Stress-tolerant","Competitive")))
  
  
  d_PL_maxsize_class=d_PL%>%
    group_by(., Species,alpha_0,Stress,Class)%>%
    summarise(.groups = "keep",Max_patch=mean(Max_patch))
  
  d_PL_maxsize_class=transform(d_PL_maxsize_class,
                               alpha_0 = factor(alpha_0, levels=c(0,.25,.5), labels=c("alpha[e] : 0", "alpha[e] : 0.25",
                                                                                      "alpha[e] : 0.5")),
                               Species=factor(Species,levels=c("+","1","2"),labels=c("Total vegetation","Stress-tolerant","Competitive")))
  
  p2=ggplot(NULL)+
    geom_point(data=d_PL2,aes(x=Stress,Max_patch,color=Class),size=.5,alpha=.3)+
    geom_line(data=d_PL_maxsize_class,aes(x=Stress,Max_patch,color=Class),lwd=.9)+
    the_theme+labs(x="Stress (S)", y="Max patch size",color="")+
    facet_grid(Species~alpha_0,labeller=labeller(alpha_0=label_parsed))+
    scale_color_manual(values=color_class_PL,
                       labels=c("COV"="Covered","UP"="Up-bent PL","DOWN"= "Down-bent PL", "PL"="PL","DEG"="Degraded"))+
    scale_y_log10()
  
  ggsave(paste0("../Figures/2_species/CA/Max_patch_size_by_PL_class_",scena,".pdf"),p2,width = 7,height = 5)
  
  
  
  #Exponent PL by species
  
  d_PL_lambda_sp=d_PL%>%
    group_by(., Species,alpha_0,Stress)%>%
    summarise(.groups = "keep",Alpha_expo=mean(Alpha_expo))
  
  p3=ggplot(NULL)+
    geom_line(data=d_PL_lambda_sp,aes(x=Stress,-Alpha_expo,color=Species),lwd=.9)+
    geom_point(data=d_PL,aes(x=Stress,-Alpha_expo,color=Species),size=.5,alpha=.3)+
    the_theme+labs(x="Stress (S)",  y = TeX(r'(PL exponent \ \  $\lambda)'),color="")+
    facet_wrap(.~alpha_0,labeller = as_labeller(appender, default = label_parsed))+
    scale_color_manual(values=color_sp,
                       labels=c("+"="All vegetation","1"="Stress-tolerant","2"= "Competitive"))
  
  ggsave(paste0("../Figures/2_species/CA/PL_exponent_by_species_",scena,".pdf"),p3,width = 7,height = 4)
  
  
  #Exponent PL by PL class
  d_PL_lambda_class=d_PL%>%
    group_by(., Species,alpha_0,Stress,Class)%>%
    summarise(.groups = "keep",Alpha_expo=mean(Alpha_expo))
  
  d_PL_lambda_class=transform(d_PL_lambda_class,
                              alpha_0 = factor(alpha_0, levels=c(0,.25,.5), labels=c("alpha[e] : 0", "alpha[e] : 0.25",
                                                                                     "alpha[e] : 0.5")),
                              Species=factor(Species,levels=c("+","1","2"),labels=c("Total vegetation","Stress-tolerant","Competitive")))
  
  d_PL2=transform(d_PL,
                  alpha_0 = factor(alpha_0, levels=c(0,.25,.5), labels=c("alpha[e] : 0", "alpha[e] : 0.25",
                                                                         "alpha[e] : 0.5")),
                  Species=factor(Species,levels=c("+","1","2"),labels=c("Total vegetation","Stress-tolerant","Competitive")))
  p4=ggplot(NULL)+
    geom_line(data=d_PL_lambda_class,aes(x=Stress,-Alpha_expo,color=Class,group=interaction(Class,Species,alpha_0)),lwd=.9)+
    geom_point(data=d_PL2%>%mutate(., alpha_0=as.character(alpha_0)),aes(x=Stress,-Alpha_expo,color=Class),alpha=.3,size=.5)+
    the_theme+labs(x="Stress (S)", y = TeX(r'(PL exponent \ \  $\lambda)'),color="")+
    facet_grid(Species~alpha_0,labeller=labeller(alpha_0=label_parsed))+
    scale_color_manual(values=color_class_PL,
                       labels=c("COV"="Covered","UP"="Up-bent PL","DOWN"= "Down-bent PL", "PL"="PL","DEG"="Degraded"))
  
  ggsave(paste0("../Figures/2_species/CA/PL_exponent_by_PL_class_",scena,".pdf"),p4,width = 7,height = 5)
  
  
} #end loop scale competition

  
#
### d) Analyzing species dynamics ----

mapping(100,100)
n_save = 25
alpha_seq = c("0.0", "0.25", "0.5")
d_vege=tibble()
#doing the loop on simulations
for (scena in c("local","global")){
  
  if (scena=="global"){
    S_seq = c(seq(0, 0.8, length.out=30),seq(0.8,.9,length.out=4)[-1])
  } else{
    S_seq = seq(0, 0.8, length.out=30)
  }
  
  for (stress in S_seq){

    stress=round(stress, 2)
    
    if (stress == 0){
      stress="0.0"
    }
    if (stress == 1){
      stress="1.0"
    }
    
    for (alpha0 in alpha_seq){
      
      for (n in 1:n_save){
        
        
        #get landscape
        landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Big_sim/Landscape_",scena,"_S_",stress,
                                              "_alpha0_", alpha0,"_nsave_", n , ".csv"),sep=","))
        
        d_vege=rbind(d_vege,tibble(Rho_1=sum(landscape==1)/length(landscape),
                               Rho_2=sum(landscape==2)/length(landscape),Rho_p=sum(landscape %in% c(1,2))/length(landscape),N_rep=n,
                               alpha_0=alpha0,Stress=as.numeric(stress)))
      } #end loop replicates
    } #end loop a12
  } # end loop stress
  
  
  
  d_vege_merged=d_vege%>%
    group_by(.,alpha_0,Stress)%>%
    summarise(.,Rho_1=mean(Rho_1),Rho_2=mean(Rho_2),Rho_p=mean(Rho_p),.groups = "keep" )
  
  color_sp=c("Rho_p"="gray70","Rho_1"=as.character(color_rho[4]),"Rho_2"= as.character(color_rho[2]))
  appender <- function(string) {
    TeX(paste("$\\alpha_e = $", string))}
  
  
  
  p=ggplot(NULL)+
    geom_path(data=d_vege_merged%>% melt(., id.vars=c("alpha_0","Stress")),
              aes(x=Stress,value,color=variable,group=variable),lwd=1)+
    geom_path(data=d_vege%>%melt(., id.vars=c("alpha_0","Stress","N_rep")),
              aes(x=Stress,value,color=variable,group=interaction(variable,N_rep)),alpha=.17)+
    the_theme+labs(x="Stress (S)", y="Patch density",color="")+
    facet_wrap(.~alpha_0,labeller = as_labeller(appender, default = label_parsed))+
    scale_color_manual(values=color_sp,
                       labels=c("Rho_p"="All vegetation","Rho_1"="Stress-toletant","Rho_2"= "Competitive"))
  
  ggsave(paste0("../Figures/2_species/CA/Species_dynamics_CA_",scena,".pdf"),p,width = 7,height = 4)
  
  
} #end loop scenarios competition




#
# Step 4) Testing EWS on spatial dynamics ----
rm(list = ls())
source("./2_Species_analysis_functions.R")

# Here we test the spatial EWS on the whole dynamics. We used the previous simulations of the PL exponent
# For each replicate, we compute the slope of variance, skewness and near-neighbor correlation (Moran I) along stress gradient


mapping(100,100)
n_save = 25
alpha_seq = c("0.0", "0.25", "0.5")


Spatial_EWS=tibble()
for (scena in c("local","global")){
  
  if (scena=="global"){
    S_seq = c(seq(0, 0.8, length.out=30),seq(0.8,.9,length.out=4)[-1])
  } else{
    S_seq = seq(0, 0.8, length.out=30)
  }
  
  for (n in 1:n_save){

    for (alpha0 in alpha_seq){
      
      for (sp in c(1,2,"+")){
        
        all_landscape=list()
        u=1
        for (stress in S_seq){
          
          
          stress=round(stress, 2)
          
          if (stress == 0){
            stress="0.0"
          }
          if (stress == 1){
            stress="1.0"
          }
          
          #get landscape
          landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Big_sim/Landscape_",scena,"_S_",stress,
                                                "_alpha0_", alpha0,"_nsave_", n , ".csv"),sep=","))
          
          if (sp=="+"){
            landscape[landscape %in% c(1,2)]=1
            landscape[landscape<1]=0
            landscape=landscape>0
          } else {
            landscape[landscape != as.numeric(sp)]=0
            landscape[landscape == as.numeric(sp)]=1
            landscape=landscape>0
          }
          
          all_landscape[[u]] = landscape # putting all landscapes in a list in order to use spatial warnings package
          
          u=u+1
        } #end loop stress
        
        # Getting the spatial EWS using spatialwarning package
        generic_sp=as_tibble(as.data.frame(generic_sews(all_landscape,subsize = 5,moranI_coarse_grain = T)))%>%
          mutate(.,matrixn=rep(S_seq,each=4))%>%
          rename(., Stress=matrixn)%>%
          add_column(., alpha_0=as.numeric(alpha0),
                     replicate=n,Species=sp)
        
        fit_spatial_ews=tibble()
        
        # Applying linear regression along stress gradient and saving the slope for each measure
        for (measure in unique(generic_sp$indic)[-4]){ #deleting the mean
          
          generic_sp_measure=filter(generic_sp,indic==measure)
          slope_measure=lm(value~Stress,data=generic_sp_measure)$coefficients[2] #getting the slope
          
          test_slope = summary(lm(value~Stress,data=generic_sp_measure))$coefficients[2,4]
          
          fit_spatial_ews=rbind(fit_spatial_ews,tibble(Slope=slope_measure,Metric=measure, #merging data
                                                       alpha_0=unique(generic_sp$alpha_0),Rep=unique(generic_sp$replicate),
                                                       Species=sp,Signif = test_slope))
          
        }
        
        Spatial_EWS=rbind(Spatial_EWS,fit_spatial_ews) #merging to the main df
        
      }#end species loop
      
    } #end loop a12
  
  } # end loop replicate
  write.table(Spatial_EWS,paste0("../Table/2_species/EWS_spatial/Spatial_EWS_",scena,".csv"),sep=";")
  
  
  #Analysing data and ploting graphs
  
  
  Spatial_EWS=as_tibble(mutate(Spatial_EWS,Metric=recode_factor(Metric,"skewness"="Skewness","moran"="Moran's I","variance"="Variance"),
                               Species=as.character(Species),
                               alpha_0=as.numeric(alpha_0),
                               Signif=as.numeric(Signif),
                               Metric=as.character(Metric)))
  
  Spatial_EWS_replicate=Spatial_EWS%>%
    group_by(Species,alpha_0,Metric)%>%
    summarise(.groups = "keep",Mean_slope=mean(Slope),Sd_slope=sd(Slope))
  
  
  
  color_sp=c("gray70",as.character(color_rho[4]),as.character(color_rho[2]))
  
  p=ggplot(NULL)+
    geom_jitter(data=Spatial_EWS,aes(x=as.numeric(alpha_0),y=Slope,color=Species,fill=Species),alpha=.5,width = 0.025,height = 0)+
    geom_errorbar(data=Spatial_EWS_replicate,aes(x=as.numeric(alpha_0),y=Mean_slope,
                                                 ymin=Mean_slope-Sd_slope,
                                                 ymax=Mean_slope+Sd_slope,fill=Species), shape=21,width=0,color="black") +
    geom_point(data=Spatial_EWS_replicate,aes(x=as.numeric(alpha_0),y=Mean_slope,fill=Species),color="black",shape=21,size=3)+
    scale_color_manual(values=color_sp,
                       labels=c("+"="Total vegetation",
                                "1"="Stress-tolerant",
                                "2"="Competitive"))+
    scale_fill_manual(values=color_sp,
                      labels=c("+"="Total vegetation",
                               "1"="Stress-tolerant",
                               "2"="Competitive"))+
    facet_wrap(.~Metric,scales = "free")+the_theme+
    labs(y="Slope",x=TeX(r'(Competition strength \ $\alpha_e)'),color="",fill="")+
    geom_hline(yintercept = 0,linetype=9,lwd=.1)+
    scale_x_continuous(breaks = c(0,.25,.5))
  
  ggsave(paste0("../Figures/2_species/CA/Spatial_EWS_",scena,".pdf"),width = 8,height = 4)
  
  
} #end loop scenarios competition







# Step 5) Testing Temporal EWS ----

rm(list = ls())
source("./2_Species_analysis_functions.R")

# Here we test the spatial EWS on the whole dynamics. We used the previous simulations of the PL exponent
# For each replicate, we compute the slope of variance, skewness and near-neighbor correlation (Moran I) along stress gradient


mapping(100,100)
n_save = 25
alpha_seq = c("0.0", "0.25", "0.5")


Temporal_EWS=tibble()
for (scena in c("local","global")){
  
  if (scena=="global"){
    S_seq = c(seq(0, 0.8, length.out=30),seq(0.8,.9,length.out=4)[-1])
  } else{
    S_seq = seq(0, 0.8, length.out=30)
  }
  
  for (alpha0 in alpha_seq){
    
    d_dyn=tibble()
    
    all_landscape=list()
    u=1
    for (stress in S_seq){
      stress=round(stress, 2)
      
      if (stress == 0){
        stress="0.0"
      }
      if (stress == 1){
        stress="1.0"
      }
      
      for (n in 1:n_save){
        for (sp in c(1,2,"+")){
          
          
          
          #get landscape
          landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Big_sim/Landscape_",scena,"_S_",stress,
                                                "_alpha0_", alpha0,"_nsave_", n , ".csv"),sep=","))
          
          if (sp=="+"){
            landscape[landscape %in% c(1,2)]=1
            landscape[landscape<1]=0
            landscape=landscape>0
          } else {
            landscape[landscape != as.numeric(sp)]=0
            landscape[landscape == as.numeric(sp)]=1
            landscape=landscape>0
          }
          
          d_dyn=rbind(d_dyn,tibble(Sp=c(as.character(sp)),
                                   Rho=length(landscape[landscape==T])/length(landscape),
                                   Rep=n,Stress=stress,
                                   alpha_0=alpha0))
          
          all_landscape[[u]] = landscape # putting all landscapes in a list in order to use spatial warnings package
          
          u=u+1
        } #end species loop
        
      }# end loop replicate 
      
      
    }#end loop stress
    
    
    t=d_dyn%>%
      group_by(Sp,Stress,alpha_0)%>%
      summarise(Rho_m=as.numeric(mean(Rho)),.groups = "keep")%>%
      mutate(Rho_m=as.numeric(Rho_m),Stress=as.numeric(Stress))%>%
      as.data.frame(.)
      
    ggplot(t)+
      geom_path(aes(x=Stress,y=Rho_m,color=Sp))
      
    
    for (a0 in unique(d_dyn$alpha_0)){
        for (stress in unique(d_dyn$Stress)){
          d_fil=filter(d_dyn,alpha_0==a0,Stress==stress)%>%
            spread(Sp,Rho)
          
          
          
        }
    }
    
    
    
  } #end loop a12
  
}


