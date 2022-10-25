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




## 2) Hysteresis size as a function of competition coefficient ----


Run_dynamics_hysteresis(plot = T, N_seq_c = 3, N_seq_S = 300)



## 3) Net effects between species ----

# we want to calculate the net effect between both species. We evaluate that by increasing slightly the recruitment rate of 1 at eq


tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim = 100;epsilon=10^(-8)
C_for_analyse = seq(0, .3, length.out = N_sim)[c(1,round((1/2)*N_sim),N_sim)]
S_seq = seq(0, 1, length.out = N_sim)
branch_seq=c("Restoration","Degradation")
param=Get_MF_parameters()  
param=c(param[1:3],"beta1"=1,"beta2"=1,param[5:10])

d2 = d3 = tibble()

for (branch in branch_seq){ #branch bifurcation diagrams
  
  if (branch =="Degradation"){ 
    state = Get_MF_initial_state(c(.4,.4,.1))
    S_seq = seq(0, 1, length.out = N_sim)
  }else {
    state = Get_MF_initial_state(c(.005,.005,.49))
    S_seq = rev(seq(0, 1, length.out = N_sim))
  }
  julia_assign("state", state)
  
  for (comp in C_for_analyse) {
    
    for (S in S_seq) {
      
      for (sp in 1:2) { #for both species
        
        param["S"] = S
        param["alpha_0"] = comp # update parameters
        julia_assign("p", param)
        param["cintra"]=.3
        
        
        # first without press perturbation
        
        prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
        colnames(mean_densities)="Densities"
        
        d2 = rbind(d2, mean_densities %>% add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],Branches=branch)) 
        d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, alpha_0 = comp,Branches=branch))
        
        
        # with press perturbation on growth rate proxy
        param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
        julia_assign("p", param)
        
        prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
        colnames(mean_densities)="Densities"
        
        d2 = rbind(d2, mean_densities %>% add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],Branches=branch))
        
        param[paste0("beta",sp)] = 1
      }
    }
  }
  
}

d2[d2 < epsilon] = 0
colnames(d2) = c("eq", "S", "alpha_0", "Type","Species","Branches")

net_effect =sapply(seq(1, nrow(d2) , by = 2),function(x){
  return((d2$Eq[x+1] - d2$Eq[x]) / epsilon)
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
           mutate(., alpha_0=as.factor(round(alpha_0,2)))%>%
           filter(., Branches=="Degradation")) +
  geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
  the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess-tolerant","2"="Competitive"))+
  geom_hline(yintercept = 0,linetype=9)+
  facet_wrap(.~alpha_0)+
  labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial \beta_j}$)'))

ggsave(paste0("../Figures/2_species/MF/Net_effects_MF.pdf"),p,width = 7,height = 4)

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
write.table(d2,"../Table/2_species/MF/Random_ini_com_MF.csv",sep=";")

d2=read.table("../Table/2_species/MF/Random_ini_com_MF.csv",sep=";")

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




## 5) Multistability ----



tspan = c(0, 5000)
t = seq(0, 5000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.3,length.out=N_sim)
branches=c("Degradation","Restoration")

for (branch in branches){
  
  if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
    state =Get_MF_initial_state(c(.4,.4,.1))
    S_seq=seq(0,1, length.out = N_sim)
  }else {
    state =Get_MF_initial_state(c(.005,.005,.49))
    S_seq=rev(seq(0,1, length.out = N_sim))
  }
  
  
  for (ccomp in c_seq){
    
    for (S in S_seq) { #varying dispersal scale
      
      julia_assign("state", state)
      param=Get_MF_parameters()
      param["cintra"]=.3
      param["S"] = S
      param["alpha_0"] = ccomp
      julia_assign("p", param)
      
      
      
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      
      
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = ccomp,
                                                 Branch=branch))
      
    }
  }
}
d2[d2 < 10^-4] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2
d2=d2[,-c(3:4)]
colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","Branches","Rho_plus")
write.table(d2,paste0("../Table/2_species/MF/Multistability_MF.csv"),sep=";")



d2=read.table(paste0("../Table/2_species/MF/Multistability_MF.csv"),sep=";")

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

d2=d2[order(d2$alpha_0,d2$Stress),]

all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
  if (d2$state[x] != d2$state[x+1]){
    return(paste0(d2$state[x],"/", d2$state[x+1]))
  }
  else {return(d2$state[x])}
})

d_state=d2%>%
  filter(., Branches=="Degradation")%>%
  select(.,-Branches)
d_state$all_state=all_state


color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")

p=ggplot(d_state) +
  geom_tile(aes(x=Stress,y=alpha_0,fill=all_state))+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Stress_tolerant" = "#7BD8D3",
                             "Stress_tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress_tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  geom_hline(yintercept = round(unique(d2$alpha_0)[c(1,50,100)],4))+
  the_theme+theme(legend.text = element_text(size=8))


for (u in 1:3){
  d_plot=filter(d2,round(alpha_0,4) == round(unique(d2$alpha_0)[c(1,50,100)[u]],4))
  assign(paste0("p_",u),
         ggplot(d_plot%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
           geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
           labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
           the_theme +
           scale_color_manual(values = color_rho[c(2, 4)]) +
           scale_linetype_manual(values=c(1,9))+
           theme(legend.text = element_text(size = 7))
         )
  
}

p2=ggarrange(p_3+xlab(""),p_2+xlab(""),p_1,nrow = 3,common.legend = T,legend = "bottom")


ggsave("../Figures/2_species/MF/Multistability.pdf",ggarrange(p,p2,ncol = 2,widths = c(1.5,1)),width = 10,height = 6)

## 6) Niche expansion ----



tspan = c(0, 7000) #to avoid long transient
t = seq(0, 7000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 100
N_rep2 = 10
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(0, .3, length.out = N_rep2)
f_seq=seq(0,1,length.out=N_rep2)

d_niche=d_RNE=d_all_dyn=tibble() #initializing the tibble

for (f in f_seq){ #varying facilitation strength
  
  for (alpha0 in alpha_seq) {
    
    
    #Setting the parameters
    param=Get_MF_parameters()
    
    param["cintra"]=.3
    param["f"]=f
    param["alpha_0"]=alpha0
    
    
    
    # 1) Competitive species alone
    
    state=Get_MF_initial_state(c(0,.8,.1))
    julia_assign("state", state)
    d2 = d3 = tibble()
    
    for (S in S_seq) {
      
      param["S"] = S
      

      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      

      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Competitive"))
      d3 = rbind(d3, d[nrow(d), ] %>% add_column(S = S,Sp="Competitive",alpha_0=alpha0))
      
    }
    d2[d2 < 10^-4] = 0
    d3[d3 < 10^-4] = 0
    
    S_critic2_alone=abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) #range of values where competitive species is
    
    
    
    # 2) Stress-tolerant species alone
    
    state=Get_MF_initial_state(c(0.8,0,.1))
    julia_assign("state", state)
    
    d2 = tibble()
    S_seq = seq(0, 1, length.out = 100)
    
    for (S in S_seq) {
      
      param["S"] = S
      julia_assign("p", param)
      
      
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      
      
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Stress_tolerant"))
      d3 = rbind(d3, d[nrow(d), ] %>% add_column(S = S,Sp="Stress_tolerant",alpha_0=alpha0))
      
    }
    
    d2[d2 < 10^-4] = 0
    d3[d3 < 10^-4] = 0
    
    S_critic1_alone=abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
    
    
    
    #3) Coexisting species
    
    state=Get_MF_initial_state()
    julia_assign("state", state)
    
    #varying the global interspecific competition
    
    d2 = tibble()
    
    for (S in S_seq) { #varying the stress 
      
      param["S"] = S
      param["alpha_0"] = alpha0
      julia_assign("p", param)
      
      
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      
      
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
      
      d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, Sp="Both"))
      d3 = rbind(d3, d[nrow(d), ] %>% add_column(S = S,Sp="Both",alpha_0=alpha0))
      
    }
    d2[d2 < 10^-4] = 0
    d3[d3 < 10^-4] = 0
    colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_0", "S", "alpha_0")
    d2$rho_plus = d2$rho_1 + d2$rho_2
    
    S_critic1_both = abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
    S_critic2_both = abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) 
    
    
    
    
    #Putting niche in the tibble
    
    d_niche=rbind(d_niche,tibble(
      Facilitation = f,  alpha_0 =alpha0,
      Delta_niche_1 = 100 * (S_critic1_both - S_critic1_alone) / S_critic1_alone,  #making it a percentage of initial niche
      Delta_niche_2 = 100 * (S_critic2_both - S_critic2_alone) / S_critic2_alone)) #making it a percentage of initial niche
    
    
    
    d3$rho_plus = d3$rho_1 + d3$rho_2
    
    
    d_compe=filter(d3,Sp=="Competitive")
    d_stesstol=filter(d3,Sp=="Stress_tolerant")
    d_both=filter(d3,Sp=="Both")
    
    d_RNE=rbind(d_RNE,tibble(S=d_both$S,alpha_0=d_both$alpha_0, Facilitation=f, 
                             RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1),
                             RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2)))
    
    
    d_all_dyn=rbind(d3,d_all_dyn)
    
    
  } #end competition loop
  
} #end facilitation loop

write.table(d_niche,"../Table/2_species/MF/Niche_expansion_MF.csv",sep=";")





d_niche=read.table("../Table/2_species/MF/Niche_expansion_MF.csv",sep=";")
alpha_seq=seq(0,.3,length.out=10)
d_niche$alpha_0=alpha_seq


p=ggplot(d_niche%>%
           melt(.,measure.vars=c("Delta_niche_2")))+
  geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
  the_theme+
  scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
  labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="Niche expansion (%)")

ggsave(paste0("../Figures/2_species/MF/Niche_expansion_MF_facilitation.pdf"),p,width = 4,height = 4)





## 7) Comparing Net-effects and RII ----
### a) Simulation ----

tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
Nsim=100
S_seq = seq(0, 1, length.out = Nsim)
c_seq = seq(0,.3, length.out = Nsim)
epsilon=10^(-9)

d_RNE=d_Net_effect=tibble()


for (S in S_seq) {
  
  for (comp in c_seq) {
    
    # 1) RII 
    
    param=Get_MF_parameters()
    param["S"] = S
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
    
    # 2) Net-effects 
    
    param=Get_MF_parameters()  
    param=c(param[1:3],"beta1"=1,"beta2"=1,param[5:10])
    
    for (sp in 1:2) { #for both species
      
      param["S"] = S
      param["alpha_0"] = comp # update parameters
      julia_assign("p", param)

      
      # first without press perturbation
      
      prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
      colnames(mean_densities)="Densities"
      
      d_Net_effect = rbind(d_Net_effect, mean_densities %>% add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp])) 

      
      # with press perturbation on growth rate proxy
      param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
      julia_assign("p", param)
      
      prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
      colnames(mean_densities)="Densities"
      
      d_Net_effect = rbind(d_Net_effect, mean_densities %>% add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp]))
      
      param[paste0("beta",sp)] = 1
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
save(file="../Table/2_species/MF/Comparing_net_effects.RData",d_all)

### b) Analysis ----
load(file="../Table/2_species/MF/Comparing_net_effects.RData")
d_net=d_all[[1]];d_RNE=d_all[[2]]


d_net$value[d_net$value>10]=NA
d_net$value[d_net$value < -10]=NA

ggplot(d_net)+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_wrap(.~Species)+
  the_theme+
  scale_fill_gradient2(low="red",mid="white",high="blue")


p1=ggplot(d_RNE%>%
         melt(., measure.vars=c("RNE_stress_tol","RNE_competitive"))%>%
         mutate(., variable=recode_factor(variable,"RNE_stress_tol"='Stress-tolerant',"RNE_competitive"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_wrap(.~variable)+
  the_theme+labs(x="",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=12))

p2=ggplot(d_RNE%>%
         melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
         #mutate(., value=rescale(value,to=c(-1,1)))%>%
         mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_wrap(.~variable)+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  theme(strip.background.x = element_blank(),strip.text.x = element_blank())

p_tot=ggarrange(p1,p2,nrow = 2,common.legend = T,legend = "bottom",labels=LETTERS[1:2])
ggsave("../Figures/2_species/MF/Comparizon_RII_NIntA.pdf",p_tot,width = 7,height = 7)



# NIntA

p2=ggplot(d_RNE%>%
            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
            #mutate(., value=rescale(value,to=c(-1,1)))%>%
            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_wrap(.~variable)+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .1))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .32,x=1.05),label="a")+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .3))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .12,x=1.05),label="b")+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")

alpha_for_plot=c(.1, .3)

for (alpha in 1:2){
  assign(paste0("p_",alpha),ggplot(d_RNE%>%
                                     filter(., alpha_0==alpha_for_plot[alpha])%>%
                                     melt(., measure.vars=c("NintA_comp")))+#%>%
           #           mutate(., value=rescale(value,to=c(-1,1))))+
           geom_hline(yintercept = 0,linetype=9,lwd=.5)+
           geom_line(aes(x=S,y=value))+
           the_theme+theme(axis.text = element_text(size =9), axis.title = element_text(size = 10))+
           labs(x="Stress (S)",y=TeX("$\\NInt_{A}$")))
}

p_right=ggarrange(p_2+theme(axis.text.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank())+
                    labs(x=""),p_1,ggplot() + theme_void(),nrow=3,labels=letters[1:2])

p_tot=ggarrange(p2,p_right,widths = c(2.2,1))
ggsave("../Figures/2_species/MF/NIntA_net_effect.pdf",p_tot,width = 7,height = 4)


# RII

p2=ggplot(d_RNE%>%
            melt(., measure.vars=c("RNE_stress_tol","RNE_competitive"))%>%
            mutate(., variable=recode_factor(variable,"RNE_stress_tol"='Stress-tolerant',"RNE_competitive"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_wrap(.~variable)+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .1))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .32,x=1.05),label="a")+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .3))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .12,x=1.05),label="b")+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")







for (alpha in 1:2){
  assign(paste0("p_",alpha),ggplot(d_RNE%>%
                                     filter(., alpha_0==alpha_for_plot[alpha])%>%
                                     melt(., measure.vars=c("RNE_competitive")))+#%>%
           #           mutate(., value=rescale(value,to=c(-1,1))))+
           geom_hline(yintercept = 0,linetype=9,lwd=.5)+
           geom_line(aes(x=S,y=value))+
           the_theme+theme(axis.text = element_text(size =9), axis.title = element_text(size = 10))+
           labs(x="Stress (S)",y=TeX("$\\NInt_{A}$")))
}

p_right=ggarrange(p_2+theme(axis.text.x = element_blank(),
                           axis.line.x = element_blank(),
                           axis.ticks.x = element_blank())+
                   labs(x=""),p_1,ggplot() + theme_void(),nrow=3,labels=letters[1:2])

p_tot=ggarrange(p2,p_right,widths = c(2.2,1))
ggsave("../Figures/2_species/MF/RII_net_effect.pdf",p_tot,width = 7,height = 4)




# Step 2) Pair approximation (PA) ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#

## 1) Exploration along competitive ability ----
### a) Do the simulation ----
dir.create("../Table/2_species/PA/Sim_PA_scales",showWarnings = F)
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
                  paste0("../Table/2_species/PA/Sim_PA_scales/Tipping_sp",name_fig,".csv"),sep=";")
      
      
      
      
      
      
      
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
      
      write.table(d2,paste0("../Table/2_species/PA/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")
      
      
      
      
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
      
      write.table(d_RNE,paste0("../Table/2_species/PA/Sim_PA_scales/2_species_RNE",name_fig,".csv"),sep=";")
      
      
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
        d_state=read.table(paste0("../Table/2_species/PA/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")
        Scritic=read.table(paste0("../Table/2_species/PA/Sim_PA_scales/Tipping_sp",name_fig,".csv"),sep=";")
        d_RNE=read.table(paste0("../Table/2_species/PA/Sim_PA_scales/2_species_RNE",name_fig,".csv"),sep=";")
  
        
        
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
        d_state=rbind(d_state,read.table(paste0("../Table/2_species/PA/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")%>%
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
write.table(d2,"../Table/2_species/PA/Random_ini_com_PA.csv",sep=";")


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




## 3) State diagram multistability ----
### a) Simulation ----




tspan = c(0, 8000)
t = seq(0, 8000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.3,length.out=N_sim)
name_scena=c("global_C_local_F","global_C_global_F","local_C_local_F","local_C_global_F")
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
    
    for (scena_ID in 1:4){ #for each scenario of species pairs
      
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
      
      d2[d2 < 10^-4] = 0
      d2$rho_plus = d2$rho_1 + d2$rho_2
      d2=d2[,c(1:2,10:15)]
      colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","Scena","Delta","Branches","Rho_plus")
      write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Fixed_traits/Multistability_mapping_Scena_",
                            name_scena[scena_ID],"_branch_",branch,"_disp_",disp,".csv"),sep=";")
      d2=tibble()
    }
  }
}













### b) Analysis ----

d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.3,length.out=N_sim)
name_scena=c("global_C_local_F","global_C_global_F","local_C_local_F","local_C_global_F")
delta_seq=c(.1,.9)
branches=c("Degradation","Restoration")

for (branch in branches){
  
  for (disp in delta_seq) {
    
    for (scena_ID in 1:4){ #for each scenario of species pairs
      
      d2=rbind(d2,read.table(paste0("../Table/2_species/PA/Multistability_PA/Fixed_traits/Multistability_mapping_Scena_",name_scena[scena_ID],
                                    "_branch_",branch,"_disp_",disp,".csv"),sep=";"))
      
    }
  }
}
write.table(d2,paste0("../Table/2_species/PA/Multistability_fixed_traits_PA.csv"),sep=";")






d2=read.table(paste0("../Table/2_species/PA/Multistability_fixed_traits_PA.csv"),sep=";")
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

d2=d2[order(d2$alpha_0,d2$Scena,d2$Delta,d2$Stress),]

all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
  if (d2$state[x] != d2$state[x+1]){
    return(paste0(d2$state[x],"/", d2$state[x+1]))
  }
  else {return(d2$state[x])}
})

d_state=d2%>%
  filter(., Branches=="Degradation")%>%
  select(.,-Branches)
d_state$all_state=all_state

for (scale in c("local","global")){
  
  d2t=transform(d_state%>% 
                  filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
                  mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6)),
                Delta = factor(Delta, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
                Scena=factor(Scena,levels=c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")),
                             labels=c("Global facilitation","Local facilitation")))
  
  color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
  
  p=ggplot(d2t%>%
             mutate(all_state=recode_factor(all_state,
                                            "Desert/Coexistence"="Coexistence/Desert",
                                            "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                            "Desert/Stress_tolerant"="Stress_tolerant/Desert"))) +
    geom_tile(aes(x=Stress,y=alpha_0,fill=all_state))+
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'), fill = "") +
    theme(legend.text = element_text(size = 11))+
    scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                               "Competitive" = "#ACD87B", 
                               "Stress_tolerant" = "#7BD8D3",
                               "Stress_tolerant/Desert" ="#0F8E87",
                               "Coexistence/Stress_tolerant"="#9BBBB9",
                               "Coexistence/Desert"="#C19E5E",
                               "Desert"=  "#696969"
                               ))+
    facet_grid(Scena~Delta,labeller=labeller(Delta=label_parsed))+
    the_theme+theme(strip.text.x = element_text(size=12))
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Multistability_",scale,"_competition.pdf"),p,width = 7,height = 6)
  
  
  
  
  d2t=transform(d2%>%filter(., alpha_0 %in% c(0,.3,unique(d2$alpha_0)[50]))%>% 
                  filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F"))),
                Delta = factor(Delta, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
                Scena=factor(Scena,levels=c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")),
                             labels=c("Global facilitation","Local facilitation")),
                alpha_0 = factor(alpha_0, levels=c(0,unique(d2$alpha_0)[50],.3), labels=c("alpha[0] : 0", "alpha[0] : 0.15",
                                                                       "alpha[0] : 0.3")))
  
  
  
  
  p=ggplot(d2t%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
    geom_line(aes(x=Stress,y=value,linetype=Branches,color=variable))+
    facet_grid(Scena+Delta~alpha_0,labeller=labeller(Delta=label_parsed,alpha_0=label_parsed))+
    the_theme+theme(strip.text.x = element_text(size=12))+
    scale_color_manual(values=as.character(color_rho[c(4,2)]))+
    labs(x="Stress (S)",y="Densities",linetype="",color="")
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Bifu_bistability_species_",scale,"_competition.pdf"),width = 8,height = 7)
  
  p=ggplot(d2t)+
    geom_line(aes(x=Stress,y=Rho_plus,linetype=Branches))+
    facet_grid(Scena+Delta~alpha_0,labeller=labeller(Delta=label_parsed,alpha_0=label_parsed))+
    the_theme+theme(strip.text.x = element_text(size=12))+
    labs(x="Stress (S)",y="Densities",linetype="",color="")
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Bifu_bistability_community_",scale,"_competition.pdf"),width = 8,height = 7)
  
  set.seed(123)
  v=runif(2)
  d2t$CSI=v[1]*d2t$Stress_tolerant+v[2]*d2t$Competitive
  
  
  p=ggplot(d2t)+
    geom_point(aes(x=Stress,y=CSI),size=.5,alpha=.5,shape=21)+
    facet_grid(Scena+Delta~alpha_0,labeller=labeller(Delta=label_parsed,alpha_0=label_parsed))+
    the_theme+theme(strip.text.x = element_text(size=12))+
    labs(x="Stress (S)",y="Community index",linetype="",color="")
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/CSI_",scale,"_competition.pdf"),width = 8,height = 7)
  
}


## 4) Scale and strength of facilitation, dispersal and competition on niche expansion ----
### a) Do the simulation ----



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
write.table(d_niche,"../Table/2_species/PA/Niche_expansion_PA.csv",sep=";")

d_RNE[,6:7][is.na(d_RNE[,6:7])] = NA

write.table(d_RNE,"../Table/2_species/PA/RNE_PA.csv",sep=";")
write.table(d_all_dyn,"../Table/2_species/PA/All_dyn_PA_RNE.csv",sep=";")









### b) Analyse + plot graphics ----

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

N_rep = 50
S_seq = c(0,.1,.25,.4)
alpha_seq = seq(0, .3, length.out = N_rep)
f_seq=.9
delta_seq=c(.1, .9)
cintra_seq=c(.3,.5,.7)


name_scena=c("local_C_local_F","global_C_global_F","local_C_global_F","global_C_local_F")

d_clustering=tibble() #initializing the tibble

for (scena_ID in 1:4){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (aii in cintra_seq){ #varying intraspecific competition strength
      
      for (f in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) { #varying competition
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          
          param["cintra"]=aii
          param["f"]=f
          param["delta"]=disp
          param["alpha_0"]=alpha0
          
          
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
            Rho_1=d2$rho_1,Rho_2=d2$rho_2,Rho_12=d2$rho_12,
            Rho_22=d2$rho_22,Rho_11=d2$rho_11,
            S   = d2$S,alpha_0 = d2$alpha_0,
            f=f,delta=disp,Scena=scena_ID,
            cintra=aii
          ))

          
        } #end competition loop
        
      } #end facilitation loop
      
    } #end h loop
    
  } #end dispersal loop
  
} #end scenario loop

write.table(d_clustering,"../Table/2_species/PA/Clustering_PA.csv",sep=";")



### b) Analysis ----
d_clustering = read.table("../Table/2_species/PA/Clustering_PA.csv",sep=";")
fig_col=colorRampPalette(c("yellow","orange","red"))


d_clustering$c12=d_clustering$c11=d_clustering$c22=0
name_scena=c("local_C_local_F","global_C_global_F","local_C_global_F","global_C_local_F")

for (nr in 1:nrow(d_clustering)){
  
  rho_2=d_clustering$Rho_2[nr]
  rho_1=d_clustering$Rho_1[nr]
  rho_12=d_clustering$Rho_12[nr]
  rho_11=d_clustering$Rho_11[nr]
  rho_22=d_clustering$Rho_22[nr]
  
  if (rho_2 < 10^-2) rho_2 = 0#if (rho_2 < 10^-4) rho_2 = 0
  if (rho_1 < 10^-2) rho_1 = 0#if (rho_1 < 10^-4) rho_1 = 0
  
  if (rho_1 > 0 && rho_2 > 0){
    q2_1 = rho_12 / rho_1 #average number of species 2 around species 1 
    c21 = q2_1 / rho_2
    q1_2 = rho_12 / rho_2 #average number of species 1 around species 2
    c12 = q1_2 / rho_1
    
  } else{
    q1_2 = 0
    c12 = NaN
    q2_1 = 0
    c21 = NaN
  }
  
  if (rho_2 > 0){
    q2_2 = rho_22 / rho_2 #average number of species 2 around species 2 
    c22 = q2_2 / rho_2
  } else{
    q2_2 = 0
    c22 = NaN
  }
  
  if (rho_1 > 0){
    q1_1 = rho_1 / rho_1 #average number of species 1 around species 1 
    c11 = q1_1 / rho_1
  } else{
    q1_1 = 0
    c11 = NaN
  }
  
  d_clustering$c12[nr]=c12
  d_clustering$c22[nr]=c22
  d_clustering$c11[nr]=c11
  d_clustering$Scena[nr]=name_scena[as.numeric(d_clustering$Scena[nr])]
  
  
}
write.table(d_clustering,"../Table/2_species/PA/Clustering_PA.csv",sep=";")


pdf("../Figures/2_species/PA/Clustering_all_figs.pdf",width = 7,height = 6)

for (scena_name in unique(d_clustering$Scena)){
  
  for (disp in c(0.1,.9)){
    
    print(ggplot(d_clustering%>%
                   filter(., delta==disp,Scena==scena_name)%>%
                   melt(., measure.vars=c("Rho_1","Rho_2")))+
            geom_point(aes(x=as.numeric(alpha_0),y=value,color=variable),size=.5,alpha=.7)+
            geom_hline(yintercept = 0)+
            facet_grid(cintra~S,labeller=label_bquote(cols = Stress == .(S),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y="Densities")+
            scale_color_manual(values=as.character(color_rho[c(2,4)]))+
            ggtitle(paste0(scena_name,", delta = ",disp)))
    
    print(ggplot(d_clustering%>%
                   filter(., delta==disp,Scena==scena_name)%>%
                   group_by(., cintra,alpha_0,S,delta)%>%
                   summarise(.,.groups ="keep",c12=mean(c12) ))+
            geom_point(aes(x=as.numeric(alpha_0),y=c12))+
            geom_hline(yintercept = 1,linetype=9,lwd=.5)+
            facet_grid(cintra~S,labeller=label_bquote(cols = Stress == .(S),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{12}$"))+
            scale_y_log10())
    
    print(ggplot(d_clustering%>%
                   filter(., delta==disp,Scena==scena_name)%>%
                   group_by(., cintra,alpha_0,S,delta)%>%
                   summarise(.,.groups ="keep",c11=mean(c11) ))+
            geom_point(aes(x=as.numeric(alpha_0),y=c11))+
            geom_hline(yintercept = 1,linetype=9,lwd=.5)+
            facet_grid(cintra~S,labeller=label_bquote(cols = Stress == .(S),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{11}$"))+
            scale_y_log10())
    
    print(ggplot(d_clustering%>%
                   filter(., delta==disp,Scena==scena_name)%>%
                   group_by(., cintra,alpha_0,S,delta)%>%
                   summarise(.,.groups ="keep",c22=mean(c22) ))+
            geom_point(aes(x=as.numeric(alpha_0),y=c22))+
            geom_hline(yintercept = 1,linetype=9,lwd=.5)+
            facet_grid(cintra~S,labeller=label_bquote(cols = Stress == .(S),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{22}$"))+
            scale_y_log10())
    
  }
}

dev.off()


#Making a clean figure of the mechanisms


d_clustering = read.table("../Table/2_species/PA/Clustering_PA.csv",sep=";")

name_mesu=c("c12","c11","c22")
yname=c(TeX("$c_{12}$"),TeX("$c_{11}$"),TeX("$c_{22}$"))

for (measu in 1:3){
  assign(paste0("p_",measu),ggplot(d_clustering%>%
                                     filter(., Scena %in% name_scena[c(1,4)],S==.25,cintra==.3,alpha_0 %in% unique(d_clustering$alpha_0)[seq(1,50,by=2)])%>%
                                     melt(., measure.vars=name_mesu[measu])%>%
                                     mutate(., Scena=recode_factor(Scena,"global_C_local_F"="Global","local_C_local_F"="Local")))+
           geom_point(aes(x=as.numeric(alpha_0),y=value,shape=as.factor(delta),color=Scena))+
           the_theme+
           labs(x=TeX("$\\alpha_{e}$"),y=yname[measu],color="Scale competition",shape=TeX("$\\delta$"))+
           scale_shape_manual(values=c(0,1))+
           scale_color_manual(values=c("#8108A9","#DAAF42"))
  )
}
p_tot=ggarrange(p_2+ylim(2,5),p_3+geom_hline(yintercept = 1),p_1+geom_hline(yintercept = 1),ncol=3,common.legend = T,legend = "bottom")
ggsave("../Figures/2_species/PA/Clustering_mecanisms.pdf",width = 7,height = 4)




d_clustering = read.table("../Table/2_species/PA/Clustering_PA.csv",sep=";")


p=ggplot(d_clustering%>%
         filter(., Scena %in% name_scena[c(1,4)],S==.25,alpha_0 %in% unique(d_clustering$alpha_0)[seq(1,50,by=2)])%>%
         melt(., measure.vars="c12")%>%
         mutate(., Scena=recode_factor(Scena,"global_C_local_F"="Global","local_C_local_F"="Local")))+
  geom_point(aes(x=as.numeric(alpha_0),y=value,shape=as.factor(delta),color=Scena))+
  the_theme+
  facet_wrap(.~cintra,labeller=label_bquote(cols = alpha[0] == .(cintra) ))+
  labs(x=TeX("$\\alpha_{e}$"),y=TeX("$\\c_{12}$"),color="Scale competition",shape=TeX("$\\delta$"))+
  scale_shape_manual(values=c(0,1))+
  scale_color_manual(values=c("#8108A9","#DAAF42"))+
  geom_hline(yintercept = 1)+
  theme(strip.text.x = element_text(size=12))


ggsave("../Figures/2_species/PA/Clustering_mecanisms_intraspecific.pdf",p,width = 7,height = 4)



## 6) Net effects from partial derivative ----
### a) Simulation ----

tspan = c(0, 30000)
t = seq(0, 30000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
julia_assign("state", state)

N_sim = 100;epsilon=10^(-8)
name_scena=c("global_C_local_F","local_C_local_F")
disp_seq=c(.1,.9)
c_seq = seq(0, .3, length.out = N_sim)
S_seq = seq(0, 1, length.out = N_sim)

d_niche=d_RNE=d_all_dyn=tibble() 
d2 = d3 = tibble()
param=Get_PA_parameters()  
param=c(param[1:3],"beta1"=1,"beta2"=1,param[5:12])
param["cintra"]=.3
C_for_analyse = seq(0, .3, length.out = N_sim)[c(1,N_sim/2,N_sim)]



for (scena_ID in 1:2){ 
  
  for (comp in C_for_analyse) {
    
    for (disp in disp_seq){
      
      for (S in S_seq) {
        
        for (sp in 1:2) { #for both species
          
          param["S"] = S
          param["alpha_0"] = comp # update parameters
          param["delta"] = disp
          julia_assign("p", param)
          
          
          # first without press perturbation
          if (scena_ID==1){
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
          }else if (scena_ID==2){  
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F_press, state, tspan, p)")
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-5000):nrow(d),c(1,2)[-sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],
                                  Scena=name_scena[scena_ID],Disp=disp)) 
          d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, alpha_0 = comp))
          
          
          
          # with press perturbation on growth rate proxy
          
          param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
          julia_assign("p", param)
          
          
          if (scena_ID==1){
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
          }else if (scena_ID==2){  
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F_press, state, tspan, p)")
          }
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-5000):nrow(d),c(1,2)[-sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],
                                  Scena=name_scena[scena_ID],Disp=disp))
          
          param[paste0("beta",sp)] = 1
        }
      }
    }
  }
}
  
d2[d2 < epsilon] = 0
colnames(d2) = c("Eq", "S", "alpha_0", "Type","Species","Scena","Disp")


write.table(d2,"../Table/2_species/PA/Net_effects_partial_deriv_PA.csv",sep=";")



### b) Analysis ----

d2=read.table("../Table/2_species/PA/Net_effects_partial_deriv_PA.csv",sep=";")



N_sim = 100;epsilon=10^(-8)
c_seq = seq(0, .3, length.out = N_sim)
C_for_analyse = seq(0, .3, length.out = N_sim)[c(1,N_sim/2,N_sim)]

net_effect =sapply(seq(1, nrow(d2) , by = 2),function(x){
  return((d2$Eq[x+1] - d2$Eq[x]) / epsilon)
})


d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect

d_net=transform(d_net,
                Disp = factor(Disp, levels=c(.1,.9), labels=c("delta : 0.1", "delta : .9")),
                Scena=factor(Scena,levels=c("global_C_local_F","local_C_local_F"),
                             labels=c("Global competition","Local competition")))

for (id_plot in 1:3){
  assign(paste0("p_",id_plot),
         ggplot(d_net %>%filter(round(alpha_0,4)==round(C_for_analyse[id_plot]))) +
           geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
           the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
           geom_hline(yintercept = 0,linetype=9)+
           facet_grid(Scena~Disp,labeller=labeller(Disp=label_parsed))+
           labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial \beta_j}$)'))  )
}

p=ggarrange(p_1,p_2,p_3,ncol=3,common.legend = T,legend = "bottom")


ggsave(paste0("../Figures/2_species/MF/Net_effects_","main","_",func_MF,".pdf"),p_tot,width = 8,height = 5)


## 7) Comparing Net-effects and RII ----
### a) Simulation ----



tspan = c(0, 20000)
t = seq(0, 20000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
Nsim=100
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



## 8) Coexistence and trait difference ----
### a) Simulation ----
tspan = c(0, 6000)
t = seq(0, 6000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=100
S_seq = seq(0,1, length.out = N_sim)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)
dispersal_scale=c(.1,.9)
branches=c("Degradation","Restoration")

for (disp in dispersal_scale){
  

  for (psi1 in psi1_seq){
  
    for (branch in branches){
      d2 = tibble()
      
      if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
        state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
        S_seq=seq(0,1, length.out = N_sim)
      }else {
        state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
        S_seq=rev(seq(0,1, length.out = N_sim))
      }
      
      for (cinter in c_inter_seq){
      
        for (psi2 in psi_seq){
          
          for (S in S_seq) { 
            
            julia_assign("state", state)
            param=Get_PA_parameters()
            param["delta"]=disp
            param["cintra"]=.3
            param["alpha_0"]=cinter
            param["S"] = S
            param["psi_1"]=psi1
            param["psi_2"]=psi2
            julia_assign("p", param)
          
            prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
            
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d), ] %>% add_column(Stress = S, Psi2 = psi2, Psi1 = psi1,alpha_0=cinter,
                                                       Branch=branch))
            
          }
        } # end trait value 2nd species
        
        d2[d2 < 10^-4] = 0
        d2$rho_plus = d2$rho_1 + d2$rho_2
        d2=d2[,c(1,2,10:15)]
        colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "Psi2","Psi1","alpha_0","Branches","Rho_plus")
        write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                              cinter,"_branch_",branch,
                              "_Psi1_",psi1,"_delta_",disp,".csv"),sep=";")
        
      } # end loop interspecific competition
      
    } # end loop branch
    
  } # end loop first species trait
  
} #end loop dispersal



### b) Analysis ----

c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)
dispersal_scale=c(.1,.9)

for (disp in dispersal_scale){

  for (Psi_sp1 in psi1_seq){
    d=tibble()  
    
    for (branch in c("Restoration","Degradation")){
      for (cinter in c_inter_seq){
        
        d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                            cinter,"_branch_",branch,
                            "_Psi1_",Psi_sp1,"_delta_",disp,".csv"),sep=";")
        
        d=rbind(d,d2)
        
      } # end loop interspecific competition
    } # end loop branch
    
    
    
    
    d2=d
    d2[,1:2][d2[,1:2] < 10^-4] = 0
    
    
    #COmputing CSI index
    set.seed(123)
    u=runif(2)
    d2$CSI = sapply(1:nrow(d2),function(x){
      return(u[1]*d2$Stress_tolerant[x]+u[2]*d2$Competitive[x])
    })
    
    
    
    if (Psi_sp1 ==1) {
      
      
      
      d2$state = sapply(1:nrow(d2), function(x) {
        if (d2[x, 1] > 0 & d2[x, 2] > 0) {
          return("Coexistence")
        }
        if (d2[x, 1] > 0 & d2[x, 2] == 0) {
          return("Stress_tolerant")
        }
        if (d2[x, 1] == 0 & d2[x, 2] > 0) {
          return("Species 2")
        }
        if (d2[x, 1] == 0 & d2[x, 2] == 0) {
          return("Desert")
        }
      })
      
      d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1),]
      
      all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
        if (d2$state[x] != d2$state[x+1]){
          return(paste0(d2$state[x],"/", d2$state[x+1]))
        }
        else {return(d2$state[x])}
      })
      
      d_state=d2%>%
        filter(., Branches=="Degradation")%>%
        select(.,-Branches)
      d_state$all_state=all_state
      
      
      color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
      
      appender <- function(string) {
        TeX(paste("$\\alpha_e = $", string))}
      
      
      
  
      p=ggplot(d_state%>%
                 mutate(all_state=recode_factor(all_state,
                                                "Species 2/Coexistence"="Coexistence/Species 2",
                                                "Desert/Species 2"="Species 2/Desert",
                                                "Desert/Coexistence"="Coexistence/Desert",
                                                "Desert/Stress_tolerant"="Stress_tolerant/Desert",
                                                "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                                "Stress_tolerant/Species 2"="Species 2/Stress_tolerant",
                                                "Species 2/Coexistence"="Coexistence/Species 2"))) +
        geom_tile(aes(x=Stress,y=as.numeric(1-Psi2),fill=all_state))+
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(x = "Stress (S)", y = TeX(r'(Trait difference \ |$\psi_1-\psi_2|)'), fill = "") +
        theme(legend.text = element_text(size = 11))+
        scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                   "Species 2" = "#ACD87B",
                                   "Coexistence/Species 2" = "#DDEFCA",
                                   "Stress_tolerant" = "#7BD8D3",
                                   "Stress_tolerant/Desert" ="#0F8E87",
                                   "Coexistence/Stress_tolerant"="#9BBBB9",
                                   "Coexistence/Desert"="#C19E5E",
                                   "Desert"=  "#696969",
                                   "Species 2/Stress_tolerant" = "#9B68A0"),
                          labels=c("Coexistence","Species 2","Coexistence/Species 2","Stress-tolerant","Stress-tolerant/Desert",
                                   "Coexistence/Stress-tolerant",
                                   "Coexistence/Desert","Desert","Species 2/Stress-tolerant"))+
        facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
        the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
      
      
      ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Multistability_trait_variation_sp_Psi1_",Psi_sp1,"_dispersal_",disp,".pdf"),p,width = 9,height = 4)
      
      
      
      d2t=transform(d2,
                    alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                            "alpha[e] == 0.5")),
                    Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                                labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
      
      
      
      
      #Bifurcation diagram
      p=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                 filter(., Psi2 %in% unique(d_state$Psi2)[c(1,75,100)])%>%
                 mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                 filter(., alpha_0!=.1))+
        geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
        labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
        the_theme +
        scale_color_manual(values = color_rho[c(2, 4)],labels=c("Species 2","Stress-tolerant")) +
        scale_linetype_manual(values=c(1,9))+
        theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
        facet_grid(Psi2~alpha_0, scales = "free",  
                   labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
      
      
      
      ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Bifurcation_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,".pdf"),p,width = 7,height = 4)
      
      
      
      
      #Community index
      p=ggplot(d2%>%
                 filter(., Psi2 %in% unique(d2$Psi2)[c(1,75,100)])%>%
                 mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                 filter(., alpha_0!=.1))+
        geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
        labs(x = "Stress (S)", y = "Community index", color = "",linetype="") +
        the_theme +
        theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
        facet_grid(Psi2~alpha_0, scales = "free",  
                   labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
      
      ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/CSI_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,".pdf"),p,width = 7,height = 4)
      
      
      
      
      
      
    } 
    if (Psi_sp1==0){
      
      d2$state = sapply(1:nrow(d2), function(x) {
        if (d2[x, 1] > 0 & d2[x, 2] > 0) {
          return("Coexistence")
        }
        if (d2[x, 1] > 0 & d2[x, 2] == 0) {
          return("Competitive")
        }
        if (d2[x, 1] == 0 & d2[x, 2] > 0) {
          return("Species 2")
        }
        if (d2[x, 1] == 0 & d2[x, 2] == 0) {
          return("Desert")
        }
      })
      
      d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1),]
      
      all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
        if (d2$state[x] != d2$state[x+1]){
          return(paste0(d2$state[x],"/", d2$state[x+1]))
        }
        else {return(d2$state[x])}
      })
      
      d_state=d2%>%
        filter(., Branches=="Degradation")%>%
        select(.,-Branches)
      d_state$all_state=all_state
      
      
      color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
      
      appender <- function(string) {
        TeX(paste("$\\alpha_e = $", string))}
      
  
      
      p=ggplot(d_state%>%
                 mutate(all_state=recode_factor(all_state,
                                                "Species 2/Coexistence"="Coexistence/Species 2",
                                                "Desert/Species 2"="Species 2/Desert",
                                                "Desert/Coexistence"="Coexistence/Desert",
                                                "Desert/Competitive"="Competitive/Desert",
                                                "Competitive/Coexistence" = "Coexistence/Competitive"
                                                ),
                        Psi2=round(Psi2,5),
                        Stress=round(Stress,5))) +
        geom_tile(aes(x=Stress,y=abs(as.numeric(Psi2)),fill=all_state))+
        theme_classic() +
        theme(legend.position = "bottom") +
        labs(x = "Stress (S)", y = TeX(r'(Trait difference \ |$\psi_1-\psi_2|)'), fill = "") +
        theme(legend.text = element_text(size = 11))+
        scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                   "Species 2" = "#7BD8D3",
                                   "Coexistence/Species 2" = "#9BBBB9",
                                   "Species 2/Desert" ="#0F8E87",
                                   "Coexistence/Desert"="#C19E5E",
                                   "Desert"=  "#696969"))+
        facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
        the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
      
      
      ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Multistability_trait_variation_sp_Psi1_",Psi_sp1,"_dispersal_",disp,".pdf"),p,width = 9,height = 4)
      
      
      
      #Bifurcation diagrams
      
      d2t=transform(d2,
                    alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                            "alpha[e] == 0.5")),
                    Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                                labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
      
      
      p=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                 filter(., Psi2 %in% unique(d_state$Psi2)[c(1,75,100)])%>%
                 mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                 filter(., alpha_0!=.1))+
        geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
        labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
        the_theme +
        scale_color_manual(values = color_rho[c(2, 4)],labels=c("Species 2","Stress-tolerant")) +
        scale_linetype_manual(values=c(1,9))+
        theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
        facet_grid(Psi2~alpha_0, scales = "free",  
                   labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
      
      
      
      ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Bifurcation_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,".pdf"),p,width = 7,height = 4)
      
      
      #COmmunity index
      p=ggplot(d2%>%
                 filter(., Psi2 %in% unique(d2$Psi2)[c(1,75,100)])%>%
                 mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                 filter(., alpha_0!=.1))+
        geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
        labs(x = "Stress (S)", y = "Community index", color = "",linetype="") +
        the_theme +
        theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
        facet_grid(Psi2~alpha_0, scales = "free",  
                   labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
      
      ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/CSI_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,".pdf"),p,width = 7,height = 4)
      
    }
  }
}



## 9) Functional diversity ----

N_rep = 50
h_seq=c(1,1.5)[1]
f_seq=c(0,.3,.9)

id_for_RNE = seq(1:N_rep)[c(1, length(seq(1:N_rep))/2, length(seq(1:N_rep)))]

name_scena=c("local_C_local_F","global_C_global_F",
             "local_C_global_F","global_C_local_F")
pdf("../Figures/2_species/PA/Functional_diversity.pdf",width = 7,height = 5)

for (scena_ID in c(1:4)[c(2,4)]){ #for each scenario of species pairs
  
  for (f in f_seq){ #varying facilitation strength
    
    for (h in h_seq){ #varying competitive advantage strength
      
      for (disp in c(.1, .9)){
        
        d_diversity=tibble()
        name_fig=paste0("_d_",disp,"_f_",f,"_h_",h,"_scena_",name_scena[scena_ID])
        
        #loading data
        d_state=read.table(paste0("../Table/2_species/PA/Sim_PA_scales/2_species_PA",name_fig,".csv"),sep=";")
        
        for (x in 1:nrow(d_state)){
          
          densi=d_state[x,1:2]
          if (any(densi==0)){
            d_diversity=rbind(d_diversity,tibble(FD_0 =0, FD_1=0,  FD_2 =0,  D_0=0,   D_1=0,   D_2=0,)%>%
                                add_column(., 
                                           Stress=d_state$S[x],
                                           alpha_0=d_state$alpha_0[x]))
          } else{
            d_diversity=rbind(d_diversity,Get_diversity_community(trait =c(1,0) ,densities =densi )%>%
                                add_column(., 
                                           Stress=d_state$S[x],
                                           alpha_0=d_state$alpha_0[x]))
          }
        }
        
        print(ggplot(d_diversity %>% 
                 melt(.,id.vars=c("Stress","alpha_0")))+
          geom_tile(aes(x=Stress,y=alpha_0,fill=value))+
          facet_wrap(.~variable,scales = "free")+
          scale_fill_viridis_c()+
          labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'))+
          the_theme+ggtitle(paste0("Figure for ",name_fig)))
        
        
      }
    }
  }
}
dev.off()

# Step 3) CA between both species : example and clustering ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
#

## 1) Clustering of species  ----



#Simulations were made on Julia 2_species_clustering.jl file
d_clustering=tibble()
count_u=1

for (disp in c(0.1,0.9)){
  
  
  list_all_files=list.files(paste0("../Table/2_species/CA/Clustering_species/",disp))
  
  for (i in list_all_files){
    landscape=read.table(paste0("../Table/2_species/CA/Clustering_species/",disp,"/",i),sep=",")
    landscape=as.matrix(landscape)
    neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    
    
    rho_1 = length(which(landscape == 1)) / length(landscape)
    rho_2 = length(which((landscape == 2))) / length(landscape)
    
    if (rho_2 < 10^-4) rho_2 = 0
    if (rho_1 < 10^-4) rho_1 = 0
    
    if (rho_1 > 0 && rho_2 > 0){
       q2_1 = mean(neigh_2[which((landscape == 1))] / 4) #average number of species 2 around species 1 
       c21 = q2_1 / rho_2
       q1_2 = mean(neigh_1[which((landscape == 2))] / 4) #average number of species 1 around species 2
       c12 = q1_2 / rho_1
       
    } else{
       q1_2 = 0
       c12 = NaN
       q2_1 = 0
       c21 = NaN
    }
    
    if (rho_2 > 0){
      q2_2 = mean(neigh_2[which((landscape == 2))] / 4) #average number of species 2 around species 2 
      c22 = q2_2 / rho_2
    } else{
      q2_2 = 0
      c22 = NaN
    }
    
    if (rho_1 > 0){
      q1_1 = mean(neigh_1[which((landscape == 1))] / 4) #average number of species 1 around species 1 
      c11 = q1_1 / rho_1
    } else{
      q1_1 = 0
      c11 = NaN
    }
    
    # clustering total vegetation
    
    binary_landscape = landscape
    binary_landscape[binary_landscape >  0] = 1
    binary_landscape[binary_landscape <= 0] = 0
    rho_p = length(which((binary_landscape == 1))) / length(binary_landscape)
  
    neigh_vege= simecol::neighbors(x =binary_landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    
  
    cpp = (mean(neigh_vege[which((binary_landscape == 1))] / 4)) / rho_p
    
    
    
    
    d_clustering=rbind(d_clustering,tibble(c12=c12,c21=c21,c11=c11,c22=c22,cpp=cpp,Stress=as.numeric(strsplit(i,"_")[[1]][5]),
                                           Rho_1=rho_1,Rho_2=rho_2,Rho_plus=rho_p,
                                           alpha_0=as.numeric(strsplit(i,"_")[[1]][7]),
                                           Scale_comp=strsplit(i,"_")[[1]][3],
                                           cintra=as.numeric(strsplit(i,"_")[[1]][9]),
                                           nsave=as.numeric(gsub(x=strsplit(i,"_")[[1]][11],replacement = "",pattern = ".csv")),
                                           Dispersal=disp))
    print(count_u)
    count_u=count_u+1
  }
}

write.table(d_clustering,"../Table/2_species/CA/Clustering_species.csv",sep=";")



d_clustering = read.table("../Table/2_species/CA/Clustering_species.csv",sep=";")
pdf("../Figures/2_species/CA/Clustering.pdf",width = 7,height = 4)
  
for (scale_comp in c("global","local")){
  for (disp in c(0.1,.9)){
    
    
    print(ggplot(d_clustering%>%
                   filter(., Dispersal==disp,Scale_comp==scale_comp)%>%
                   melt(., measure.vars=c("Rho_1","Rho_2")))+
            geom_point(aes(x=as.numeric(alpha_0),y=value,color=variable),size=.5,alpha=.7)+
            geom_hline(yintercept = 0)+
            facet_grid(cintra~Stress,labeller=label_bquote(cols = Stress == .(Stress),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{12}$"))+
            scale_color_manual(values=rev(as.character(color_rho[c(2,4)])))+
            ggtitle(paste0("Scale competition = ",scale_comp,", fraction local dispersal = ",100  * (1-disp))))
    
    
    
    print(ggplot(d_clustering%>%
                   filter(., Dispersal==disp,Scale_comp==scale_comp)%>%
                   group_by(., cintra,alpha_0,Stress,Dispersal)%>%
                   summarise(.,.groups ="keep",c12=mean(c12) ))+
            geom_point(aes(x=as.numeric(alpha_0),y=c12))+
            geom_hline(yintercept = 1,linetype=9,lwd=.5)+
            facet_grid(cintra~Stress,labeller=label_bquote(cols = Stress == .(Stress),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{12}$"))+
      scale_y_log10())
    
    print(ggplot(d_clustering%>%
                   filter(., Dispersal==disp,Scale_comp==scale_comp)%>%
                   group_by(., cintra,alpha_0,Stress,Dispersal)%>%
                   summarise(.,.groups ="keep",c21=mean(c21) ))+
            geom_point(aes(x=as.numeric(alpha_0),y=c21))+
            geom_hline(yintercept = 1,linetype=9,lwd=.5)+
            facet_grid(cintra~Stress,labeller=label_bquote(cols = Stress == .(Stress),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{21}$"))+ylim(0,15)+
      scale_y_log10())
    
    print(ggplot(d_clustering%>%
                   filter(., Dispersal==disp,Scale_comp==scale_comp)%>%
                   group_by(., cintra,alpha_0,Stress,Dispersal)%>%
                   summarise(.,.groups ="keep",cpp=mean(cpp) ))+
            geom_point(aes(x=as.numeric(alpha_0),y=cpp))+
            geom_hline(yintercept = 1,linetype=9,lwd=.5)+
            facet_grid(cintra~Stress,labeller=label_bquote(cols = Stress == .(Stress),rows = alpha[ii] == .(cintra) ))+
            the_theme+labs(x=TeX("$\\alpha_e$"),y=TeX("$\\c_{++}$"))+
      scale_y_log10())
  }
}

dev.off()



d_clustering = read.table("../Table/2_species/CA/Clustering_species.csv",sep=";")
d_clustering$c22[which(d_clustering$Rho_2<.01)]=NA
name_mesu=c("c12","c11","c22")
yname=c(TeX("$c_{12}$"),TeX("$c_{11}$"),TeX("$c_{22}$"))

for (measu in 1:3){
  assign(paste0("p_",measu),ggplot(d_clustering%>%
                                     filter(., Stress==.25,cintra==.3)%>%
                                     melt(., measure.vars=name_mesu[measu])%>%
                                     mutate(., Scale_comp=recode_factor(Scale_comp,"global"="Global","local"="Local"))%>%
                                     group_by(., Scale_comp,Dispersal,alpha_0)%>%
                                     summarise(., .groups = "keep",value=mean(value)))+
           geom_line(aes(x=as.numeric(alpha_0),y=value,color=Scale_comp,group=interaction(Scale_comp,Dispersal)),stat="smooth",alpha=.4,lwd=.8)+
           geom_point(aes(x=as.numeric(alpha_0),y=value,shape=as.factor(Dispersal),color=Scale_comp))+
           the_theme+
           labs(x=TeX("$\\alpha_{e}$"),y=yname[measu],color="Scale competition",shape=TeX("$\\delta$"))+
           scale_shape_manual(values=c(0,1))+
           scale_color_manual(values=c("#8108A9","#DAAF42"))
  )
}
p_tot=ggarrange(p_2+geom_hline(yintercept = 1),
                p_3+geom_hline(yintercept = 1)+ylim(.9,5),
                p_1+geom_hline(yintercept = 1),ncol=3,common.legend = T,legend = "bottom")
ggsave("../Figures/2_species/CA/Clustering_mecanisms_CA.pdf",width = 7,height = 4)





## 2) Patch size distribution ----
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
    
    landscape=as.matrix(read.table(paste0("../Table/2_species/CA/Patch_size/Example_sim/Landscape_size_100_stress_",stress,"_",rep,".csv"),sep=","))
    
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
ggsave("../Figures/2_species/CA/Patch_size/Example_PSD.pdf",p_tot,width = 7,height = 10)



### b) Analyzing the simulations for different competition regimes ----


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
        
        landscape=as.matrix(read.table(paste0("../Table/2_species/CA/Patch_size/Competition_regime/Landscape_size_100_stress_",stress,"_",scena,"_",rep,"_relat_compet_",alpha0,".csv"),sep=","))
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
    ggsave(paste0("../Figures/2_species/CA/Patch_size/Patch_size_distribution_",scena,"_Intersp_comp_",alpha0,".pdf"),p_tot,width = 7,height = 8)
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

    write.table(d_PL,paste0("../Table/2_species/CA/PL_summary/PL_",scena,"_S_",stress,".csv"),sep=";")

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
  
  ggsave(paste0("../Figures/2_species/CA/Patch_size/Max_patch_size_by_species_",scena,".pdf"),p,width = 7,height = 4)
  
  
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
  
  ggsave(paste0("../Figures/2_species/CA/Patch_size/Max_patch_size_by_PL_class_",scena,".pdf"),p2,width = 7,height = 5)
  
  
  
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
  
  ggsave(paste0("../Figures/2_species/CA/Patch_size/PL_exponent_by_species_",scena,".pdf"),p3,width = 7,height = 4)
  
  
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
  
  ggsave(paste0("../Figures/2_species/CA/Patch_size/PL_exponent_by_PL_class_",scena,".pdf"),p4,width = 7,height = 5)
  
  
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
  write.table(Spatial_EWS,paste0("../Table/2_species/CA/EWS_spatial/Spatial_EWS_",scena,".csv"),sep=";")
  
  
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


