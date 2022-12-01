# Other) Recruitment rate & 4 grid matrices----
rm(list = ls())
source("./Dryland_shift_functions.R")
d = expand_grid(emax = 1, e = .1, S = seq(0, 1, length.out = 200), psi = c(0, .5, 1))

p = ggplot(d) +
    geom_line(aes(x = S, y = emax * (1 - S * (1 - e * psi)), color = as.factor(psi), group = psi),size=1) +
    scale_color_manual(values = color_Nsp(9)[c(2,5,8)]) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
        x = "Stress, S",
        y = TeX(r'($Recruitment rate,  (1-S (1-e \psi_i)))'), color = TeX("$\\psi_i$")
    )



ggsave("../Figures/Recruitment_rate.pdf", p, width = 6, height = 4)



# Step 1) Two species MF model ----
rm(list = ls())
source("./Dryland_shift_functions.R")
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
colnames(d2) = c("Eq", "S", "alpha_0", "Type","Species","Branches")

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

ggarrange(
  ggplot(d3 %>%melt(.,measure.vars=c("V1","V2"))%>%
         mutate(., alpha_0=as.factor(round(alpha_0,2)))) +
  geom_line(aes(x = S, y = value, color = variable),lwd=1,alpha=.5) +
  the_theme+scale_color_manual(values=c("blue","green"),labels=c("V1"= "Stess-tolerant","V2"="Competitive"))+
  geom_hline(yintercept = 0,linetype=9)+
  facet_grid(Branches~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
  labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial \beta_j}$)'))
  ,
  ggplot(d_net %>%
           mutate(., alpha_0=as.factor(round(alpha_0,2)))) +
  geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
  the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess-tolerant","2"="Competitive"))+
  geom_hline(yintercept = 0,linetype=9)+
  facet_grid(Branches~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
  labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial \beta_j}$)'))
  ,nrow=2)

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
c_seq=seq(0,.4,length.out=N_sim)
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
    
    for (S in S_seq) { 
      
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





## 7) Neighboring effects ----
### a) Simulation ----

tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
Nsim=100
S_seq = seq(0, 1, length.out = Nsim)
c_seq = seq(0,.3, length.out = Nsim)
epsilon=10^(-9)

d_RNE=tibble()


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



d_all=d_RNE
save(file="../Table/2_species/MF/Comparing_net_effects.RData",d_all)

### b) Analysis ----
load(file="../Table/2_species/MF/Comparing_net_effects.RData")
d_net=d_all



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
branches=c("Degradation","Restoration")


  
  
for (psi1 in psi1_seq){
  
  for (branch in branches){
    
    if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
      state = Get_MF_initial_state(c(.4,.4,.1))
      S_seq=seq(0,1, length.out = N_sim)
    }else {
      state = Get_MF_initial_state(c(.005,.005,.49))
      S_seq=rev(seq(0,1, length.out = N_sim))
    }
    
    for (cinter in c_inter_seq){
      d2 = tibble()
      
      
      for (psi2 in psi_seq){
        
        for (S in S_seq) { 
          
          julia_assign("state", state)
          param=Get_MF_parameters()
          param["cintra"]=.3
          param["alpha_0"]=cinter
          param["S"] = S
          param["psi_1"]=psi1
          param["psi_2"]=psi2
          julia_assign("p", param)
          
          prob = julia_eval("ODEProblem(MF_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(Stress = S, Psi2 = psi2, Psi1 = psi1,alpha_0=cinter,
                                                     Branch=branch))
          
        }
      } # end trait value 2nd species
      
      d2[d2 < 10^-4] = 0
      d2$rho_plus = d2$rho_1 + d2$rho_2
      colnames(d2) = c("Stress_tolerant", "Competitive","Degraded","Fertile", "Stress", "Psi2","Psi1","alpha_0","Branches","Rho_plus")
      write.table(d2,paste0("../Table/2_species/MF/Varying_traits/Multistability_varying_trait_interspe_comp_",
                            cinter,"_branch_",branch,
                            "_Psi1_",psi1,".csv"),sep=";")
      
    } # end loop interspecific competition
    
  } # end loop branch
  
} # end loop first species trait
  



### b) Analysis ----

c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)


  
for (Psi_sp1 in psi1_seq){
  d=tibble()  
  
  for (branch in c("Restoration","Degradation")){
    for (cinter in c_inter_seq){
      
      d2=read.table(paste0("../Table/2_species/MF/Varying_traits/Multistability_varying_trait_interspe_comp_",
                           cinter,"_branch_",branch,
                           "_Psi1_",Psi_sp1,".csv"),sep=";")
      
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
                                 "Species 2/Stress_tolerant" = "#C998CE"),
                        labels=c("Coexistence","Species 2","Coexistence/Species 2","Stress-tolerant","Stress-tolerant/Desert",
                                 "Coexistence/Stress-tolerant",
                                 "Coexistence/Desert","Desert","Species 2/Stress-tolerant"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
    
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/Multistability_trait_variation_sp_Psi1_",Psi_sp1,".pdf"),p,width = 9,height = 4)
    
    
    
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
      scale_color_manual(values = as.character(color_rho[c(2, 4)]),labels=c("Species 2","Stress-tolerant")) +
      scale_linetype_manual(values=c(1,9))+
      theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
      facet_grid(Psi2~alpha_0, scales = "free",  
                 labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
    
    
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/Bifurcation_varying_traits_Psi1_",Psi_sp1,".pdf"),p,width = 7,height = 4)
    
    
    
    
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
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/CSI_varying_traits_Psi1_",Psi_sp1,".pdf"),p,width = 7,height = 4)
    
    
    
    #Restoration and degradation points
    d_tipping=tibble()
    for (a0 in unique(d2$alpha_0)){
      for (psi in unique(d2$Psi2)){
        for (direction in unique(d2$Branches)){
          
          d_a0=filter(d2,alpha_0==a0,Psi2==psi,Branches==direction)
          
          if (direction=="Restoration") d_a0=d_a0[order(d_a0$Stress,decreasing = T),]
          
          density_C=d_a0$Competitive
          stress_C=stress_ST=d_a0$Stress
          
          if (direction=="Restoration"){ #initially absent
            
            stress_C=stress_C[-c(1:min(which(abs(diff(density_C))>0)))]
            density_C=density_C[-c(1:min(which(abs(diff(density_C))>0)))]
            Tipping_C = stress_C[1]
            
          } else{
            
            Tipping_C = stress_C[min(which(density_C==0))-1]
            
          }
          
          d_tipping=rbind(d_tipping,tibble(Competition=a0,Psi2=psi,Branch=direction,
                                           Tipping_C = Tipping_C))
          
        }
      }
    }
    
    p=ggplot(d_tipping)+
      geom_smooth(aes(x=Psi2,y=Tipping_C,color=Branch,group=interaction(Competition,Branch)),se = F)+
      the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
      labs(x=TeX(r'(Trait species 2, \ $\psi_2)'),y="Tipping point",color="")+
      scale_color_manual(values=c("black","blue"))+
      scale_x_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)'),
                                             breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1")),
                         breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1"))
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/Tipping_points_varying_traits_Psi1_",Psi_sp1,".pdf"),p,width = 7,height = 4)
    
    
    
    
    
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
                                 "Competitive" = "#ACD87B",
                                 "Coexistence/Species 2" = "#9BBBB9",
                                 "Species 2/Desert" ="#0F8E87",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Coexistence/Competitive" = "#DDEFCA",
                                 "Coexistence/Species 2"="#9BBBB9",
                                 "Desert"=  "#696969"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
    
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/Multistability_trait_variation_sp_Psi1_",Psi_sp1,".pdf"),p,width = 9,height = 4)
    
    
    
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
      scale_color_manual(values = rev(color_rho[c(2, 4)]),labels=c("Species 2","Competitive")) +
      scale_linetype_manual(values=c(1,9))+
      theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
      facet_grid(Psi2~alpha_0, scales = "free",  
                 labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
    
    
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/Bifurcation_varying_traits_Psi1_",Psi_sp1,".pdf"),p,width = 7,height = 4)
    
    
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
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/CSI_varying_traits_Psi1_",Psi_sp1,".pdf"),p,width = 7,height = 4)
    
    
    
    #Restoration and degradation points
    d_tipping=tibble()
    for (a0 in unique(d2$alpha_0)){
      for (psi in unique(d2$Psi2)){
        for (direction in unique(d2$Branches)){
          
          d_a0=filter(d2,alpha_0==a0,Psi2==psi,Branches==direction)
          
          if (direction=="Restoration") d_a0=d_a0[order(d_a0$Stress,decreasing = T),]
          
          density_ST=d_a0$Stress_tolerant
          stress_ST=stress_ST=d_a0$Stress
          
          if (direction=="Restoration"){ #initially absent
            
            stress_ST=stress_ST[-c(1:min(which(abs(diff(density_ST))>0)))]
            density_ST=density_ST[-c(1:min(which(abs(diff(density_ST))>0)))]
            Tipping_ST = stress_ST[1]
            
          } else{
            
            Tipping_ST = stress_ST[min(which(density_ST==0))-1]
            
          }
          
          d_tipping=rbind(d_tipping,tibble(Competition=a0,Psi2=psi,Branch=direction,
                                           Tipping_ST = Tipping_ST))
          
        }
      }
    }
    
    p=ggplot(d_tipping)+
      geom_smooth(aes(x=Psi2,y=Tipping_ST,color=Branch,group=interaction(Competition,Branch)),se = F)+
      the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
      labs(x=TeX(r'(Trait species 1, \ $\psi_1)'),y="Tipping point",color="")+
      scale_color_manual(values=c("black","blue"))+
      scale_x_continuous(sec.axis = sec_axis(trans = ~ (.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)'),
                                             breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1")),
                         breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1"))
    
    ggsave(paste0("../Figures/2_species/MF/Multistability/Tipping_points_varying_traits_Psi1_",Psi_sp1,".pdf"),p,width = 7,height = 4)
    
  }
}




## 9) Niche species, varying traits ----
### a) Simulation ----

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 100
S_seq = seq(0,1,length.out=N_rep)
alpha_seq = seq(0, .3, length.out = 3)
f_seq=c(.3,.9)
trait_sp2_seq=seq(0,1,length.out=50)
trait1_seq=c(0,1)



d_niche=tibble() #initializing the tibble

for (psi_sp1 in trait1_seq){ #varying the trait of sp1
  
  for (facil in f_seq){ #varying facilitation strength
    
    for (alpha0 in alpha_seq) { #varying competition
      
      for (trait_2 in trait_sp2_seq){
        
        #Setting the parameters
        param=Get_MF_parameters()
        
        param["f"]=facil
        param["alpha_0"]=alpha0
        param["psi_1"]=psi_sp1 #we fixed to stress-tolerant species 
        param["psi_2"]=trait_2 #and vary the other species
        
        
        
        # 1) Competitive species alone
        
        state=Get_MF_initial_state(c(0,.8,.1))
        julia_assign("state", state)
        d2 = tibble()
        
        for (S in S_seq) {
          
          param["S"] = S
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(MF_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m","rho_0")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Competitive"))
          
          
        }
        d2[d2 < 10^-4] = 0
        
        
        S_critic2_alone=abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) #range of values where competitive species is
        
        
        
        # 2) Stress-tolerant species alone
        
        state=Get_MF_initial_state(c(0.8,0,.1))
        julia_assign("state", state)
        d2 = tibble()
        
        for (S in S_seq) {
          
          param["S"] = S
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(MF_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m","rho_0")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Stress-tolerant"))
          
          
        }
        d2[d2 < 10^-4] = 0
        
        
        S_critic1_alone=abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
        
        
        
        #3) Coexisting species
        
        state=Get_MF_initial_state(c(0.4,0.4,.1))
        julia_assign("state", state)
        d2 = tibble()
        
        for (S in S_seq) { #varying the stress 
          
          param["S"] = S
          julia_assign("p", param)
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(MF_two_species_varying_trait, state, tspan, p)")
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          
          colnames(d) = c("rho_1", "rho_2", "rho_m","rho_0")
          
          d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, Sp="Both"))
          
        }
        d2[d2 < 10^-4] = 0
        d2$rho_plus = d2$rho_1 + d2$rho_2
        
        S_critic1_both = abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
        S_critic2_both = abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) 
        
        
        d_niche=rbind(d_niche,tibble(
          Facilitation = facil, alpha_0 = alpha0, Psi1 = psi_sp1, Psi2 = trait_2,
          Delta_niche_1 = 100 * (S_critic1_both - S_critic1_alone) / S_critic1_alone,  #making it a percentage of initial niche
          Delta_niche_2 = 100 * (S_critic2_both - S_critic2_alone) / S_critic2_alone)) #making it a percentage of initial niche
        
        
      } #end varying trait loop
      
    } #end competition loop
    
  } #end facilitation loop
  
} #end sp1 trait loop

write.table(d_niche,"../Table/2_species/MF/Niche_expansion_varying_traits.csv",sep=";")




### b) Analysis ----

d_niche = read.table("../Table/2_species/MF/Niche_expansion_varying_traits.csv",sep=";")
d_niche = do.call(data.frame,lapply(d_niche,function(x) replace(x, is.infinite(x), NA)))


p=ggplot(NULL)+
  geom_path(data=d_niche%>%
              melt(., measure.vars=c("Delta_niche_1","Delta_niche_2"))%>%
              mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
              mutate(., Delta_psi=abs(Psi1-Psi2),
                     Facilitation=as.factor(Facilitation),
                     alpha_0=as.factor(alpha_0)),
            aes(x=Delta_psi,y=value,group=interaction(alpha_0,Facilitation),color=alpha_0),lwd=1)+
  geom_point(data=d_niche%>%
               melt(., measure.vars=c("Delta_niche_1","Delta_niche_2"))%>%
               mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
               mutate(., Delta_psi=abs(Psi1-Psi2),
                      Facilitation=as.factor(Facilitation),
                      alpha_0=as.factor(alpha_0))%>%
               filter(., Psi2 %in% unique(.$Psi2)[seq(1,length(unique(.$Psi2)),by =5)]), #to only have 10 points
             aes(x=Delta_psi,y=value,shape=Facilitation,color=alpha_0),fill="white",size=2)+
  facet_grid(variable~Psi1,scales="free",labeller = label_bquote(cols = psi[1] == .(Psi1)))+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=c("#91D4DE","#26A8D2","#06247B"))+
  labs(y = "Niche change (%)", x = TeX(r'(Trait difference \ |$\psi_1-\psi_2|)'),
       color=TeX(r'(Competition strength \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f)')) +
  the_theme

ggsave("../Figures/2_species/MF/Niche_expansion_varying_traits.pdf",p,width = 7,height = 6)




# Step 2) Pair approximation (PA) ----
rm(list = ls())
source("./Dryland_shift_functions.R")
julia_setup()
de = diffeq_setup()
#

## 1) Exploration along competitive ability ----
### a) Do the simulation ----
dir.create("../Table/2_species/PA/Sim_PA_scales",showWarnings = F)
# As we can get oscillations, function to get the mean densities

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
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

tspan = c(0, 2000)
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


n_random_ini=20
d2 = tibble()
S_seq = seq(0,1, length.out = 30)[1]
c_seq=seq(0,.3,length.out=30)[30]
name_scena=c("global_C_local_F","global_C_global_F")[1]
disp=.1

for (a0 in c_seq) {
  
  for (scena_ID in 1:length(name_scena)){ #for each scenario of species pairs
    
    for (S in S_seq) { #varying dispersal scale
      
      for (n_ini in 1:n_random_ini){ # for each combination, we draw random communities
        
        type_random=ifelse(n_ini<11,"random_D","random_R")
        
        state =Get_PA_initial_state(Get_MF_initial_state(type=type_random)[-4])
        julia_assign("state", state)
        param=Get_PA_parameters()
        param["cintra"]=.3
        param["S"] = S
        param["alpha_0"] = a0
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
        
        d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = a0,n_ini=n_ini,Scena=name_scena[scena_ID],Delta=disp))
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




tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.4,length.out=N_sim)
name_scena=c("local_C_local_F","global_C_global_F","local_C_global_F","global_C_local_F")
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
c_seq=seq(0,.4,length.out=N_sim)
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
  
  if (scale=="local"){
    d2t=d_state%>%
      filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
      mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6))%>%
      mutate(., Scena=recode_factor(Scena,"local_C_global_F"="Global facilitation","local_C_local_F"="Local facilitation"))
    
    
  } else {
    d2t=d_state%>%
      filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
      mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6))%>%
      mutate(., Scena=recode_factor(Scena,"global_C_global_F"="Global facilitation","global_C_local_F"="Local facilitation"))
  }

  color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
  
  p=ggplot(d2t%>%
             mutate(all_state=recode_factor(all_state,
                                            "Desert/Coexistence"="Coexistence/Desert",
                                            "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                            "Desert/Stress_tolerant"="Stress_tolerant/Desert",
                                            "Competitive/Stress_tolerant"="Competitive/Stress-tolerant"))) +
    geom_tile(aes(x=Stress,y=alpha_0,fill=all_state))+
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'), fill = "") +
    theme(legend.text = element_text(size = 11))+
    scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                               "Competitive/Stress-tolerant"="#C998CE",
                               "Competitive" = "#ACD87B", 
                               "Competitive/Coexistence" = "#DDEFCA",
                               "Stress_tolerant" = "#7BD8D3",
                               "Stress_tolerant/Desert" ="#0F8E87",
                               "Coexistence/Stress_tolerant"="#9BBBB9",
                               "Coexistence/Desert"="#C19E5E",
                               "Desert"=  "#696969"
                               ))+
    facet_grid(Scena~Delta,labeller=label_bquote(cols = delta == .(Delta)))+
    the_theme+theme(strip.text.x = element_text(size=12))
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Multistability_",scale,"_competition.pdf"),p,width = 7,height = 6)
  
  
  if (scale=="local"){
    d2t=d2%>%
      filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
      mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6))%>%
      mutate(., Scena=recode_factor(Scena,"local_C_global_F"="Global facilitation","local_C_local_F"="Local facilitation"))
    
    
  } else {
    d2t=d2%>%
      filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
      mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6))%>%
      mutate(., Scena=recode_factor(Scena,"global_C_global_F"="Global facilitation","global_C_local_F"="Local facilitation"))
  }
  
  
  
  p=ggplot(d2t%>%filter(., alpha_0 %in% c(0,.4,unique(d2$alpha_0)[50]))
           %>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
    geom_line(aes(x=Stress,y=value,linetype=Branches,color=variable))+
    facet_grid(Scena+Delta~alpha_0,labeller=label_bquote(rows = delta == .(Delta),cols = alpha[e] == .(alpha_0)))+
    the_theme+theme(strip.text.x = element_text(size=12))+
    scale_color_manual(values=as.character(color_rho[c(4,2)]))+
    labs(x="Stress (S)",y="Densities",linetype="",color="")
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Bifu_bistability_species_",scale,"_competition.pdf"),width = 8,height = 7)
  
  p=ggplot(d2t%>%filter(., alpha_0 %in% c(0,.4,unique(d2$alpha_0)[50])))+
    geom_line(aes(x=Stress,y=Rho_plus,linetype=Branches))+
    facet_grid(Scena+Delta~alpha_0,labeller=labeller(Delta=label_parsed,alpha_0=label_parsed))+
    the_theme+theme(strip.text.x = element_text(size=12))+
    labs(x="Stress (S)",y="Densities",linetype="",color="")
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Bifu_bistability_community_",scale,"_competition.pdf"),width = 8,height = 7)
  
  set.seed(123)
  v=runif(2)
  d2t$CSI=v[1]*d2t$Stress_tolerant+v[2]*d2t$Competitive
  
  
  p=ggplot(d2t%>%filter(., alpha_0 %in% c(0,.4,unique(d2$alpha_0)[50])))+
    geom_point(aes(x=Stress,y=CSI),size=.5,alpha=.5,shape=21)+
    facet_grid(Scena+Delta~alpha_0,labeller=labeller(Delta=label_parsed,alpha_0=label_parsed))+
    the_theme+theme(strip.text.x = element_text(size=12))+
    labs(x="Stress (S)",y="Community index",linetype="",color="")
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/CSI_",scale,"_competition.pdf"),width = 8,height = 7)
  
  
  #Hysteresis size
  
  d_hysteresis=tibble()
  for (a0 in unique(d2t$alpha_0)){
    for (disp in unique(d2t$Delta)){
      for (scena in unique(d2t$Scena)){
        
        d_fil=filter(d2t,alpha_0==a0,Delta==disp,Scena==scena)
        d_fil=d_fil[order(d_fil$Branches,d_fil$Stress),]
        d_hysteresis=rbind(d_hysteresis,tibble(Competition=a0,Delta=disp,Scena=scena,
                                               hysteresis_com_not_scaled = abs(sum(d_fil$Rho_plus[1:100]-d_fil$Rho_plus[101:200])),
                                               hysteresis_com_scaled = abs(sum(d_fil$Rho_plus[1:100]-d_fil$Rho_plus[101:200]))/(sum(d_fil$Rho_plus[1:100]+d_fil$Rho_plus[101:200])/2),
                                               hysteresis_STol_not_scaled = abs(sum(d_fil$Stress_tolerant[1:100]-d_fil$Stress_tolerant[101:200])),
                                               hysteresis_STol_scaled = abs(sum(d_fil$Stress_tolerant[1:100]-d_fil$Stress_tolerant[101:200]))/(sum(d_fil$Stress_tolerant[1:100]+d_fil$Stress_tolerant[101:200])/2),
                                               hysteresis_Comp_not_scaled = abs(sum(d_fil$Competitive[1:100]-d_fil$Competitive[101:200])),
                                               hysteresis_Comp_scaled = abs(sum(d_fil$Competitive[1:100]-d_fil$Competitive[101:200]))/(sum(d_fil$Competitive[1:100]+d_fil$Competitive[101:200])/2)
                                               ))
      }
    }
  }
  

  p=ggplot(d_hysteresis%>%melt(., measure.vars=c("hysteresis_STol_scaled","hysteresis_Comp_scaled","hysteresis_com_scaled"))%>%
             mutate(., variable=recode_factor(variable,"hysteresis_STol_scaled"="Stress-tolerant","hysteresis_Comp_scaled"="Competitive",
                                              "hysteresis_com_scaled"="Community" )))+
    geom_point(aes(x=Competition,y=value,color=variable),size=1.25,alpha=.5)+
    facet_grid(Scena~Delta,labeller=label_bquote(cols = delta == .(Delta)),scales = "free")+
    the_theme+theme(strip.text.x = element_text(size=12))+
    labs(x=TeX("$\\alpha_e$"),y="Hysteresis size",linetype="",color="")+
    scale_color_manual(values=c(as.character(color_rho[c(4,2)]),"black"))
  
  
  ggsave(paste0("../Figures/2_species/PA/Multistability/Fixed_traits/Hysteresis_",scale,"_competition.pdf"),width = 7,height = 5)
  
  
}


## 4) Scale and strength of facilitation, dispersal and competition on niche expansion ----
### a) Do the simulation ----



tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 100
N_rep2 = 10
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(0, .3, length.out = N_rep2)
f_seq=seq(0,1,length.out=N_rep2)
delta_seq=c(.1, .9)


name_scena=c("global_C_local_F","global_C_global_F")

d_niche=d_RNE=d_all_dyn=tibble() #initializing the tibble

for (scena_ID in 1:2){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    

      
    for (f in f_seq){ #varying facilitation strength
      
      for (alpha0 in alpha_seq) {
        
        
        #Setting the parameters
        param=Get_PA_parameters()
        
        param["cintra"]=.3
        param["f"]=f
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
          Facilitation = f, Local_disp = disp, H = 1, Scena = name_scena[scena_ID], alpha_0 =alpha0,
          Delta_niche_1 = 100 * (S_critic1_both - S_critic1_alone) / S_critic1_alone,  #making it a percentage of initial niche
          Delta_niche_2 = 100 * (S_critic2_both - S_critic2_alone) / S_critic2_alone)) #making it a percentage of initial niche
        
        
        
        d3$rho_plus = d3$rho_1 + d3$rho_2
        
        
        d_compe=filter(d3,Sp=="Competitive")
        d_stesstol=filter(d3,Sp=="Stress_tolerant")
        d_both=filter(d3,Sp=="Both")
        
        d_RNE=rbind(d_RNE,tibble(S=d_both$S,alpha_0=d_both$alpha_0, Facilitation=f,H=1,Local_disp=disp,Scena=name_scena[scena_ID], 
                                 RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1),
                                 RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2)))
        
        
        d_all_dyn=rbind(d3,d_all_dyn)
        
        
      } #end competition loop
      
    } #end facilitation loop
    
  } #end dispersal loop
  
} #end scenario loop
write.table(d_niche,"../Table/2_species/PA/Niche_expansion_PA.csv",sep=";")

d_RNE[,6:7][is.na(d_RNE[,6:7])] = NA

write.table(d_RNE,"../Table/2_species/PA/RNE_PA.csv",sep=";")
write.table(d_all_dyn,"../Table/2_species/PA/All_dyn_PA_RNE.csv",sep=";")









### b) Analyse + plot graphics ----

#First, niche of species
d_niche=read.table("../Table/2_species/PA/Niche_expansion_PA.csv",sep=";")
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


## 5) Net effects fixed traits ----
### a) Simulation ----

tspan = c(0, 30000)
t = seq(0, 30000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
julia_assign("state", state)

N_sim = 100;epsilon=10^(-8)
C_for_analyse = c(0,.15,.3)
S_seq = seq(0, 1, length.out = N_sim)
disp_seq=c(.1,.9)

d_niche=d_RNE=d_all_dyn=tibble() 
d2 = d3 = tibble()
param=Get_PA_parameters()  
param=c(param[1:3],"beta1"=1,"beta2"=1,param[5:12])





  
for (comp in C_for_analyse) {
  
  for (disp in disp_seq){
    
    for (S in S_seq) {
      
      for (sp in 1:2) { #for both species
        
        param["S"] = S
        param["alpha_0"] = comp # update parameters
        param["delta"] = disp
        julia_assign("p", param)
        
        
        # first without press perturbation
        
        prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
        
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-3000):nrow(d),c(1,2)[-sp]])) %>%
                     add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],
                                Disp=disp)) 
        d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, alpha_0 = comp))
        
        # with press perturbation on growth rate proxy
        
        param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
        julia_assign("p", param)
        
        
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-3000):nrow(d),c(1,2)[-sp]])) %>%
                     add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],
                                Disp=disp))
        
        param[paste0("beta",sp)] = 1
      }
    }
  }
}

  
d2[d2 < epsilon] = 0
colnames(d2) = c("Eq", "S", "alpha_0", "Type","Species","Disp")


write.table(d2,"../Table/2_species/PA/Net_effects_partial_deriv_PA.csv",sep=";")



### b) Analysis ----

d2=read.table("../Table/2_species/PA/Net_effects_partial_deriv_PA.csv",sep=";")




net_effect =sapply(seq(1, nrow(d2) , by = 2),function(x){
  return((d2$Eq[x+1] - d2$Eq[x]) / epsilon)
})


d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect


for (id_plot in 1:3){
  assign(paste0("p_",id_plot),
         ggplot(d_net %>%filter(round(alpha_0,4)==round(C_for_analyse[id_plot]))) +
           geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
           the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
           geom_hline(yintercept = 0,linetype=9)+
           facet_wrap(.~Disp,labeller=label_bquote(cols = delta == .(Disp)))+
           labs(x="Stress (S)",alpha=TeX('$\\alpha_e$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial \beta_j}$)'))  )
}

p=ggarrange(p_1,p_2,p_3,ncol=3,common.legend = T,legend = "bottom")


ggsave(paste0("../Figures/2_species/PA/Net_effects_PA.pdf"),p_tot,width = 8,height = 5)


## 6) Neighboring effects fixed traits ----
### a) Simulation ----

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
Nsim=50
S_seq = seq(0, 1, length.out = Nsim)
c_seq = seq(0,.3, length.out = Nsim)
epsilon=10^(-9)
disp_seq=c(.1,.9)
scale_seq=c("LocalF_GlobalC","LocalF_LocalC","GlobalF_GlobalC","GlobalF_LocalC")

d_RNE=tibble()


for (scale in scale_seq){
  
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
        
        if (scale == "LocalF_GlobalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        if (scale == "LocalF_LocalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }
        if (scale == "GlobalF_GlobalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
          
        }
        if (scale == "GlobalF_LocalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
          
        }
        
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d_RNE = rbind(d_RNE, d[nrow(d),c(1,2) ] %>% add_column(S = S,
                                                               alpha_0 = comp,
                                                               Sp="Competitive",
                                                               Scale_comp=scale,
                                                               Dispersal=disp))
        
        
        
        #only stress tolerant one
        state =Get_PA_initial_state(c(0.8,0,.1,.1))
        
        julia_assign("state", state)
        
        if (scale == "LocalF_GlobalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        if (scale == "LocalF_LocalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }
        if (scale == "GlobalF_GlobalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
          
        }
        if (scale == "GlobalF_LocalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
          
        }
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d_RNE = rbind(d_RNE, d[nrow(d),c(1,2) ] %>% add_column(S = S,
                                                               alpha_0 = comp,
                                                               Sp="Stress-tolerant",
                                                               Scale_comp=scale,
                                                               Dispersal=disp))


        
        
        #both species
        
        state =Get_PA_initial_state(c(0.4,.4,.1,.1))
        julia_assign("state", state)
        
        if (scale == "LocalF_GlobalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
        }
        if (scale == "LocalF_LocalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
          
        }
        if (scale == "GlobalF_GlobalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
          
        }
        if (scale == "GlobalF_LocalC"){ 
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
          
        }
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
        
        d_RNE = rbind(d_RNE, d[nrow(d),c(1,2) ] %>% add_column(S = S, 
                                                         alpha_0 = comp,
                                                         Sp="Both",
                                                         Scale_comp=scale,
                                                         Dispersal=disp))
        
        
      }
    }
  }
}
d_RNE[d_RNE<10^(-4)]=0

d_RNE$rho_plus = d_RNE$rho_1 + d_RNE$rho_2

d_compe=filter(d_RNE,Sp=="Competitive")
d_stesstol=filter(d_RNE,Sp=="Stress-tolerant")
d_both=filter(d_RNE,Sp=="Both")

d_RNE=tibble(S=d_both$S,alpha_0=d_both$alpha_0, 
             RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1), #actually its RII
             RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2),
             NintA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
             NintA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)),
             Scale=d_both$Scale_comp,Disp=d_both$Dispersal)

d_RNE[,3:6][is.na(d_RNE[,3:6])] = NA



d_all=d_RNE
save(file="../Table/2_species/PA/Comparing_net_effects.RData",d_all)

### b) Analysis ----
load(file="../Table/2_species/PA/Comparing_net_effects.RData")
d_RNE=d_all



d_RNE[,3:6][d_RNE[,3:6] > 2]=NA
d_RNE[,3:6][d_RNE[,3:6] < -1]=NA


p1=ggplot(d_RNE%>%
            melt(., measure.vars=c("RNE_stress_tol","RNE_competitive"))%>%
            mutate(., variable=recode_factor(variable,"RNE_stress_tol"='Stress-tolerant',"RNE_competitive"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_grid(Scale~Disp+variable)+
  the_theme+labs(x="",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  theme(axis.text.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank(),
        strip.text.x = element_text(size=12))
ggsave("../Figures/2_species/PA/RII_varying_scales_disp.pdf",p1,width = 7,height = 7)

p2=ggplot(d_RNE%>%
            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
            #mutate(., value=rescale(value,to=c(-1,1)))%>%
            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_grid(Scale~Disp+variable)+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")

ggsave("../Figures/2_species/PA/NIntA_varying_scales_disp.pdf",p2,width = 7,height = 7)





## 7) Multistability and trait difference ----
### a) Simulation ----
tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=100
S_seq = seq(0,1, length.out = N_sim)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=c(0,.1, .2, .3,.4)[5]
psi1_seq=c(1,0)
f_seq=c(0,.45,.9)
dispersal_scale=c(.1,.9)
branches=c("Degradation","Restoration")

for (facil in f_seq){

  for (disp in dispersal_scale){
    
  
    for (psi1 in psi1_seq){
    
      for (branch in branches){
        
        if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
          state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
          S_seq=seq(0,1, length.out = N_sim)
        }else {
          state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
          S_seq=rev(seq(0,1, length.out = N_sim))
        }
        
        for (cinter in c_inter_seq){
          d2 = tibble()
          
        
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
              param["f"]=facil
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
                                "_Psi1_",psi1,"_delta_",disp,"_facilitation_",facil,".csv"),sep=";")
          
        } # end loop interspecific competition
        
      } # end loop branch
      
    } # end loop first species trait
    
  } #end loop dispersal

}#end facilitation loop


### b) Analysis ----

c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)
dispersal_scale=c(.1,.9)
f_seq=c(0,.45,.9)

for (f in f_seq){

  for (disp in dispersal_scale){
  
    for (Psi_sp1 in psi1_seq){
      d=tibble()  
      
      for (branch in c("Restoration","Degradation")){
        for (cinter in c_inter_seq){
          
          d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                               cinter,"_branch_",branch,
                               "_Psi1_",Psi_sp1,"_delta_",disp,"_facilitation_",f,".csv"),sep=";")
          d=rbind(d,d2)
          
        } # end loop interspecific competition
      } # end loop branch
      
      
      
      
      d[,1:2][d[,1:2] < 10^-4] = 0
      
      
      #COmputing CSI index
      set.seed(123)
      u=runif(2)
      d$CSI = sapply(1:nrow(d),function(x){
        return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
      })
      
      
      
      
      if (Psi_sp1 ==1) {
        
        d2=filter(d,Psi1==Psi_sp1)
        
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
                                     "Species 2/Stress_tolerant" = "#C998CE"),
                            labels=c("Coexistence","Species 2","Coexistence/Species 2","Stress-tolerant","Stress-tolerant/Desert",
                                     "Coexistence/Stress-tolerant",
                                     "Coexistence/Desert","Desert","Species 2/Stress-tolerant"))+
          facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
          the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Multistability_trait_variation_sp_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 9,height = 4)
        
        
        
        d2t=transform(d2,
                      alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                              "alpha[e] == 0.5")),
                      Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                                  labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
        
        
        
        
        #Bifurcation diagram species
        p=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                   filter(., Psi2 %in% unique(d_state$Psi2)[c(2,75,100)])%>%
                   mutate(Psi2=round(abs(Psi2),3))%>%
                   filter(., alpha_0!=.1))+
          geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
          labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
          the_theme +
          scale_color_manual(values = color_rho[c(2, 4)],labels=c("Species 2","Stress-tolerant")) +
          scale_linetype_manual(values=c(1,9))+
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = psi[2] == .(Psi2)))
        
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Bifurcation_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        

        #Community index
        p=ggplot(d2%>%
                   filter(., Psi2 %in% unique(d2$Psi2)[c(2,75,100)])%>%
                   mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                   filter(., alpha_0!=.1))+
          geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
          labs(x = "Stress (S)", y = "Community index", color = "",linetype="") +
          the_theme +
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/CSI_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        
        
        
        #Restoration and degradation points
        d_tipping=tibble()
        for (a0 in unique(d2$alpha_0)){
          for (psi in unique(d2$Psi2)){
            for (direction in unique(d2$Branches)){
              
              d_a0=filter(d2,alpha_0==a0,Psi2==psi,Branches==direction)
              
              if (direction=="Restoration") d_a0=d_a0[order(d_a0$Stress,decreasing = T),]
              
              density_C=d_a0$Competitive
              stress_C=stress_ST=d_a0$Stress
                
              if (direction=="Restoration"){ #initially absent
                
                stress_C=stress_C[-c(1:min(which(abs(diff(density_C))>0)))]
                density_C=density_C[-c(1:min(which(abs(diff(density_C))>0)))]
                Tipping_C = stress_C[1]
                
              } else{
                
                Tipping_C = stress_C[min(which(density_C==0))-1]
                
              }
              
              d_tipping=rbind(d_tipping,tibble(Competition=a0,Psi2=psi,Branch=direction,
                                               Tipping_C = Tipping_C))
            
            }
          }
        }
        
        p=ggplot(d_tipping)+
          geom_smooth(aes(x=Psi2,y=Tipping_C,color=Branch,group=interaction(Competition,Branch)),se = F)+
          the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
          labs(x=TeX(r'(Trait species 2, \ $\psi_2)'),y="Tipping point",color="")+
          scale_color_manual(values=c("black","blue"))+
          scale_x_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                                 name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)'),
                                                 breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1")),
                             breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1"))
          
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Tipping_points_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)        
        
        
        
        
        
      } 
      if (Psi_sp1==0){
        
        
        d2=filter(d,Psi1==Psi_sp1)
        
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
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Multistability_trait_variation_sp_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 9,height = 4)
        
        
        
        #Bifurcation diagrams
        
        d2t=transform(d2,
                      alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                              "alpha[e] == 0.5")),
                      Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                                  labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
        
        
        p=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                   filter(., Psi2 %in% unique(d_state$Psi2)[c(2,75,100)])%>%
                   mutate(Psi2=round(Psi2,3))%>%
                   filter(., alpha_0!=.1))+
          geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
          labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
          the_theme +
          scale_color_manual(values = (as.character(color_rho[c(2, 4)])),labels=c("Competitive","Species 2")) +
          scale_linetype_manual(values=c(1,9))+
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = psi[1] == .(Psi2)))
        
        
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Bifurcation_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        
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
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/CSI_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        

        
        
        #Restoration and degradation points
        d_tipping=tibble()
        for (a0 in unique(d2$alpha_0)){
          for (psi in unique(d2$Psi2)){
            for (direction in unique(d2$Branches)){
              
              d_a0=filter(d2,alpha_0==a0,Psi2==psi,Branches==direction)
              
              if (direction=="Restoration") d_a0=d_a0[order(d_a0$Stress,decreasing = T),]
              
              density_ST=d_a0$Stress_tolerant
              stress_ST=stress_ST=d_a0$Stress
              
              if (direction=="Restoration"){ #initially absent
                
                stress_ST=stress_ST[-c(1:min(which(abs(diff(density_ST))>0)))]
                density_ST=density_ST[-c(1:min(which(abs(diff(density_ST))>0)))]
                Tipping_ST = stress_ST[1]
                
              } else{
                
                Tipping_ST = stress_ST[min(which(density_ST==0))-1]
                
              }
              
              d_tipping=rbind(d_tipping,tibble(Competition=a0,Psi2=psi,Branch=direction,
                                               Tipping_ST = Tipping_ST))
              
            }
          }
        }
        
        p=ggplot(d_tipping)+
          geom_smooth(aes(x=Psi2,y=Tipping_ST,color=Branch,group=interaction(Competition,Branch)),se = F)+
          the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
          labs(x=TeX(r'(Trait species 1, \ $\psi_1)'),y="Tipping point",color="")+
          scale_color_manual(values=c("black","blue"))+
          scale_x_continuous(sec.axis = sec_axis(trans = ~ (.x) ,
                                                 name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)'),
                                                 breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1")),
                             breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1"))
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Tipping_points_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)        
        
        
        
      }
    }
  }
}


## 8) Clustering between species fixed traits  ----
### a) Simulation ----


tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 50
S_seq = c(0,.1,.73,.77)
alpha_seq = c(.2)
f_seq=.9
delta_seq=seq(0,1,length.out=N_rep)
cintra_seq=c(.3)


name_scena=c("global_C_global_F","global_C_local_F")

d_clustering=tibble() #initializing the tibble

for (scena_ID in 1:2){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (aii in cintra_seq){ #varying intraspecific competition strength
      
      for (f in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) { #varying competition
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          param["r"]=0
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
            
            if (scena_ID==1){  #global C, global F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }else if (scena_ID==2){ #global C, local F
              
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
            Rho_22=d2$rho_22,Rho_11=d2$rho_11,Rho_10=d2$rho_1 - d2$rho_11 - d2$rho_12 - d2$rho_1m,
            Rho_20=d2$rho_2 - d2$rho_22 - d2$rho_12 - d2$rho_2m,
            Rho_0=1-d2$rho_1-d2$rho_1-d2$rho_m,
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


d_clustering$c12=d_clustering$c11=d_clustering$c22=0
name_scena=c("global_C_global_F","global_C_local_F")

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
    q1_1 = rho_11 / rho_1 #average number of species 1 around species 1 
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



d_clustering=mutate(d_clustering, Scena=recode_factor(Scena,"global_C_global_F"="Global facilitation","global_C_local_F"="Local facilitation"))



p=ggplot(d_clustering%>%
               filter(., alpha_0==.2)%>%
               group_by(., cintra,alpha_0,S,delta,Scena)%>%
               summarise(.,.groups ="keep",c11=mean(c11) ))+
  geom_point(aes(x=as.numeric(delta),y=c11),size=1.5,alpha=.7,shape=1)+
  facet_grid(Scena~S,labeller=label_bquote(cols = Stress == .(S)),scales = "free")+
  the_theme+labs(x=TeX("$\\delta$"),y=expression(paste("Stress-tolerant clustering (c"[11],")")))+
  scale_color_manual(values=as.character(color_rho[c(2,4)]))+
  geom_hline(yintercept = 1)




ggsave("../Figures/2_species/PA/Clustering_11_dispersal.pdf",p,width = 7,height = 4)






## 9) Niche species, varying traits ----
### a) Simulation ----

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 100
S_seq = seq(0,1,length.out=N_rep)
alpha_seq = seq(0, .3, length.out = 3)
f_seq=c(.3,.9)
trait_sp2_seq=seq(0,1,length.out=50)
trait1_seq=c(0,1)



d_niche=tibble() #initializing the tibble

for (psi_sp1 in trait1_seq){ #varying the trait of sp1
  
  for (f in f_seq){ #varying facilitation strength
    
    for (alpha0 in alpha_seq) { #varying competition
      
      for (trait_2 in trait_sp2_seq){
        
        #Setting the parameters
        param=Get_PA_parameters()
        
        param["f"]=f
        param["alpha_0"]=alpha0
        param["psi_1"]=psi_sp1 #we fixed to stress-tolerant species 
        param["psi_2"]=trait_2 #and vary the other species
        
        
        
        # 1) Competitive species alone
        
        state=Get_PA_initial_state(c(0,.8,.1))
        julia_assign("state", state)
        d2 = tibble()
        
        for (S in S_seq) {
          
          param["S"] = S
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Competitive"))
          
          
        }
        d2[d2 < 10^-4] = 0
        
        
        S_critic2_alone=abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) #range of values where competitive species is
        
        
        
        # 2) Stress-tolerant species alone
        
        state=Get_PA_initial_state(c(0.8,0,.1))
        julia_assign("state", state)
        
        d2 = tibble()
        S_seq = seq(0, 1, length.out = 100)
        
        for (S in S_seq) {
          
          param["S"] = S
          julia_assign("p", param)
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S,Sp="Stress_tolerant"))
          
        }
        
        d2[d2 < 10^-4] = 0
        
        S_critic1_alone=abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
        
        
        
        #3) Coexisting species
        
        state=Get_PA_initial_state(c(0.4,.4,.1))
        julia_assign("state", state)
        
        
        d2 = tibble()
        
        for (S in S_seq) { #varying the stress 
          
          param["S"] = S
          julia_assign("p", param)
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, Sp="Both"))
          
        }
        d2[d2 < 10^-4] = 0
        d2$rho_plus = d2$rho_1 + d2$rho_2
        
        S_critic1_both = abs(diff(range(d2$S[which(d2$rho_1 !=0 )]))) 
        S_critic2_both = abs(diff(range(d2$S[which(d2$rho_2 !=0 )]))) 
        
        
        d_niche=rbind(d_niche,tibble(
          Facilitation = f, alpha_0 = alpha0, Psi1 = psi_sp1, Psi2 = trait_2,
          Delta_niche_1 = 100 * (S_critic1_both - S_critic1_alone) / S_critic1_alone,  #making it a percentage of initial niche
          Delta_niche_2 = 100 * (S_critic2_both - S_critic2_alone) / S_critic2_alone)) #making it a percentage of initial niche
        
        
      } #end varying trait loop
      
    } #end competition loop
    
  } #end facilitation loop
  
} #end sp1 trait loop
write.table(d_niche,"../Table/2_species/PA/Niche_expansion_varying_traits.csv",sep=";")




### b) Analysis ----


d_niche = read.table(paste0("../Table/2_species/PA/Niche_expansion_varying_traits.csv"),sep=";")
d_niche = do.call(data.frame,lapply(d_niche,function(x) replace(x, is.infinite(x), NA)))
d_niche$Delta_niche_1[102]=NA


p=ggplot(NULL)+
  geom_path(data=d_niche%>%
              melt(., measure.vars=c("Delta_niche_1","Delta_niche_2"))%>%
              mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
              mutate(., Delta_psi=abs(Psi1-Psi2),
                     Facilitation=as.factor(Facilitation),
                     alpha_0=as.factor(alpha_0)),
            aes(x=Delta_psi,y=value,group=interaction(alpha_0,Facilitation),color=alpha_0),lwd=1)+
  geom_point(data=d_niche%>%
               melt(., measure.vars=c("Delta_niche_1","Delta_niche_2"))%>%
               mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
               mutate(., Delta_psi=abs(Psi1-Psi2),
                      Facilitation=as.factor(Facilitation),
                      alpha_0=as.factor(alpha_0))%>%
               filter(., Psi2 %in% unique(.$Psi2)[seq(1,length(unique(.$Psi2)),by =5)]), #to only have 10 points
             aes(x=Delta_psi,y=value,shape=Facilitation,color=alpha_0),fill="white",size=2)+
  facet_grid(variable~Psi1,scales="free",labeller = label_bquote(cols = psi[1] == .(Psi1)))+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=c("#91D4DE","#26A8D2","#06247B"))+
  labs(y = "Niche change (%)", x = TeX(r'(Trait difference \ |$\psi_1-\psi_2|)'),
       color=TeX(r'(Competition strength \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f)')) +
  the_theme

ggsave("../Figures/2_species/PA/Niche_expansion_varying_traits.pdf",p,width = 7,height = 6)






## 10) Neighboring effects, varying traits ----
### a) Simulation ----



tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 50
S_seq = seq(0,1,length.out=N_rep)
alpha_seq = seq(0, .3, length.out = 3)
f_seq=c(.9)
trait_sp2_seq=seq(0,1,length.out=50)
trait1_seq=c(0,1)


d_RII=tibble() #initializing the tibble

for (psi_sp1 in trait1_seq){ #varying the trait of sp1
  
  for (f in f_seq){ #varying facilitation strength
    
    for (alpha0 in alpha_seq) { #varying competition
      
      for (trait_2 in trait_sp2_seq){
        
        #Setting the parameters
        param=Get_PA_parameters()
        
        param["f"]=f
        param["alpha_0"]=alpha0
        param["psi_1"]=psi_sp1 #we fixed to stress-tolerant species 
        param["psi_2"]=trait_2 #and vary the other species
        
        
        # 1) Competitive species alone
      
        state=Get_PA_initial_state(c(0,.8,.1))
        
        julia_assign("state", state)
        d2 = tibble()
        
        for (S in S_seq) {
          
          param["S"] = S
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d_RII = rbind(d_RII,d[nrow(d), ] %>% add_column(S = S,Sp="Competitive",alpha_0=alpha0,Psi1=psi_sp1,Psi2=trait_2,Facilitation=f))
          
        }
        
        
        
        # 2) Stress-tolerant species alone
        
        
        state=Get_PA_initial_state(c(.8,0,.1))
        julia_assign("state", state)
        
        d2 = tibble()
        
        for (S in S_seq) {
          
          param["S"] = S
          julia_assign("p", param)
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d_RII = rbind(d_RII,d[nrow(d), ] %>% add_column(S = S,Sp="Stress_tolerant",alpha_0=alpha0,Psi1=psi_sp1,Psi2=trait_2,Facilitation=f))
          
        }
        
        
        
        
        #3) Coexisting species
        
        state=Get_PA_initial_state(c(.4,.4,.1))
        julia_assign("state", state)
        
        d2 = tibble()
        
        for (S in S_seq) { #varying the stress 
          
          param["S"] = S
          julia_assign("p", param)
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d_RII = rbind(d_RII,d[nrow(d), ] %>% add_column(S = S,Sp="Both",alpha_0=alpha0,Psi1=psi_sp1,Psi2=trait_2,Facilitation=f))
          
        }
        d_RII[d_RII < 10^-4] = 0
        
        
      } #end varying trait loop
      
    } #end competition loop
    
  } #end facilitation loop
  
} #end sp1 trait loop


write.table(d_RII,"../Table/2_species/PA/RII_varying_traits.csv",sep=";")


 
### b) Analysis ----
d_RII=read.table(paste0("../Table/2_species/PA/RII_varying_traits.csv"),sep=";")

d_compe=filter(d_RII,Sp=="Competitive")
d_stesstol=filter(d_RII,Sp=="Stress_tolerant")
d_both=filter(d_RII,Sp=="Both")

d_RII=tibble(S=d_both$S,alpha_0=d_both$alpha_0, Psi1=d_both$Psi1,Psi2=d_both$Psi2,Facilitation=d_both$Facilitation,
             RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1), #acutally its RII
             RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2),
             NintA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
             NintA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)))
d_RII[,3:6][is.na(d_RII[,3:6])] = NA


p1=ggplot(d_RII %>%filter(., Psi1==1)%>%
            melt(., measure.vars=c("NintA_st","NintA_comp")) %>%
            mutate(., variable = recode_factor(variable, "NintA_st" = "Species 1", "NintA_comp" = "Species 2"))) +
  geom_tile(aes(x=round(S,5),y=round(Psi2,5),fill=as.numeric(value)))+
  facet_grid(variable~alpha_0,labeller = label_bquote(cols= alpha[e] == .(alpha_0)))+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  the_theme+ggtitle(TeX("$\\psi_1 = 1$"))+
  labs(x="Stress, S",y=TeX(r'(Trait species 2, \ $\psi_2)'),fill=TeX("$NInt_A$"))

p2=ggplot(d_RII %>%filter(., Psi1==0)%>%
            melt(., measure.vars=c("NintA_st","NintA_comp")) %>%
            mutate(., variable = recode_factor(variable, "NintA_st" = "Species 2", "NintA_comp" = "Species 1"))) +
  geom_tile(aes(x=round(S,5),y=round(Psi2,5),fill=as.numeric(value)))+
  facet_grid(variable~alpha_0,labeller = label_bquote(cols= alpha[e] == .(alpha_0)))+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  the_theme+ggtitle(TeX("$\\psi_2 = 0$"))+
  labs(x="Stress, S",y=TeX(r'(Trait species 1, \ $\psi_1)'),fill=TeX("$NInt_A$"))

p_tot=ggarrange(p1+theme(axis.title.x = element_blank(),axis.line.x = element_blank(),axis.ticks.x = element_blank()),
                p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
                nrow=2,labels = letters[1:2],common.legend = T,legend = "bottom")

ggsave("../Figures/2_species/PA/NInt_A_varying_traits.pdf",p_tot,width = 7,height = 8)

  
## 11) Net effects varying traits ----
### a) Simulation ----


tspan = c(0, 20000)
t = seq(0, 20000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
julia_assign("state", state)

N_sim = 50;epsilon=10^(-6)
psi2_seq= c(.1,.3,.5,.7,.9)
psi1_seq=c(0,1)
S_seq = seq(0, 1, length.out = N_sim)
C_seq=c(.1,.2,.3)

d2 = d3 = tibble()




for (psi1 in psi1_seq){
  
  for (comp in C_seq) {
    
    for (psi2 in psi2_seq){
      
      for (S in S_seq) {
        
        for (sp in 1:2) { #for both species
          
          param=Get_PA_parameters()  
          param["S"] = S
          param["alpha_0"] = comp # update parameters
          param["psi_2"] = psi2
          param["psi_1"] = psi1
          julia_assign("p", param)
          
          
          # first without press perturbation
          
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
          
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-2500):nrow(d),c(1,2)[sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],
                                  Psi1=psi1,Psi2=psi2)) 
          d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, alpha_0 = comp,Psi1=psi1,Psi2=psi2))
          
          # with press perturbation on growth rate proxy
          
          param["psi_2"]=param["psi_2"]+epsilon
          julia_assign("p", param)
          
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_varying_trait_press, state, tspan, p)")
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-2500):nrow(d),c(1,2)[sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],
                                  Psi1=psi1,Psi2=psi2))
          
        }
      }
    }
  }
}


d2[d2 < epsilon] = 0
colnames(d2) = c("Eq", "S", "alpha_0", "Type","Species","Psi1","Psi2")

write.table(d2,"../Table/2_species/PA/Net_effect_varying_traits.csv",sep=";")



### b) Analysis ----


epsilon=10^(-6)
d2=read.table("../Table/2_species/PA/Net_effect_varying_traits.csv",sep=";")


net_effect =sapply(seq(1, nrow(d2) , by = 2),function(x){
  return((d2$Eq[x+1] - d2$Eq[x]) / epsilon)
})


d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect


p1=ggplot(d_net%>%
            filter(.,Psi1==1,alpha_0%in% c(.1,.3),Species==1)%>%
            mutate(.,Species=recode_factor(Species,"1" ="Species 1","2" = "Species 2")))+
  ggtitle(TeX("$\\psi_1=1$"))+
  geom_line(aes(x=S,y=value,color=Psi2,group=interaction(Psi1,Psi2,alpha_0,Species)))+
  facet_grid(.~alpha_0,scales = "free",labeller = label_bquote(cols = alpha[e]==.(alpha_0)))+
  labs(x="Stress, S",y=TeX(r'($\frac{\partial \rho_{+_1}}{\partial \psi_2})'),
       color=TeX("$\\psi_2 \ \ $"))+
  the_theme+
  scale_color_stepsn(colours=color_Nsp(6),breaks=seq(0,1,by=.2))+
  geom_hline(yintercept = 0,linetype=9)

p2=ggplot(d_net%>%
            filter(.,Psi1==0,Species==2,alpha_0 %in% c(.1,.2))%>%
            mutate(.,Species=recode_factor(Species,"1" ="Species 1","2" = "Species 2")))+ #here we inverse to be coherent with the analysis made 
  ggtitle(TeX("$\\psi_2=0$"))+
  geom_line(aes(x=S,y=value,color=Psi2,group=interaction(Psi1,Psi2,alpha_0,Species)))+
  facet_grid(.~alpha_0,scales = "free",labeller = label_bquote(cols = alpha[e]==.(alpha_0)))+
  labs(x="Stress, S",y=TeX(r'($\frac{\partial \rho_{+_2}}{\partial \psi_1})'),
       color=TeX("$\\psi_1 \ \ $"))+
  the_theme+
  scale_color_stepsn(colours=color_Nsp(6),breaks=seq(0,1,by=.2))+
  geom_hline(yintercept = 0,linetype=9)


p_tot=ggarrange(p1+theme(strip.text.x = element_text(size=13),strip.background.x = element_blank()),
                p2+theme(strip.text.x = element_text(size=13),strip.background.x = element_blank()),nrow=2,labels = letters[1:2])


ggsave("../Figures/2_species/PA/Net_effects_varying_traits.pdf",p_tot,width = 7,height = 8)


## 12) Varying the trade-off shape ----

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=50
S_seq = seq(0,1, length.out = N_sim)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=c(.3)
psi1_seq=c(1,0)
f_seq=c(.9)
dispersal_scale=c(.1)
branches=c("Degradation","Restoration")
shape_trade_off=c(.5,.75,1,1.5,2)

  
for (disp in dispersal_scale){
  
  
  for (psi1 in psi1_seq){
    
    for (branch in branches){
      
      if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
        state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
        S_seq=seq(0,1, length.out = N_sim)
      }else {
        state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
        S_seq=rev(seq(0,1, length.out = N_sim))
      }
      
      for (tradeoff in shape_trade_off){
        d2 = tibble()
        
        
        for (psi2 in psi_seq){
          
          for (S in S_seq) { 
            
            julia_assign("state", state)
            param=Get_PA_parameters()
            param["delta"]=disp
            param["cintra"]=.3
            param["alpha_0"]=.3
            param["S"] = S
            param["psi_1"]=psi1
            param["psi_2"]=psi2
            param["f"]=.9
            param["shape"]=tradeoff
            julia_assign("p", param)
            
            prob = julia_eval("ODEProblem(PA_two_species_varying_trait_trade_off, state, tspan, p)")
            
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d), ] %>% add_column(Stress = S, Psi2 = psi2, Psi1 = psi1,alpha_0=.3,
                                                       Branch=branch,Shape=tradeoff))
            
          }
        } # end trait value 2nd species
        
        d2[d2 < 10^-4] = 0
        d2$rho_plus = d2$rho_1 + d2$rho_2
        d2=d2[,c(1,2,10:16)]
        colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "Psi2","Psi1","alpha_0","Branches","Shape","Rho_plus")
        write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Multistability_varying_tradeoff_",
                              tradeoff,"_branch_",branch,
                              "_Psi1_",psi1,"_delta_",disp,"_facilitation_",.9,".csv"),sep=";")
        
      } # end loop interspecific competition
      
    } # end loop branch
    
  } # end loop first species trait
  
} #end loop dispersal





for (Psi_sp1 in psi1_seq){
  d=tibble()  
  
  for (branch in c("Restoration","Degradation")){
    for (tradeoff in shape_trade_off){
      
      d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Multistability_varying_tradeoff_",
                           tradeoff,"_branch_",branch,
                           "_Psi1_",Psi_sp1,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";")
      colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "Psi2","Psi1","alpha_0","Branches","Shape","Rho_plus")
      d=rbind(d,d2)
      
    } # end loop interspecific competition
  } # end loop branch
  
  
  
  
  d[,1:2][d[,1:2] < 10^-4] = 0
  
  
  #COmputing CSI index
  set.seed(123)
  u=runif(2)
  d$CSI = sapply(1:nrow(d),function(x){
    return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
  })
  
  
  
  
  if (Psi_sp1 ==1) {
    
    d2=filter(d,Psi1==Psi_sp1)
    
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
    
    d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1,d2$Shape),]
    
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
    
    
    
    
    p1=ggplot(d_state%>%
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
                                 "Species 2/Stress_tolerant" = "#C998CE"),
                        labels=c("Coexistence","Species 2","Coexistence/Species 2","Stress-tolerant","Stress-tolerant/Desert",
                                 "Coexistence/Stress-tolerant",
                                 "Coexistence/Desert","Desert","Species 2/Stress-tolerant"))+
      facet_grid(.~Shape,labeller=label_bquote(cols = gamma == .(Shape)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
    
 
    
    
  } 
  if (Psi_sp1==0){
    
    
    d2=filter(d,Psi1==Psi_sp1)
    
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
    
    d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1,d2$Shape),]
    
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
    
    
    
    p2=ggplot(d_state%>%
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
      facet_grid(.~Shape,labeller=label_bquote(cols = gamma == .(Shape)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
    

  }
}
p_tot=ggarrange(p1,p2,labels = letters[1:2],nrow = 2)

ggsave("../Figures/2_species/PA/Varying_trade_off_shape.pdf",p_tot,width = 10,height = 8)

## 13) Which invades ----

rm(list = ls())
source("./Dryland_shift_functions.R")
julia_setup()
de = diffeq_setup()

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=100
S_seq = c(0,.3,.75)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=c(0,.1, .2, .3,.4)[c(3,4,5)][2]
psi1_seq=seq(0,1,length.out=N_sim)
f_seq=c(.9)
dispersal_scale=c(.1)
branches=c("Degradation","Restoration")

for (facil in f_seq){
  
  for (disp in dispersal_scale){
    
    for (S in S_seq) { 
      
      
      for (cinter in c_inter_seq){
        
        
        
        for (branch in branches){
          
          if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
            state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
            S_seq=seq(0,1, length.out = N_sim)
          }else {
            state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
            S_seq=rev(seq(0,1, length.out = N_sim))
          }
          
          d2 = tibble()
          
          
          for (psi2 in psi_seq){
            
            for (psi1 in psi1_seq){
              
              julia_assign("state", state)
              param=Get_PA_parameters()
              param["delta"]=disp
              param["cintra"]=.3
              param["alpha_0"]=cinter
              param["S"] = S
              param["psi_1"]=psi1
              param["psi_2"]=psi2
              param["f"]=facil
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
          write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Invasion/Test_interspe_comp_",
                                cinter,"_branch_",branch,
                                "_stress_",S,"_delta_",disp,"_facilitation_",facil,".csv"),sep=";")
          
        } # end loop interspecific competition
        
      } # end loop branch
      
    } # end loop first species trait
    
  } #end loop dispersal
  
}#end facilitation loop





c_inter_seq=c(.2,.3,.4)
stress_seq=c(0,.3,.75)


d=tibble()  
for (stress in stress_seq){
  
  for (cinter in c_inter_seq){
    for (branch in c("Restoration","Degradation")){
      
      d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Invasion/Test_interspe_comp_",
                           cinter,"_branch_",branch,
                           "_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";")
      d=rbind(d,d2)
      
    }
    
  } # end loop interspecific competition
} # end loop branch




d[,1:2][d[,1:2] < 10^-4] = 0


#COmputing CSI index
set.seed(123)
u=runif(2)
d$CSI = sapply(1:nrow(d),function(x){
  return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
})

d2=d
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

d_state=d_state%>%
  mutate(all_state=recode_factor(all_state,
                                 "Desert/Species 2"="Sp2/Desert",
                                 "Desert/Coexistence"="Coexistence/Desert",
                                 "Stress_tolerant"="Sp1",
                                 "Coexistence/Stress_tolerant"="Coexistence/Sp1",
                                 "Desert/Stress_tolerant"="Sp1/Desert",
                                 "Stress_tolerant/Coexistence"="Coexistence/Sp1",
                                 "Stress_tolerant/Species 2"="Sp1/Sp2",
                                 "Species 2/Stress_tolerant"="Sp1/Sp2",
                                 "Species 2/Coexistence"="Coexistence/Sp2",
                                 "Coexistence/Species 2"="Coexistence/Sp2",
                                 "Species 2" = "Sp2"))

d_state$all_state2=d_state$all_state

for (nr in 1:nrow(d_state)){
  if (d_state$Psi2[nr]>d_state$Psi1[nr]){
    d_state$all_state[nr]=NA
  }
}

color_multistability=c("Coexistence" = "#D8CC7B",
                       "Sp2" = "#ACD87B",
                       "Coexistence/Sp2" = "#DDEFCA",
                       "Sp1" = "#7BD8D3",
                       "Sp1/Desert" ="#0F8E87",
                       "Coexistence/Sp1"="#9BBBB9",
                       "Coexistence/Desert"="#C19E5E",
                       "Desert"=  "#696969",
                       "Sp1/Sp2" = "#C998CE")

p=ggplot(d_state%>%filter(., Stress !=.75)) +
  geom_tile(aes(x=Psi1,y=as.numeric(Psi2),fill=all_state2))+
  theme_classic() +
  theme(legend.position = "bottom") +
  theme(legend.text = element_text(size = 11))+
  the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
  facet_grid(Stress~alpha_0,labeller = label_bquote(cols= alpha[e]==.(alpha_0),rows = Stress==.(Stress) ))+
  scale_fill_manual(values=color_multistability,na.value = "white")+
  labs(x=TeX("$\\psi_1$"),y=TeX("$\\psi_2$"),fill="")
ggsave("../Figures/2_species/PA/Testing_relative_strategies_twofaces.pdf",width = 8,height = 7)






## 14) Fraction gradient with bistability ----

rm(list = ls())
source("./Dryland_shift_functions.R")
julia_setup()
de = diffeq_setup()

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=30
S_seq =seq(0,1,length.out=N_sim)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=rev(c(0,.1, .2, .3,.4))
psi1_seq=seq(0,1,length.out=N_sim)
f_seq=c(.9)
dispersal_scale=c(.1)
branches=c("Degradation","Restoration")

for (facil in f_seq){
  
  for (disp in dispersal_scale){
    
    for (S in S_seq) { 
      
      
      for (cinter in c_inter_seq){
        
        
        
        for (branch in branches){
          
          if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
            state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
            S_seq=seq(0,1, length.out = N_sim)
          }else {
            state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
            S_seq=rev(seq(0,1, length.out = N_sim))
          }
          
          d2 = tibble()
          
          
          for (psi2 in psi_seq){
            
            for (psi1 in psi1_seq){
              
              julia_assign("state", state)
              param=Get_PA_parameters()
              param["delta"]=disp
              param["cintra"]=.3
              param["alpha_0"]=cinter
              param["S"] = S
              param["psi_1"]=psi1
              param["psi_2"]=psi2
              param["f"]=facil
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
          write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                                cinter,"_branch_",branch,
                                "_stress_",S,"_delta_",disp,"_facilitation_",facil,".csv"),sep=";")
          
        } # end loop interspecific competition
        
      } # end loop branch
      
    } # end loop first species trait
    
  } #end loop dispersal
  
}#end facilitation loop



c_inter_seq=c(0,.1,.2,.3,.4)
stress_seq=seq(0,1,length.out=30)


d=tibble()  
d_bistab=tibble()
for (stress in stress_seq){
  
  for (cinter in c_inter_seq){
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                               cinter,"_branch_Degradation_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                               cinter,"_branch_Restoration_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"))
    d=rbind(d,d2)
    
    
    
  } # end loop interspecific competition
} # end loop branch

d2=d
d2[,1:2][d2[,1:2] < 10^-4] = 0
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

d_state$multistab=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c( "Species 2/Stress_tolerant","Species 2/Desert","Coexistence/Desert",
                               "Stress_tolerant/Species 2","Stress_tolerant/Desert")){
    return(1)
  } else {return(0)}
  
})

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in c(.3,.4)){#unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert"))))
    }
  }
}

multistability=d_state%>%
  group_by(., Psi1,Psi2,alpha_0)%>%
  summarise(., .groups = "keep",frac_multi=sum(multistab)/length(unique(d_state$Stress)))

ggplot(multistability)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=frac_multi))+
  facet_wrap(.~alpha_0)+
  the_theme+
  scale_fill_viridis_c()


ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi))+
  facet_wrap(.~alpha_0)+
  the_theme+
  scale_fill_viridis_c()


color_plot=c("Coexistence" = "#D8CC7B",
             "Species 2" = "#ACD87B",
             "Coexistence/Species 2" = "#DDEFCA",
             "Species 2/Coexistence" = "#DDEFCA",
             "Stress_tolerant" = "#7BD8D3",
             "Stress_tolerant/Desert" ="#0F8E87",
             "Desert/Stress_tolerant" ="#0F8E87",
             "Coexistence/Stress_tolerant"="#9BBBB9",
             "Stress_tolerant/Coexistence"="#9BBBB9",
             "Coexistence/Desert"="#C19E5E",
             "Desert/Coexistence"="#C19E5E",
             "Desert"=  "#696969",
             "Stress_tolerant/Species 2" = "#C998CE",
             "Species 2/Stress_tolerant" = "#C998CE")

pdf("./All_plot.pdf",width = 7,height = 5)
for (i in unique(d_state$Psi2)){
  print(d_state%>%
    filter(., Psi2==i)%>%
    ggplot(.)+
    geom_tile(aes(x=Stress,y=Psi1,fill=all_state))+
    the_theme+
    facet_wrap(.~alpha_0,labeller = label_bquote(cols = alpha[e]==.(alpha_0)))+
    scale_fill_manual(values=color_plot)+
    theme(legend.position = "none")+
    ggtitle(paste(round(i,3))))
  
}
dev.off()











#testing hysteresis size


c_inter_seq=c(0,.1,.2,.3,.4)
stress_seq=seq(0,1,length.out=30)


d=tibble()  
d_bistab=tibble()
for (stress in stress_seq){
  
  for (cinter in c_inter_seq){
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                               cinter,"_branch_Degradation_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                               cinter,"_branch_Restoration_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"))
    d=rbind(d,d2)
    
    
    
  } # end loop interspecific competition
} # end loop branch
d=d[order(d$Branches),]
d[,1:2][d[,1:2] < 10^-2] = 0


d_h=tibble()
for (a0 in unique(d$alpha_0)){
  for (p1 in unique(d$Psi1)){
    for (p2 in unique(d$Psi2)){
      
      if (p1>p2){ #as there is a symmetry
        d_fil=filter(d,Psi1==p1,Psi2==p2,alpha_0==a0)
        
        d_h=rbind(d_h,tibble(Psi1=p1,Psi2=p2,alpha_0=a0,
                             hyst=d_fil$Stress[min(which(d_fil$Competitive[1:26]==0))-1]-
                               d_fil$Stress[min(which(d_fil$Competitive[27:52]==0))]))
        
      }
    }
  }
}

ggplot(d_h)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=hyst))+
  facet_wrap(.~alpha_0)+
  the_theme+
  scale_fill_viridis_c()





# Step 3) Nspecies analysis ----
rm(list = ls())
source("./Dryland_shift_functions.R")


## 1) first exploration with CSI & rho_+ ----
d=tibble()
Nsp=15
list_csv=list.files('../Table/N_species/MF/Big_sim_Nsp/Random_ini/')

for ( k in 1:length(list_csv)){
  d2=read.table(paste0("../Table/N_species/MF/Big_sim_Nsp/Random_ini/",list_csv[k]),sep=",")
  colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
  d=rbind(d,d2)
}
d$Rho_p=rowSums(d[,1:Nsp])

d[d<10^(-4)]=0

d$CSI = sapply(1:nrow(d),function(x){
  set.seed(432)
  u=runif(Nsp)
  return(sum(d[x,1:Nsp]*u))})

d$Psi_normalized = sapply(1:nrow(d),function(x){
  trait=rev(seq(0,1,length.out=Nsp))
  if (d$CSI[x]==0){return(NA)
  }else{  return(sum(d[x,1:Nsp]*trait)/rowSums(d[x,1:Nsp]))
  }
})

d$Branch=sapply(1:nrow(d),function(x){
  return(ifelse(d$Branch[x]==1,"Degradation","Restoration"))
})

write.table(d,'../Table/N_species/MF/Random_ini.csv',sep=";")


d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")
Nsp=15

pD=ggplot(d%>%filter(., Branch=="Degradation"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+ggtitle("Degradation")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")

pR=ggplot(d%>%filter(., Branch=="Restoration"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+ggtitle("Restoration")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")
ggsave("../Figures/N_species/MF/CSI_random_ini.pdf",ggarrange(pD,pR,nrow = 2),width = 7,height = 10)

pD=ggplot(d%>%filter(., Branch=="Degradation"))+
  geom_line(aes(x=Stress,y=Rho_p,color=Psi_normalized,group=Random_ini),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Species density",color=expression(paste(bar(psi),"    ")))+ggtitle("Degradation")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")

pR=ggplot(d%>%filter(., Branch=="Restoration"))+
  geom_line(aes(x=Stress,y=Rho_p,color=Psi_normalized,group=Random_ini),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Species density",color=expression(paste(bar(psi),"    ")))+ggtitle("Restoration")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")
ggsave("../Figures/N_species/MF/Global_cover_random_ini.pdf",ggarrange(pD,pR,nrow = 2),width = 7,height = 10)




## 2) Species level : occurrence, size of tipping points  ----

treshold_sp=0.1
d_ASS_sp=tibble()

for (i in c(5,10,15,20,25,30,35)){
  
  d_t=tibble()
  
  Nsp=i
  traits=rev(seq(0,1,length.out=Nsp))
  list_csv=list.files(paste0('../Table/N_species/MF/',i,'_sp/'))
  
  for ( k in 1:length(list_csv)){ #for each replicate
    
    d2=read.table(paste0("../Table/N_species/MF/",i,"_sp/",list_csv[k]),sep=",")
    colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
    
    d2[d2<10^(-4)]=0
    d2=melt(d2,measure.vars=paste0("Sp_",1:Nsp))
    
    for (sp in 1:Nsp){
      
      d_sp=filter(d2,variable==paste0("Sp_",sp)) 
      
      
      if (any(d_sp$value[1:(nrow(d_sp)/2)]>0)){# we only keep species that are present in the system
        
        tipping=ifelse(any(abs(diff(d_sp$value[1:(nrow(d_sp)/2)]))>treshold_sp),1,0) #Is there a tipping point in the system for species sp
        
        d_ASS_sp=rbind(d_ASS_sp,tibble(Tipping=tipping,Species=sp,
                                       Trait=traits[sp],Branch="Degradation",
                                       Random_ini=k,Richness=Nsp))
        
      }
      
      if (any(d_sp$value[((nrow(d_sp)/2)+1):nrow(d_sp)] >0)){# we only keep species that are present in the system
        
        tipping=ifelse(any(abs(diff(d_sp$value[((nrow(d_sp)/2)+1):nrow(d_sp)]))>treshold_sp),1,0) #Is there a tipping point in the system for species sp
        
        d_ASS_sp=rbind(d_ASS_sp,tibble(Tipping=tipping,Species=sp,
                                       Trait=traits[sp],Branch="Restoration",
                                       Random_ini=k,Richness=Nsp))
        
      }
    }
    
  }
  
}

write.table(d_ASS_sp,"../Table/N_species/MF/d_ASS_sp.csv",sep=";")

d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=150


#Species-specific tipping points size & frequency
p=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .1)+
  facet_grid(Branch~Richness)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Fraction of tipping points")+
  scale_fill_gradientn(colours = color_Nsp(100))

ggsave("../Figures/N_species/MF/Fraction_tipping_points.pdf",p,width = 9,height = 5)





## 3) varying species diversity ----

threshold=.1
d_tipping=d_richness=d_richness2=d_tot=tibble()

for (i in c(5,10,15,20,25,30,35)){
  
  d_t=tibble()
  
  # pdf(paste0("../Figures/N_species/MF/Bifu_Nsp/Dyn_Nspecies_",i,".pdf"),width = 6,height = 4)
  
  Nsp=i
  list_csv=list.files(paste0('../Table/N_species/MF/',i,'_sp/'))
  
  for ( k in 1:length(list_csv)){
    
    d2=read.table(paste0("../Table/N_species/MF/",i,"_sp/",list_csv[k]),sep=",")
    colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
    
    d2[d2<10^(-4)]=0
    
    d_t=rbind(d_t,d2)
    
    
    d_richness=rbind(d_richness,
                     tibble(Nsp=Nsp,Random_ini=k,
                            Branch="Degradation",
                            Richness=length(which(colSums(d2[1:(nrow(d2)/2),])[1:Nsp] !=0))))
    
    d_richness=rbind(d_richness,
                     tibble(Nsp=Nsp,Random_ini=k,
                            Branch="Restoration",
                            Richness=length(which(colSums(d2[((nrow(d2)/2)+1):(nrow(d2)),])[1:Nsp] !=0))))
    
    
    # print(
    #   ggplot(d2%>%melt(., measure.vars=paste0("Sp_",1:Nsp)))+
    #     geom_line(aes(x=Stress,y=value,color=variable),size=.8)+
    #     facet_wrap(.~Branch)+
    #     the_theme+labs(y="",color=expression(paste(bar(psi),"    ")))+ggtitle("Degradation")
    # )
    
    
    d2$CSI = sapply(1:nrow(d2),function(x){
      set.seed(432)
      u=runif(Nsp)
      return(sum(d2[x,1:Nsp]*u))
    })
    
    d2$Psi_normalized = sapply(1:nrow(d2),function(x){
      trait=rev(seq(0,1,length.out=Nsp))
      if (d2$CSI[x]==0){return(NA)
      }else{  return(sum(d2[x,1:Nsp]*trait)/rowSums(d2[x,1:Nsp]))
      }
    })
    
    
    
    d_tot=rbind(d_tot,d2[,-c(1:Nsp)]%>%add_column(., Nsp=Nsp))
    
    for (sp in 1:Nsp){
      
      
      
      if (any(abs(diff(d2[1:(nrow(d2)/2),sp]))>threshold)){ #if species shift
        nb_shift=length(which(round((abs(diff(d2[1:(nrow(d2)/2),sp]))),4)>threshold))
        
        if (nb_shift %% 2==1){ #case for the first species shifting as it is already present in the system
          nb_shift=nb_shift+1
        }
        
        
        d_tipping=rbind(d_tipping,tibble(Species=sp,Nsp=Nsp,Random_ini=k,Tipping=nb_shift/2,Branch="Degradation"))
      }
      if (any(abs(diff(d2[((nrow(d2)/2)+1):(nrow(d2)),sp]))>threshold)){
        nb_shift=length(which(round((abs(diff(d2[((nrow(d2)/2)+1):(nrow(d2)),sp]))),4)>threshold))
        
        if (nb_shift %% 2==1){ #case for the first species shifting as it is already present in the system
          nb_shift=nb_shift+1
        }
        
        d_tipping=rbind(d_tipping,tibble(Species=sp,Nsp=Nsp,Random_ini=k,Tipping=nb_shift/2,Branch="Restoration"))
      }
      
      
    }
  }
  # dev.off()
  
  d_t=d_t[order(d_t$Branch),]
  d_richness2=rbind(d_richness2,tibble(Nsp=Nsp,
                                       Richness_D=length(which(colSums(d_t[1:(nrow(d_t)/2),])[1:Nsp] !=0)),
                                       Richness_R=length(which(colSums(d_t[((nrow(d_t)/2)+1):nrow(d_t),])[1:Nsp] !=0))))
}

d_tipping=d_tipping%>%
  group_by(., Nsp,Random_ini,Branch)%>%
  summarise(., Nb_different_species=length(unique(Species)),.groups = "keep")

write.table(d_tipping,"../Table/N_species/MF/Multistability_tipping.csv",sep=";")
write.table(d_richness,"../Table/N_species/MF/Multistability_richness.csv",sep=";")
write.table(d_richness2,"../Table/N_species/MF/Multistability_richness2.csv",sep=";")
write.table(d_tot,"../Table/N_species/MF/Multistability_CSI.csv",sep=";")


p1=ggplot(d_tipping%>%filter(., Branch=="Degradation"))+
  geom_bar(aes(x=Nsp,fill=as.factor(Nb_different_species),y=Nb_different_species),
           position = "fill",stat="identity",alpha=.75)+
  labs(x="Number of species",y="Fraction of random initial conditions",fill="Number of shifts")+
  scale_x_continuous(breaks = seq(5,35,by=5))+
  the_theme

p2=ggplot(d_richness%>%filter(., Branch=="Degradation"))+
  geom_bar(aes(x=Nsp,fill=as.factor(Richness),y=Richness),
           position = "fill",stat="identity",alpha=.75)+
  geom_point(data=d_richness2,aes(x=Nsp,y=Richness_D/Nsp),shape=0,size=2)+
  geom_line(data=d_richness2,aes(x=Nsp,y=Richness_D/Nsp))+
  scale_y_continuous(sec.axis=sec_axis(~., name="Fraction of species observed \n across replicates"))+
  labs(x="Number of species",y="Fraction of random initial conditions",fill="Number of species")+
  scale_x_continuous(breaks = seq(5,35,by=5))+
  the_theme


ggsave("../Figures/N_species/MF/Nb_shift_diversity.pdf",ggarrange(p1,p2,ncol = 2,labels = letters[1:2],align = "hv"),width=10,height=5)
