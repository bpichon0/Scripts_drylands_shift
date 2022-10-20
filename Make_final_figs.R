rm(list = ls())
source("./2_Species_analysis_functions.R")
source("./N_Species_analysis_functions.R")
dir.create("../Figures/Clean_figs/",showWarnings = F)



# Fig 1 : State diagram, stability and RNE ----

d_state=read.table("../Table/2_species/2_species_PA_global_comp.csv",sep=";")


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
a_values_bifu =  unique(d_state$alpha_0)[c(1, round(length( unique(d_state$alpha_0)) / 2), length( unique(d_state$alpha_0)))]


color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress-tolerant" = "#7BD8D3")

# state at equilibrium
p1 = ggplot(d_state) +
  geom_tile(aes(x = S, y = alpha_0 , fill = state)) +
  the_theme +
  annotate("text", x = rep(1.05, 3), y = a_values_bifu+max(a_values_bifu)/35 , label = c("e, f", "c, d", "a, b"), color = "black",size=5) +
  scale_fill_manual(values = color_rho[-2]) +
  theme(legend.position = "bottom") +
  labs(x = "Stress (S)", y = TeX(r'(Competition strength \ $\alpha_0)'), fill = "") +
  geom_hline(yintercept = a_values_bifu , lwd = .1, color = "gray40") +
  theme(legend.text = element_text(size = 12))



#bifurcation diagrams
S_critic1=read.table("../Table/2_species/Species_along_critical_stress_PA.csv",sep=";")[1,2]
S_critic2=read.table("../Table/2_species/Species_along_critical_stress_PA.csv",sep=";")[2,2]

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

#RNE
drne=read.table(paste0("../Table/2_species/2_species_RNE_PA.csv"),sep=";")

d_compe=filter(drne,Sp=="Competitive")
d_stesstol=filter(drne,Sp=="Stress-tolerant")
d_both=filter(drne,Sp=="Both")

d_RNE=tibble(S=d_both$S,alpha_0=d_both$alpha_0, 
             RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1),
             RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2))
d_RNE[,3:4][is.na(d_RNE[,3:4])] = NA
# d_RNE$RNE_stress_tol[246] = NA #to avoidnumeric problems



for (a0 in 1:3) {
  assign(
    paste0("p_", a0),
    ggplot(d_RNE %>% filter(., alpha_0 == c(0,.15,.3)[a0]) %>%melt(., id.vars=c("S","alpha_0")) %>%
             mutate(., variable = recode_factor(variable, "RNE_stress_tol" = "Stress-tolerant", "RNE_competitive" = "Competitive"))) +
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
      labs(x="Stress (S)",y="RNE",color="")+ylim(-1,1)+
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

ggsave("../Figures/Clean_figs/Fig_1.pdf",width = 12,height = 6)

