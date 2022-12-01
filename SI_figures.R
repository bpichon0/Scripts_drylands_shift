source('./Dryland_shift_functions.R')
library(grid)


# 2-species ----
## Multistability along competition gradient with local competition ----
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


scale="local"
d2t=transform(d_state%>% 
                filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
                mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6)),
              Scena=factor(Scena,levels=c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")),
                           labels=c("Global facilitation","Local facilitation")))

color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")


d2t$Multistability=sapply(1:nrow(d2t),function(x){
  if (d2t$all_state[x] %in% c("Competitive","Coexistence","Desert","Stress_tolerant")){
    return("no_multi")
  } else{
    return("multi")
  }
})


p1_pattern=ggplot(d2t%>%
                    mutate(.,all_state=recode_factor(all_state,
                                                     "Competitive/Coexistence"="Coexistence/Competitive",
                                                     "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
                                                     "Stress_tolerant"="Stress-tolerant",
                                                     "Stress_tolerant/Desert"="Stress-tolerant/Desert"))) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',
                    pattern_spacing = 0.05, 
                    pattern_density = 0.01, 
                    pattern_fill    = alpha("black",.7),
                    pattern_colour  = alpha("black",.7))+
  
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Coexistence/Competitive" = "#DDEFCA",
                             "Stress-tolerant" = "#7BD8D3",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(Scena~Delta,labeller=label_bquote(cols = delta == .(Delta)))+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9))+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))


Fig_2_SI=p1_pattern+guides(pattern=F)

ggsave("../Figures/Final_figs/SI/Fig_2_with_local_competition.pdf",Fig_2_SI,width = 8,height = 6)


















## Multistability but with mean-field model ----

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

d_state$Multistability=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c("Competitive","Coexistence","Desert","Stress_tolerant")){
    return("no_multi")
  } else{
    return("multi")
  }
})

p=ggplot(d_state) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',
                    pattern_spacing = 0.05, 
                    pattern_density = 0.01, 
                    pattern_fill    = alpha("black",.7),
                    pattern_colour  = alpha("black",.7))+
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
                             "Desert"=  "#696969",
                             "Competitive/Stress_tolerant"="#C998CE"
  ),labels=c("Stress_tolerant"="Stress-tolerant","Stress_tolerant/Desert"="Stress-tolerant/Desert",
             "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
             "Competitive/Stress_tolerant"="Competitive/Stress-tolerant"))+
  geom_hline(yintercept = round(unique(d2$alpha_0)[c(1,100)],4))+
  geom_text(data=tibble(x=c(1.025,1.025),y=c(0.01,.31),lab=c("c",'b')),aes(x=x,y=y,label=lab))+
  the_theme+theme(legend.text = element_text(size=10))+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))


for (u in 1:2){
  d_plot=filter(d2,round(alpha_0,4) == round(unique(d2$alpha_0)[c(1,100)[u]],4))
  assign(paste0("p_1_",u),
         ggplot(d_plot%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
           geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
           labs(x = "Stress (S)", y = "Cover", color = "",linetype="") +
           the_theme +
           scale_color_manual(values = color_rho[c(2, 4)]) +
           scale_linetype_manual(values=c(1,9))+
           theme(legend.text = element_text(size = 11))+
           scale_x_continuous(breaks = c(0,.5,1))+guides(color=F)
  )
  
  assign(paste0("p_2_",u),
         ggplot(d_plot%>%melt(., measure.vars=c("Rho_plus")))+
           geom_line(aes(x = Stress, y = value, linetype=Branches),lwd=.8) +
           labs(x = "Stress (S)", y = "Cover", color = "",linetype="") +
           the_theme +
           scale_color_manual(values = color_rho[c(2, 4)]) +
           scale_linetype_manual(values=c(1,9))+
           theme(legend.text = element_text(size = 11))+
           scale_x_continuous(breaks = c(0,.5,1))+guides(color=F)
  )
  
  
}

p2=ggarrange(p_1_2+xlab(""),p_2_2+xlab(""),p_1_1+xlab(""),p_2_1,nrow = 4,common.legend = T,legend = "bottom",labels = c(letters[2],"",letters[3],""))


ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_MF.pdf",ggarrange(p+guides(pattern=F),p2,ncol = 2,widths = c(2,1),labels=c(letters[1],"")),
       width = 10,height = 6)




















## Multistability varying traits with mean-field model ----

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
    
    d_state$Multistability=sapply(1:nrow(d_state),function(x){
      if (d_state$all_state[x] %in% c("Species 2","Coexistence","Desert","Stress_tolerant")){
        return("no_multi")
      } else{
        return("multi")
      }
    })
    
    color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
    
    appender <- function(string) {
      TeX(paste("$\\alpha_e = $", string))}
    
    
    
    
    p1=ggplot(d_state%>%
                mutate(all_state=recode_factor(all_state,
                                               "Desert/Species 2"="Competitive/Desert",
                                               "Desert/Coexistence"="Coexistence/Desert",
                                               "Species 2"="Competitive",
                                               "Desert/Stress_tolerant"="Stress-tolerant/Desert",
                                               "Coexistence/Species 2"="Coexistence/Competitive",
                                               "Stress_tolerant"= "Stress-tolerant",
                                               "Stress_tolerant/Coexistence"="Coexistence/Stress-tolerant",
                                               "Stress_tolerant/Species 2"="Competitive/Stress-tolerant",
                                               "Species 2/Coexistence"="Coexistence/Competitive"))) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),              
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = alpha("black",.7),
                        pattern_colour  = alpha("black",.7))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = "Stress (S)", y = TeX(r'(Trait competitive sp., \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Competitive" = "#ACD87B",
                                 "Coexistence/Competitive" = "#DDEFCA",
                                 "Stress-tolerant" = "#7BD8D3",
                                 "Stress-tolerant/Desert" ="#0F8E87",
                                 "Coexistence/Stress-tolerant"="#9BBBB9",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Desert"=  "#696969",
                                 "Competitive/Stress-tolerant" = "#C998CE"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
      ggtitle(TeX("$\\psi_1 = 1$"))+
      scale_pattern_manual(values=rev(c("none" ,"stripe")))
    
    
    
    
    
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
    
    d_state$Multistability=sapply(1:nrow(d_state),function(x){
      if (d_state$all_state[x] %in% c("Species 2","Coexistence","Desert","Competitive")){
        return("no_multi")
      } else{
        return("multi")
      }
    })
    
    
    p2=ggplot(d_state%>%
                mutate(all_state=recode_factor(all_state,
                                               "Species 2/Coexistence"="Coexistence/Stress-tolerant",
                                               "Desert/Species 2"="Stress-tolerant/Desert",
                                               "Desert/Coexistence"="Coexistence/Desert",
                                               "Desert/Competitive"="Competitive/Desert",
                                               "Competitive" = "Competitive",
                                               "Species 2" = "Stress-tolerant",
                                               "Coexistence/Competitive" = "Coexistence/Competitive"
                ),
                Psi2=round(Psi2,5),
                Stress=round(Stress,5))) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),              
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = alpha("black",.7),
                        pattern_colour  = alpha("black",.7))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = "Stress (S)", y = TeX(r'(Trait stress-tolerant sp., \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Stress-tolerant" = "#7BD8D3",
                                 "Competitive" = "#ACD87B",
                                 "Coexistence/Stress-tolerant" = "#9BBBB9",
                                 "Stress-tolerant/Desert" ="#0F8E87",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Coexistence/Competitive" = "#DDEFCA",
                                 "Coexistence/Stress-tolerant"="#9BBBB9",
                                 "Desert"=  "#696969"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      scale_y_continuous(sec.axis = sec_axis(trans = ~ (.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
      ggtitle(TeX("$\\psi_2 = 0$"))+
      scale_pattern_manual(values=rev(c("none" ,"stripe")))
    
    
    
  }
}

ggsave("../Figures/Final_figs/SI/Multistability_varying_traits_MF.pdf",
       ggarrange(p1+guides(pattern=F)+theme(axis.line.x = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank()),
                 p2+guides(pattern=F)+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
                 nrow = 2,heights=c(1,1),labels=letters[1:2],common.legend = T,legend = "bottom"),
       width = 10,height = 7)



## Tipping points MF model varying trait ----


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
    
    p1=ggplot(d_tipping)+
      geom_smooth(aes(x=Psi2,y=Tipping_C,color=Branch,group=interaction(Competition,Branch)),se = F)+
      the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
      labs(x=TeX(r'(Trait competitive sp., \ $\psi_2)'),y="Threshold stress",color="")+
      scale_color_manual(values=c("black","blue"))+
      scale_x_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)'),
                                             breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1")),
                         breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1"))
    
  } 
  if (Psi_sp1==0){
    
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
    
    p2=ggplot(d_tipping)+
      geom_smooth(aes(x=Psi2,y=Tipping_ST,color=Branch,group=interaction(Competition,Branch)),se = F)+
      the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
      labs(x=TeX(r'(Trait stress-tolerant sp., \ $\psi_1)'),y="Threshold stress",color="")+
      scale_color_manual(values=c("black","blue"))+
      scale_x_continuous(sec.axis = sec_axis(trans = ~ (.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)'),
                                             breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1")),
                         breaks = c(0,.25,.5,.75,1),labels = c("0","0.25","0.5","0.75","1"))
    
    
  }
}

p_tot=ggarrange(p1+theme(legend.position = "none")+ggtitle(TeX("$\\psi_1 = 1$")),
                p2+ggtitle(TeX("$\\psi_2 = 0$"))+theme(strip.background.x = element_blank(),
                                                       strip.text.x = element_blank()),nrow=2,labels = letters[1:2])

ggsave("../Figures/Final_figs/SI/Tipping_points_traits_MF.pdf",width = 7,height = 7)


## Net-effects, balance competition-facilitation ----

epsilon=10^(-6)
d2=read.table("../Table/2_species/PA/Net_effect_varying_traits.csv",sep=";")

net_effect =sapply(seq(1, nrow(d2) , by = 2),function(x){
  return((d2$Eq[x+1] - d2$Eq[x]) / epsilon)
})


d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect


p11=ggplot(d_net%>%
             filter(.,Psi1==1,alpha_0%in% c(.1),Species==1)%>%
             mutate(.,Species=recode_factor(Species,"1" ="Species 1","2" = "Species 2")))+
  ggtitle(TeX("$\\psi_1=1, \\alpha_e=0.1$"))+
  geom_line(aes(x=S,y=value,color=Psi2,group=interaction(Psi1,Psi2,alpha_0,Species)))+
  labs(x="Stress, S",y=TeX(r'($\frac{\partial \rho_{+_2}}{\partial \psi_2})'),
       color=TeX(r'($\ \psi_2 \ \ )'))+
  the_theme+
  scale_color_stepsn(colours=color_Nsp(6),breaks=seq(0,1,by=.2))+
  geom_hline(yintercept = 0,linetype=9)

p12=ggplot(d_net%>%
             filter(.,Psi1==1,alpha_0%in% c(.1),Species==2)%>% #see d2, species are inverted for \psi_1=1
             mutate(.,Species=recode_factor(Species,"1" ="Species 1","2" = "Species 2")))+
  ggtitle(TeX("$\\psi_1=1, \\alpha_e=0.1$"))+
  geom_line(aes(x=S,y=value,color=Psi2,group=interaction(Psi1,Psi2,alpha_0,Species)))+
  labs(x="Stress, S",y=TeX(r'($\frac{\partial \rho_{+_1}}{\partial \psi_2})'),
       color=TeX(r'($\ \psi_2 \ \ )'))+
  the_theme+
  scale_color_stepsn(colours=color_Nsp(6),breaks=seq(0,1,by=.2))+
  geom_hline(yintercept = 0,linetype=9)

p21=ggplot(d_net%>%
             filter(.,Psi1==0,Species==2,alpha_0 %in% c(.1))%>%
             mutate(.,Species=recode_factor(Species,"1" ="Species 1","2" = "Species 2")))+ #here we inverse to be coherent with the analysis made 
  ggtitle(TeX("$\\psi_2=0, \\alpha_e=0.1$"))+
  geom_line(aes(x=S,y=value,color=Psi2,group=interaction(Psi1,Psi2,alpha_0,Species)))+
  labs(x="Stress, S",y=TeX(r'($\frac{\partial \rho_{+_2}}{\partial \psi_1})'),
       color=TeX(r'($\ \psi_1 \ \ )'))+
  the_theme+
  scale_color_stepsn(colours=color_Nsp(6),breaks=seq(0,1,by=.2))+
  geom_hline(yintercept = 0,linetype=9)


p22=ggplot(d_net%>%
             filter(.,Psi1==0,Species==2,alpha_0 %in% c(.2))%>%
             mutate(.,Species=recode_factor(Species,"1" ="Species 1","2" = "Species 2")))+ #here we inverse to be coherent with the analysis made 
  ggtitle(TeX("$\\psi_2=0, \\alpha_e=0.2$"))+
  geom_line(aes(x=S,y=value,color=Psi2,group=interaction(Psi1,Psi2,alpha_0,Species)))+
  labs(x="Stress, S",y=TeX(r'($\frac{\partial \rho_{+_2}}{\partial \psi_1})'),
       color=TeX(r'($\ \psi_1 \ \ )'))+
  the_theme+
  scale_color_stepsn(colours=color_Nsp(6),breaks=seq(0,1,by=.2))+
  geom_hline(yintercept = 0,linetype=9)


p_tot=ggarrange(ggarrange(p11+theme(),
                          p12,ncol=2,labels = c(letters[1:2]),common.legend = T,legend = "bottom"),
                ggarrange(p21,
                          p22+ylab(""),ncol = 2,labels = c(letters[3:4]),common.legend = T,legend = "bottom"),
                nrow=2)

ggsave("../Figures/Final_figs/SI/Net_effects_varying_traits.pdf",p_tot,width = 7,height = 7)




## Clustering of species 1 along dispersal gradient & species/fertile sites pairs ----
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
           filter(., alpha_0==.2, S %in% c(.73, .77))%>%
           group_by(., cintra,alpha_0,S,delta,Scena)%>%
           summarise(.,.groups ="keep",c11=mean(c11) ))+
  geom_point(aes(x=as.numeric(delta),y=c11),size=1.5,alpha=.7,shape=1)+
  facet_grid(Scena~S,labeller=label_bquote(cols = Stress == .(S)),scales = "free")+
  the_theme+labs(x=TeX("$\\delta$"),y=expression(paste("Stress-tolerant clustering (c"[11],")")))+
  scale_color_manual(values=as.character(color_rho[c(2,4)]))+
  geom_hline(yintercept = 1)+ylim(0,3)
ggsave("../Figures/Final_figs/SI/Clustering_11_dispersal.pdf",p,width = 7,height = 4)

p=ggplot(d_clustering%>%
         filter(., alpha_0==.2, S %in% c(0))%>%
         melt(., measure.vars=c("Rho_20","Rho_10"))%>%
         mutate(.,variable=recode_factor(variable,"Rho_10"="Pair Stress-tolerant sp./Fertile",
                                "Rho_20"="Pair Competitive sp./Fertile")))+
  geom_point(aes(x=as.numeric(delta),y=value,color=variable),size=1.5)+
  ggtitle("Stress (S) = 0")+
  the_theme+labs(x=TeX("$\\delta$"),y=expression(paste("Cover")))+
  facet_grid(Scena~.,scales = "free")+
  
  scale_color_manual(values=c("black","grey"))+
  labs(color="")

ggsave("../Figures/Final_figs/SI/Pair_10_20.pdf",p,width = 5,height = 4)





## NintA for MF model and PA model ----

load(file="../Table/2_species/MF/Comparing_net_effects.RData")
d_net=d_all[[1]];d_RNE=d_all[[2]]


d_net$value[d_net$value>10]=NA
d_net$value[d_net$value < -10]=NA

p1=ggplot(d_RNE%>%
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
                          variable == "Competitive"), aes(y = .32,x=1.05),label="b")+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .3))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .12,x=1.05),label="c")+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")

alpha_for_plot=c(.1, .3)

for (alpha in 1:2){
  assign(paste0("p1_",alpha),ggplot(d_RNE%>%
                                      filter(., alpha_0==alpha_for_plot[alpha])%>%
                                      melt(., measure.vars=c("NintA_comp")))+#%>%
           #           mutate(., value=rescale(value,to=c(-1,1))))+
           geom_hline(yintercept = 0,linetype=9,lwd=.5)+
           geom_line(aes(x=S,y=value))+
           the_theme+theme(axis.text = element_text(size =9), axis.title = element_text(size = 10))+
           labs(x="Stress (S)",y=TeX("$\\NInt_{A}$")))
}

p_right=ggarrange(p1_2+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                             axis.line.x = element_blank(),
                             axis.ticks.x = element_blank())+
                    labs(x=""),p1_1,ggplot() + theme_void(),nrow=3,labels=letters[2:3],heights = c(1,1,.1))

p_tot_1=ggarrange(p1+theme(legend.position = "none",
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),axis.line.x = element_blank(),
                           axis.ticks.x = element_blank()),p_right,widths = c(2.2,1),labels = c(letters[1],""))



load(file="../Table/2_species/PA/Comparing_net_effects.RData")
d_RNE=d_all



d_RNE[,3:6][d_RNE[,3:6] > 2]=NA
d_RNE[,3:6][d_RNE[,3:6] < -1]=NA

p2=ggplot(d_RNE%>%
            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
            filter(.,Disp==.1,Scale=="LocalF_GlobalC")%>%
            #mutate(., value=rescale(value,to=c(-1,1)))%>%
            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_grid(.~factor(variable,levels=c("Stress-tolerant","Competitive")))+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  geom_hline(data=tibble(variable="Competitive",yint=c(0.1,.3)),aes(yintercept=yint,variable="Competitive"))+
  geom_text(data=tibble(variable="Competitive",y=c(0.13,.33),x=c(1.04,1.04),label=c("f","e")),aes(y=y,x=x,label=label,variable="Competitive"))


alpha_for_plot=c(unique(d_RNE$alpha_0)[(50/3)+1], .3)

for (alpha in 1:2){
  assign(paste0("p2_",alpha),ggplot(d_RNE%>%
                                      filter(., alpha_0==alpha_for_plot[alpha],Scale=="LocalF_GlobalC",Disp==.1)%>%
                                      melt(., measure.vars=c("NintA_comp")))+#%>%
           #           mutate(., value=rescale(value,to=c(-1,1))))+
           geom_hline(yintercept = 0,linetype=9,lwd=.5)+
           geom_line(aes(x=S,y=value))+
           the_theme+theme(axis.text = element_text(size =9), axis.title = element_text(size = 10))+
           labs(x="Stress (S)",y=TeX("$\\NInt_{A}$")))
}

p_right=ggarrange(p2_2+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                             axis.line.x = element_blank(),
                             axis.ticks.x = element_blank())+
                    labs(x=""),p2_1,ggplot() + theme_void(),nrow=3,labels=letters[5:6])

p_tot_2=ggarrange(p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),p_right,widths = c(2.2,1),labels = c(letters[4],""))

p_tot=ggarrange(p_tot_1,p_tot_2,
                nrow=2,heights = c(1,1.4))
ggsave("../Figures/Final_figs/SI/NintA_MF_PA.pdf",p_tot,width = 7,height = 6)




## Niche expansion in MF model ----

d_niche = read.table("../Table/2_species/MF/Niche_expansion_varying_traits.csv",sep=";")
d_niche = do.call(data.frame,lapply(d_niche,function(x) replace(x, is.infinite(x), NA)))


d_niche=read.table("../Table/2_species/MF/Niche_expansion_MF.csv",sep=";")
alpha_seq=seq(0,.3,length.out=10)
d_niche$alpha_0=alpha_seq
p1=ggplot(d_niche%>%
            melt(.,measure.vars=c("Delta_niche_2")))+
  geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
  the_theme+
  scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
  labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="% of competitive species niche change")


#niche varying traits
d_niche = read.table("../Table/2_species/MF/Niche_expansion_varying_traits.csv",sep=";")
d_niche = do.call(data.frame,lapply(d_niche,function(x) replace(x, is.infinite(x), NA)))

d_niche$Delta_niche_2[1:(nrow(d_niche)/2)] = d_niche$Delta_niche_1[1:(nrow(d_niche)/2)] 



p3_1=ggplot(NULL)+
  geom_path(data=d_niche%>%filter(., Psi1==0)%>%
              melt(., measure.vars=c("Delta_niche_2"))%>%
              mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
              mutate(., Delta_psi=abs(Psi1-Psi2),
                     Facilitation=as.factor(Facilitation),
                     alpha_0=as.factor(alpha_0)),
            aes(x=Psi2,y=value,group=interaction(alpha_0,Facilitation),color=alpha_0),lwd=1)+
  geom_point(data=d_niche%>%filter(., Psi1==0)%>%
               melt(., measure.vars=c("Delta_niche_2"))%>%
               mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
               mutate(., Delta_psi=abs(Psi1-Psi2),
                      Facilitation=as.factor(Facilitation),
                      alpha_0=as.factor(alpha_0))%>%
               filter(., Psi2 %in% unique(.$Psi2)[seq(1,length(unique(.$Psi2)),by =5)]), #to only have 10 points
             aes(x=Psi2,y=value,shape=Facilitation,color=alpha_0),fill="white",size=2)+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  labs(y = "% of species 1 niche change", x = TeX(r'(Trait species 2, \ $\psi_2)'),
       color=TeX(r'(\ \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f)'))+ggtitle(TeX("$\\psi_1 = 0$"))+
  the_theme

p3_2=ggplot(NULL)+
  geom_path(data=d_niche%>%filter(., Psi1==1)%>%
              melt(., measure.vars=c("Delta_niche_2"))%>%
              mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
              mutate(., Delta_psi=abs(Psi1-Psi2),
                     Facilitation=as.factor(Facilitation),
                     alpha_0=as.factor(alpha_0)),
            aes(x=Psi2,y=value,group=interaction(alpha_0,Facilitation),color=alpha_0),lwd=1)+
  geom_point(data=d_niche%>%filter(., Psi1==1)%>%
               melt(., measure.vars=c("Delta_niche_2"))%>%
               mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
               mutate(., Delta_psi=abs(Psi1-Psi2),
                      Facilitation=as.factor(Facilitation),
                      alpha_0=as.factor(alpha_0))%>%
               filter(., Psi2 %in% unique(.$Psi2)[seq(1,length(unique(.$Psi2)),by =5)]), #to only have 10 points
             aes(x=Psi2,y=value,shape=Facilitation,color=alpha_0),fill="white",size=2)+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  labs(y = "% of species 2 niche change", x = TeX(r'(Trait species 2, \ $\psi_2)'),
       color=TeX(r'(\ \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f)'))+ggtitle(TeX("$\\psi_1 = 1$"))+
  the_theme
p3=ggarrange(p3_2,p3_1,labels = letters[2:3],common.legend = T,legend = "bottom")

p_tot = ggarrange(ggarrange(ggplot() + theme_void(),p1,ggplot() + theme_void(),widths = c(.7,2,.7),ncol=3,labels=c("",letters[1],"")),
                  p3,nrow=2,labels=c("",""))

ggsave("../Figures/Final_figs/SI/Niche_expansion_MF.pdf",p_tot,width = 8,height = 8)



## NIntA for varying traits ----

d_RII=read.table("../Table/2_species/PA/RII_varying_traits.csv",sep=";")

d_compe=filter(d_RII,Sp=="Competitive")
d_stesstol=filter(d_RII,Sp=="Stress_tolerant")
d_both=filter(d_RII,Sp=="Both")

d_RII=tibble(S=d_both$S,alpha_0=d_both$alpha_0, Psi1=d_both$Psi1,Psi2=d_both$Psi2,Facilitation=d_both$Facilitation,
             RNE_stress_tol = (d_both$rho_1-d_stesstol$rho_1)/(d_both$rho_1+d_stesstol$rho_1), #acutally its RII
             RNE_competitive = (d_both$rho_2-d_compe$rho_2)/(d_both$rho_2+d_compe$rho_2),
             NintA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
             NintA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)),
             NImpA_comp = 2*(d_both$rho_2-d_compe$rho_2)/(2*max(d_compe$rho_2)-d_compe$rho_2+abs(d_both$rho_2-d_compe$rho_2)), #using metrics from Diaz-Sierre MEE 2017
             NImpA_st = 2*(d_both$rho_1-d_stesstol$rho_1)/(2*max(d_stesstol$rho_1)-d_stesstol$rho_1+abs(d_both$rho_1-d_stesstol$rho_1)))
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

ggsave("../Figures/Final_figs/SI/NInt_A_varying_traits.pdf",p_tot,width = 7,height = 8)


## Advantage of CSI compared to Rho+----

c_inter_seq=c(.2)
psi1_seq=c(1)
disp=.1
f=0.9

d=tibble()  
for (Psi_sp1 in psi1_seq){
  for (branch in c("Restoration","Degradation")){
    for (cinter in c_inter_seq){
      
      d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                           cinter,"_branch_",branch,
                           "_Psi1_",Psi_sp1,"_delta_",disp,"_facilitation_",f,".csv"),sep=";")
      
      d=rbind(d,d2)
      
    } # end loop interspecific competition
  } # end loop branch
}



d[,1:2][d[,1:2] < 10^-4] = 0


#COmputing CSI index
set.seed(123)
u=runif(2)
d$CSI = sapply(1:nrow(d),function(x){
  return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
})


d2=filter(d,Psi1==1)

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


d_state$Multistability=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c("Species 2","Coexistence","Desert","Stress_tolerant")){
    return("no_multi")
  } else{
    return("multi")
  }
})

d_state=d_state%>%
  mutate(all_state=recode_factor(all_state,
                                 "Species 2/Coexistence"="Coexistence/Species 2",
                                 "Desert/Species 2"="Species 2/Desert",
                                 "Desert/Coexistence"="Coexistence/Desert",
                                 "Stress_tolerant"="Species 1",
                                 "Desert/Stress_tolerant"="Species 1/Desert",
                                 "Stress_tolerant/Coexistence"="Coexistence/Species 1",
                                 "Stress_tolerant/Species 2"="Species 2/Species 1",
                                 "Species 2/Coexistence"="Coexistence/Species 2"))

color_multistability=c("Coexistence" = "#D8CC7B",
                       "Species 2" = "#ACD87B",
                       "Coexistence/Species 2" = "#DDEFCA",
                       "Species 1" = "#7BD8D3",
                       "Species 1/Desert" ="#0F8E87",
                       "Coexistence/Species 1"="#9BBBB9",
                       "Coexistence/Desert"="#C19E5E",
                       "Desert"=  "#696969",
                       "Species 2/Species 1" = "#C998CE")

#Bifurcation diagram

p=ggplot(d2%>%melt(., measure.vars=c("Rho_plus","CSI"))%>%
           filter(., Psi2 %in% unique(d_state$Psi2)[c(80)])%>%
           mutate(Psi2=round(abs(Psi2),2))%>%
           filter(., alpha_0==.2))+
  geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
  labs(x = TeX(r'(Stress, \ $S)'), y = "", color = "",linetype="") +
  the_theme +
  scale_color_manual(values = c("black","red"),labels=c("Global cover","Community index")) +
  scale_linetype_manual(values=c(1,9))+
  theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))+
  new_scale_color() +
  geom_point(data=d_state%>%
               filter(., Psi2==unique(d_state$Psi2)[c(80)],alpha_0==.2)%>%
               mutate(., y=.83),
             aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=2)+
  scale_fill_manual(values=color_multistability)+guides(color=F,fill=F)+
  scale_color_manual(values=color_multistability)


ggsave("../Figures/Final_figs/SI/Using_CSI.pdf",p,width = 6,height = 4)

## Competition experienced ----
Nsp=15
trait=rev(seq(0,1,length.out=Nsp))
comp=comp2=c()
for (i in 1:Nsp){
  comp=c(comp,sum(sapply(1:Nsp,function(x){
    exp(-abs(trait[i]-trait[x]))
  })))
  comp2=c(comp2,sum(sapply(1:Nsp,function(x){
    trait[i]*exp(-abs(trait[i]-trait[x]))
  })))
}


p1=ggplot(tibble(Com=comp,Trait=trait,Comp2=comp2))+
  geom_line(aes(x=Trait,y=Com))+
  geom_point(aes(x=Trait,y=Com,color=Trait),size=3)+
  geom_text(data=tibble(x=c(0,trait[8],1),y=c(9.75,11.75,9.75),text=c("i","ii","iii")),aes(x=x,y=y,label=text),size=4)+
  the_theme+
  labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $\sum_j e^{-|\psi_i-\psi_j|})'),color="")+
  scale_color_gradientn(colours = color_Nsp(Nsp))+
  theme(legend.position = "none")

trait_subplot=c(0,trait[8],1)
for (i in 1:3){
  
  comp=sapply(1:Nsp,function(x){
    exp(-abs(trait_subplot[i]-trait[x]))
  })
  
  if (i != 3){
    assign(paste0("p_",i),
     ggplot(tibble(Com=comp,Trait=trait))+
       geom_line(aes(x=Trait,y=Com))+
       labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $e^{-|\psi_i-\psi_j|})'))+the_theme+
       theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),
             axis.ticks.x = element_blank())+
       scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))
    )
  }else {
    assign(paste0("p_",i),
           ggplot(tibble(Com=comp,Trait=trait))+
             geom_line(aes(x=Trait,y=Com))+
             labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $e^{-|\psi_i-\psi_j|})'))+the_theme+
             scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))
    )
  }
  
}



p2=ggplot(tibble(Com=comp,Trait=trait,Comp2=comp2))+
  geom_line(aes(x=Trait,y=Comp2))+
  geom_point(aes(x=Trait,y=Comp2,color=Trait),size=3)+
  the_theme+
  labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $\psi_i \sum_j e^{-|\psi_i-\psi_j|})'),color="")+
  scale_color_gradientn(colours = color_Nsp(15))

ggsave("../Figures/Final_figs/SI/Competition_experienced.pdf",
       ggarrange(
         ggarrange(p1,
                   ggarrange(p_1,p_2,p_3,nrow = 3,labels = c('i','ii','iii'),heights = c(1,1,1.4)),labels = c("a","","",""),widths = c(2,1)),
         p2,ncol=2,common.legend = T,legend = "none",labels = c("",letters[2]),widths = c(1.5,1)),
       width = 11,height = 4)



## Comparing CA & PA ----

#First layer = Simulations with PA (~Fig 2)

d2=read.table(paste0("../Table/2_species/PA/Multistability_fixed_traits_PA.csv"),sep=";")%>%
  filter(., Delta==.1)

d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Stress-tolerant")
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


scale="global"
d2t=transform(d_state%>% 
                filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
                mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6)),
              Scena=factor(Scena,levels=c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")),
                           labels=c("Global facilitation","Local facilitation")))

#Second layer : CA simulations

list_f=list.files("../Table/2_species/CA/Justifying_use_PA/")
d=tibble()
for (i in list_f){
  d2=read.table(paste0("../Table/2_species/CA/Justifying_use_PA/",i),sep=",")
  mean_density=c(colMeans(d2[-c(1:1000),-1]),
                 str_split(i,pattern = "_")[[1]][3],
                 str_split(i,pattern = "_")[[1]][5],
                 gsub(pattern = ".csv",replacement = "",x=str_split(i,pattern = "_")[[1]][7]))
  d=rbind(d,mean_density)  
}
colnames(d)=c("Rho_1","Rho_2","Rho_0","Rho_d","Stress","Branches","alpha_0")
d=as_tibble(d)%>%
  mutate(., Rho_1=as.numeric(Rho_1),Rho_2=as.numeric(Rho_2),Rho_0=as.numeric(Rho_0),Rho_d=as.numeric(Rho_d),
         alpha_0=as.numeric(alpha_0),Stress=as.numeric(Stress))

d$state = sapply(1:nrow(d), function(x) {
  if (d[x, 1] > 0 & d[x, 2] > 0) {
    return("Coexistence")
  }
  if (d[x, 1] > 0 & d[x, 2] == 0) {
    return("Stress-tolerant")
  }
  if (d[x, 1] == 0 & d[x, 2] > 0) {
    return("Competitive")
  }
  if (d[x, 1] == 0 & d[x, 2] == 0) {
    return("Desert")
  }
})

d=d[order(d$alpha_0,d$Stress),]

all_state =sapply(seq(1, nrow(d) , by = 2),function(x){
  if (d$state[x] != d$state[x+1]){
    return(paste0(d$state[x],"/", d$state[x+1]))
  }
  else {return(d$state[x])}
})


d_state=d%>%
  filter(., Branches=="Degradation")%>%
  select(.,-Branches)
d_state$all_state=all_state

p=ggplot(NULL) +
  geom_tile(data=d2t,aes(x=Stress,y=alpha_0,fill=all_state))+
  geom_point(data=d_state,aes(x=Stress,y=alpha_0,fill=all_state),size=5,shape=21)+
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Coexistence/Competitive" = "#DDEFCA",
                             "Stress-tolerant" = "#7BD8D3",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9))
ggsave("../Figures/Final_figs/SI/Comparizon_CA_PA.pdf",p,width = 7,height = 6)


# N-species ----
## Species specific tipping points ----


d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=150


#Species-specific tipping points size & frequency
p=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .1)+
  facet_grid(Branch~Richness)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Frequency of abrupt shift",fill=TeX("$\\psi$"))+
  scale_fill_gradientn(colours = color_Nsp(100))+
  scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))

ggsave("../Figures/Final_figs/SI/Specific_specific_tipping_points.pdf",p,width = 9,height = 5)


## CSI bifu diagrams----

d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")

Nsp=15
p1=ggplot(d%>%filter(.,Branch=="Degradation"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.8,shape=1)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f[0]==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(Nsp),na.value = "black")+
  guides(fill="none")+
  theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=12))

p2=ggplot(d%>%filter(.,Branch=="Restoration"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.8,shape=1)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f[0]==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(Nsp),na.value = "black")+
  guides(fill="none")+
  theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=12))


ggsave("../Figures/Final_figs/SI/CSI_competition_facilitation_15species.pdf",
       ggarrange(p1,
                 p2,
                 nrow = 2,labels = letters[1:2]),
       width = 7,height = 10)

## Varying the threshold for shifts ----

threshold_seq=c(.05,.075,.125,.15)

for (thresh in 1:length(threshold_seq)){
  
  threshold=threshold_seq[thresh]
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
  
  assign(paste0("p_",thresh),ggplot(d_tipping%>%filter(., Branch=="Degradation"))+
    geom_bar(aes(x=Nsp,fill=as.factor(Nb_different_species),y=Nb_different_species),
             position = "fill",stat="identity",alpha=.75)+
    labs(x="Number of species",y="Fraction of random initial conditions",fill="Number of shifts")+
    scale_x_continuous(breaks = seq(5,35,by=5))+
    the_theme+
      ggtitle(paste0("Threshold = ",thresh))
    )
    
    
}

p_tot=ggarrange(p_1,p_2,p_3,p_4,labels = letters[1:4],ncol = 2,nrow = 2,align = "hv")
ggsave("../Figures/Final_figs/SI/Threshold_tipping.pdf",p_tot,width = 9,height = 10)


## Figure N_species restoration ----


d_tipping=read.table("../Table/N_species/MF/Multistability_tipping.csv",sep=";")
d_richness2=read.table("../Table/N_species/MF/Multistability_richness2.csv",sep=";")
d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")


p1=ggplot(d_richness%>%filter(., Branch=="Restoration"))+
  geom_bar(aes(x=Nsp,fill=as.factor(Richness),y=Richness),
           position = "fill",stat="identity",alpha=.75)+
  labs(x="Community richness (N)",y="Fraction of random initial conditions",fill="Number of species")+
  scale_x_continuous(breaks = seq(5,35,by=5))+
  geom_point(data=d_richness2,aes(x=Nsp,y=Richness_D/Nsp),shape=1,size=3)+
  geom_line(data=d_richness2,aes(x=Nsp,y=Richness_D/Nsp))+
  scale_y_continuous(sec.axis=sec_axis(~., name="Fraction of species observed \n all across replicates"))+
  the_theme



d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")
colnames(d_tot)[11]="Community richness (N)"

p2=ggplot(d_tot%>%filter(.,Branch==2,`Community richness (N)` %in% c(5,35)))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.8,shape=21)+
  facet_wrap(.~`Community richness (N)`,labeller = label_both)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())

d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=150
colnames(d_ASS_sp)[6]="Community richness (N)"

#Species-specific tipping points size & frequency
p3_1=ggplot(filter(d_ASS_sp, `Community richness (N)` %in% c(5),Branch=="Restoration"), 
            aes(x=Trait,y=Tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .1)+
  facet_grid(.~`Community richness (N)`)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Frequency of abrupt shift",fill="")+
  scale_fill_gradientn(colours = color_Nsp(100))+
  scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))+
  theme(strip.text.x =element_blank(),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())

p3_2=ggplot(filter(d_ASS_sp, `Community richness (N)` %in% c(35),Branch=="Restoration"), 
            aes(x=Trait,y=Tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .025)+
  facet_grid(.~`Community richness (N)`)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Frequency of abrupt shift",fill="")+
  scale_fill_gradientn(colours = color_Nsp(100))+
  scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))+
  theme(strip.text.x =element_blank(),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank(),
        axis.line.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank(),
        axis.title.y = element_blank())

p3=ggarrange(p3_1,p3_2,ncol = 2,legend = "none",widths = c(1.2,1),labels = c(letters[3],""),hjust = -3,vjust = -2)

Fig6SI=ggarrange(p1,ggarrange(p2,p3,nrow=2,labels = c(letters[2],""),hjust=-2.3,
                            common.legend = T,legend = "bottom",heights = c(1.2,1)),
               ncol=2,labels = c(letters[1],"",""),widths = c(1.3,1))

ggsave("../Figures/Final_figs/SI/N_species_restoration.pdf",Fig6SI,width=10,height=5)


## All CSI panel ----


d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")
colnames(d_tot)[11]="Community richness (N)"

p=ggplot(d_tot%>%filter(.,Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  facet_wrap(.~`Community richness (N)`,labeller = label_both)+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())




ggsave("../Figures/Final_figs/SI/CSI_Nsp_all.pdf",p,width=7,height=5)

## Test: nb ASS along stress gradient ----

d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

d_ASS=tibble()
for (branch in 1:2){
  for (s in unique(d_tot$Stress)){
    for (com in unique(d_tot$Nsp)){
      
      d_fil=filter(d_tot,Stress==s,Nsp==com,Branch==branch)
      d_ASS=rbind(d_ASS,tibble(Nsp=com,Stress=s,Nb_ASS=length(unique(round(d_fil$CSI,2))),Branch=ifelse(branch==1,"Degradation","Restoration")))
      
    }
  }
}

ggplot(d_ASS)+
  geom_smooth(aes(x=Stress,y=Nb_ASS,color=as.factor(Nsp)),se = F)+
  the_theme+
  facet_grid(.~Branch)+
  scale_color_manual(values = colorRampPalette(c("orange","black","green"))(7))

