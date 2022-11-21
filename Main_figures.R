source('./Dryland_shift_functions.R')
library(grid)




# Figure 2: Scale facilitation, dispersal, multistability ----

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


scale="global"
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

set.seed(432)
rand_vec=runif(2)

d2$CSI=d2$Stress_tolerant*rand_vec[1]+d2$Competitive*rand_vec[2]

p1=ggplot(d2t%>%
            mutate(.,all_state=recode_factor(all_state,
                                           "Desert/Coexistence"="Coexistence/Desert",
                                           "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
                                           "Stress_tolerant"="Stress-tolerant",
                                           "Stress_tolerant/Desert"="Stress-tolerant/Desert"))) +
  geom_tile(aes(x=Stress,y=alpha_0,fill=all_state))+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of interspecific competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Competitive/Coexistence" = "#DDEFCA",
                             "Stress-tolerant" = "#7BD8D3",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(Scena~Delta,labeller=label_bquote(cols = delta == .(Delta)))+
  the_theme+theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=11))+
  geom_hline(data=tibble(Scena="Local facilitation",Delta=.1,y=c(0,.3)),aes(yintercept=y))+
  geom_text(data= tibble(Scena="Local facilitation",Delta=.1,y=c(0.015,.315),x=1.05,labels=c("b","c")), 
            aes(x=x, y=y,label=labels), alpha=1, colour="black")


p1_pattern=ggplot(d2t%>%
                    mutate(.,all_state=recode_factor(all_state,
                                                     "Desert/Coexistence"="Coexistence/Desert",
                                                     "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
                                                     "Stress_tolerant"="Stress-tolerant",
                                                     "Stress_tolerant/Desert"="Stress-tolerant/Desert"))) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
            colour          = 'transparent',
            pattern_spacing = 0.05, 
            pattern_density = 0.01, 
            pattern_fill    = "black",
            pattern_colour  = "black")+
  
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of interspecific competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Competitive/Coexistence" = "#DDEFCA",
                             "Stress-tolerant" = "#7BD8D3",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(Scena~Delta,labeller=label_bquote(cols = delta == .(Delta)))+
  the_theme+theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=11))+
  geom_hline(data=tibble(Scena="Local facilitation",Delta=.1,y=c(0,.3)),aes(yintercept=y))+
  geom_text(data= tibble(Scena="Local facilitation",Delta=.1,y=c(0.015,.315),x=1.05,labels=c("b","c")), 
            aes(x=x, y=y,label=labels), alpha=1, colour="black")+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))


y_line_color=c(0.48,.14,.80,.22)
for (i in 1:4){
  if (i%%2 ==0){
    
    assign(paste0("p2_",i),
           ggplot(d2%>%melt(., measure.vars=c("CSI"))%>%
                    filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,.3)[abs(i/2)])%>%
                    mutate(.,variable=recode_factor(variable,"CSI"="Global cover")))+
             
             geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
             facet_grid(.~alpha_0,labeller=label_bquote(cols=alpha[e]==.(alpha_0)))+
             the_theme+theme(strip.text.x = element_blank(),strip.background.x = element_blank())+
             scale_color_manual(values="black")+
             labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+guides(color=F)+
             scale_x_continuous(breaks = c(0,.5,1))+
             new_scale_color() +
             geom_point(data=d2t%>%
                          mutate(all_state=recode_factor(all_state,
                                                         "Desert/Coexistence"="Coexistence/Desert",
                                                         "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                                         "Desert/Stress_tolerant"="Stress_tolerant/Desert"))%>%
                          filter(., Scena=="Local facilitation",Delta==.1,alpha_0==c(0,.3)[abs(i/2)])%>%
                          mutate(., y=y_line_color[i]),
                        aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=3)+
             scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                                        "Competitive" = "#ACD87B", 
                                        "Competitive/Coexistence" = "#DDEFCA",
                                        "Stress_tolerant" = "#7BD8D3",
                                        "Stress_tolerant/Desert" ="#0F8E87",
                                        "Coexistence/Stress_tolerant"="#9BBBB9",
                                        "Coexistence/Desert"="#C19E5E",
                                        "Desert"=  "#696969"
             ))+
             scale_color_manual(values=c("Coexistence" = "#D8CC7B", 
                                         "Competitive" = "#ACD87B", 
                                         "Competitive/Coexistence" = "#DDEFCA",
                                         "Stress_tolerant" = "#7BD8D3",
                                         "Stress_tolerant/Desert" ="#0F8E87",
                                         "Coexistence/Stress_tolerant"="#9BBBB9",
                                         "Coexistence/Desert"="#C19E5E",
                                         "Desert"=  "#696969"
             ))
           
    )

    
  }else{
    
    
    assign(paste0("p2_",i),
           ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                    filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,0,.3)[i])%>%
                    mutate(.,variable=recode_factor(variable,"Stress_tolerant"="Stress-tolerant")))+
             
             geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
             facet_grid(.~alpha_0,labeller=label_bquote(cols=alpha[e]==.(alpha_0)))+
             the_theme+theme(strip.text.x = element_blank(),strip.background.x = element_blank())+
             scale_color_manual(values=c(as.character(color_rho[c(4)]),as.character(color_rho[c(2)])))+
             labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+guides(color=F)+
             scale_x_continuous(breaks = c(0,.5,1))+
             new_scale_color() +
             geom_point(data=d2t%>%
                          mutate(all_state=recode_factor(all_state,
                                                         "Desert/Coexistence"="Coexistence/Desert",
                                                         "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                                         "Desert/Stress_tolerant"="Stress_tolerant/Desert"))%>%
                          filter(., Scena=="Local facilitation",Delta==.1,alpha_0==c(0,0,.3)[i])%>%
                          mutate(., y=y_line_color[i]),
                        aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=3)+
             scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                                        "Competitive" = "#ACD87B", 
                                        "Competitive/Coexistence" = "#DDEFCA",
                                        "Stress_tolerant" = "#7BD8D3",
                                        "Stress_tolerant/Desert" ="#0F8E87",
                                        "Coexistence/Stress_tolerant"="#9BBBB9",
                                        "Coexistence/Desert"="#C19E5E",
                                        "Desert"=  "#696969"
             ))+
             scale_color_manual(values=c("Coexistence" = "#D8CC7B", 
                                         "Competitive" = "#ACD87B", 
                                         "Competitive/Coexistence" = "#DDEFCA",
                                         "Stress_tolerant" = "#7BD8D3",
                                         "Stress_tolerant/Desert" ="#0F8E87",
                                         "Coexistence/Stress_tolerant"="#9BBBB9",
                                         "Coexistence/Desert"="#C19E5E",
                                         "Desert"=  "#696969"
             ))
           
           
    )

    
  }
}


p_bottom=ggarrange(p2_1+ylab("Cover")+guides(color="none",fill="none"),p2_2+guides(color="none",fill="none")+ylab("Community index"),
                   p2_3+ylab("Cover")+guides(color="none",fill="none"),p2_4+guides(color="none",fill="none")+ylab("Community index"),
                   common.legend = T,legend = "bottom",ncol = 4,align="hv",labels=c(letters[2],"",letters[3],""))

Fig_2=ggarrange(p1,p_bottom,nrow=2,heights = c(3,1),labels = c(letters[1],""))

ggsave("../Figures/Final_figs/Figure_2.pdf",Fig_2,width = 8,height = 9)


Fig_2=ggarrange(p1_pattern+guides(pattern=F),p_bottom,nrow=2,heights = c(3,1),labels = c(letters[1],""))

ggsave("../Figures/Final_figs/Figure_2_pattern.pdf",Fig_2,width = 8,height = 9)


# Figure 3: Trait variation and multistability ----

c_inter_seq=c(0, .2, .3)
psi1_seq=c(1,0)
disp=.1
f=0.9

d=tibble()  
for (Psi_sp1 in c(0,1)){
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


for (Psi_sp1 in c(0,1)){

  if (Psi_sp1 ==1) {
    
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
                                     "Species 2/Coexistence"="Coexistence/Competitive",
                                     "Desert/Species 2"="Competitive/Desert",
                                     "Desert/Coexistence"="Coexistence/Desert",
                                     "Stress_tolerant"="Stress-tolerant",
                                     "Desert/Stress_tolerant"="Stress-tolerant/Desert",
                                     "Stress_tolerant/Coexistence"="Coexistence/Stress-tolerant",
                                     "Stress_tolerant/Species 2"="Competitive/Stress-tolerant",
                                     "Species 2/Coexistence"="Coexistence/Competitive",
                                     "Coexistence/Species 2"="Coexistence/Competitive",
                                     "Species 2" = "Competitive"))
      
    color_multistability=c("Coexistence" = "#D8CC7B",
                           "Competitive" = "#ACD87B",
                           "Coexistence/Competitive" = "#DDEFCA",
                           "Stress-tolerant" = "#7BD8D3",
                           "Stress-tolerant/Desert" ="#0F8E87",
                           "Coexistence/Stress-tolerant"="#9BBBB9",
                           "Coexistence/Desert"="#C19E5E",
                           "Desert"=  "#696969",
                           "Competitive/Stress-tolerant" = "#C998CE")
    #State multistability
    p1=ggplot(d_state) +
      geom_tile(aes(x=Stress,y=as.numeric(Psi2),fill=all_state))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait competitive sp., \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=color_multistability)+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      geom_hline(data=tibble(alpha_0=.2,y=c(.75)),aes(yintercept=y))+
      geom_text(data= tibble(alpha_0=.2,y=c(.8),x=1.05,labels=c("c")), 
                aes(x=x, y=y,label=labels), alpha=1, colour="black")
      
    
    #State multistability same but with patterns
    
    p1_pattern=ggplot(d_state) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = "gray20",
                        pattern_colour  = "gray20")+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait competitive sp., \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=color_multistability)+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      geom_hline(data=tibble(alpha_0=.2,y=c(.8)),aes(yintercept=y))+
      geom_text(data= tibble(alpha_0=.2,y=c(.85),x=1.05,labels=c("c")), 
                aes(x=x, y=y,label=labels), alpha=1, colour="black")+
      scale_pattern_manual(values=rev(c("none" ,"stripe")))

    
    #Bifurcation diagram
    
    p2_1=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
               filter(., Psi2 %in% unique(d_state$Psi2)[c(80)])%>%
               mutate(Psi2=round(abs(Psi2),2))%>%
               filter(., alpha_0==.2))+
      geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Cover", color = "",linetype="") +
      the_theme +
      scale_color_manual(values = color_rho[c(2, 4)],labels=c("Species 2","Stress-tolerant")) +
      scale_linetype_manual(values=c(1,9))+guides(color=F)+
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))+
      new_scale_color() +
      geom_point(data=d_state%>%
                   filter(., Psi2==unique(d_state$Psi2)[c(80)],alpha_0==.2)%>%
                   mutate(., y=.83),
                 aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=3)+
      scale_fill_manual(values=color_multistability)+
      scale_color_manual(values=color_multistability)
    
    
    
    #Community index
    p2_2=ggplot(d2%>%
                 filter(., Psi2 %in% unique(d_state$Psi2)[c(80)])%>%
                 mutate(Psi2=round(abs(Psi2),2))%>%
                 filter(., alpha_0==.2))+
      geom_line(aes(x = Stress, y = CSI,linetype=Branches)) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Community index", color = "",linetype="") +
      the_theme +
      theme(axis.text.x = element_text(size = 9),
            axis.text.y = element_text(size = 9),
            strip.text.y = element_text(size=9))+
      new_scale_color() +
      geom_point(data=d_state%>%
                   filter(., Psi2==unique(d_state$Psi2)[c(80)],alpha_0==.2)%>%
                   mutate(., y=.65),
                 aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=3)+
      scale_fill_manual(values=color_multistability)+
      scale_color_manual(values=color_multistability)
    
    
  
    
    
    
    
    
  } 
  if (Psi_sp1==0){
    
    d2=filter(d,Psi1==0)
    
    
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
    
    d_state=d_state%>%
      mutate(all_state=recode_factor(all_state,
                                     "Species 2/Coexistence"="Coexistence/Stress-tolerant",
                                     "Desert/Species 2"="Stress-tolerant/Desert",
                                     "Desert/Coexistence"="Coexistence/Desert",
                                     "Desert/Competitive"="Competitive/Desert",
                                     "Species 2" = "Stress-tolerant",
                                     "Competitive/Coexistence" = "Coexistence/Competitive"))
                                     
    
    color_multistability=c("Coexistence" = "#D8CC7B",
                           "Stress-tolerant" = "#7BD8D3",
                           "Coexistence/Stress-tolerant" = "#9BBBB9",
                           "Stress-tolerant/Desert" ="#0F8E87",
                           "Coexistence/Desert"="#C19E5E",
                           "Desert"=  "#696969")
    p2=ggplot(d_state%>%
                mutate(.,Psi2=round(Psi2,5),
                Stress=round(Stress,5))) +
      geom_tile(aes(x=Stress,y=abs(as.numeric(Psi2)),fill=all_state))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait stress-tolerant sp., \ $\psi_1)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=color_multistability)+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      geom_hline(data=tibble(alpha_0=.3,y=c(.5)),aes(yintercept=y))+
      geom_text(data= tibble(alpha_0=.3,y=c(0.55),x=1.05,labels=c("d")), 
                aes(x=x, y=y,label=labels), alpha=1, colour="black")
    
    
    p2_pattern=ggplot(d_state%>%
                mutate(Psi2=round(Psi2,5),Stress=round(Stress,5))) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = "gray20",
                        pattern_colour  = "gray20")+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait stress-tolerant sp., \ $\psi_1)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=color_multistability)+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      geom_hline(data=tibble(alpha_0=.3,y=c(.5)),aes(yintercept=y))+
      geom_text(data= tibble(alpha_0=.3,y=c(0.55),x=1.05,labels=c("d")), 
                aes(x=x, y=y,label=labels), alpha=1, colour="black")+
      scale_pattern_manual(values=rev(c("none" ,"stripe")))

    
    
    
    #Bifurcation diagrams
    d2t=transform(d2,
                  alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                          "alpha[e] == 0.5")),
                  Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                              labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
    
    
    p2_3=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
               filter(., Psi2 %in% unique(d_state$Psi2)[c(50)])%>%
               mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
               filter(., alpha_0==.3))+
      geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Cover", color = "",linetype="") +
      the_theme +
      scale_color_manual(values = (as.character(color_rho[c(2, 4)])),labels=c("Competitive","Species 2")) +
      scale_linetype_manual(values=c(1,9))+guides(color=F)+
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))+
      new_scale_color() +
      geom_point(data=d_state%>%
                   filter(., Psi2==unique(d_state$Psi2)[c(50)],alpha_0==.3)%>%
                   mutate(., y=.80),
                 aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=3)+
      scale_fill_manual(values=color_multistability)+
      scale_color_manual(values=color_multistability)
    
    
    
    
    
    
    #COmmunity index
    p2_4=ggplot(d2%>%
               filter(., Psi2 %in% unique(d2$Psi2)[c(50)])%>%
               mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
               filter(., alpha_0==.3))+
      geom_line(aes(x = Stress, y = CSI,linetype=Branches)) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Community index", color = "",linetype="") +
      the_theme +
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))+
      new_scale_color() +
      geom_point(data=d_state%>%
                   filter(., Psi2==unique(d_state$Psi2)[c(50)],alpha_0==.3)%>%
                   mutate(., y=.62),
                 aes(x=Stress,y=y,fill=all_state,color=all_state),alpha=.5,shape=22,size=3)+
      scale_fill_manual(values=color_multistability)+
      scale_color_manual(values=color_multistability)
    
    
    
  }
}

Figure_3_top=ggarrange(p1+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank(),axis.title.y.right = element_text(size = 11))+
            scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                                   name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
              ggtitle(TeX(r'($\psi_1 = 1)')),
          p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank(),axis.title.y.right = element_text(size = 11))+
            scale_y_continuous(sec.axis = sec_axis(trans = ~ .x ,
                                                   name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
            ggtitle(TeX(r'($\psi_2 = 0)')),common.legend = T,legend = "bottom",nrow=2,
          labels = letters[1:2])

Figure_3_bottom=ggarrange(p2_1+guides(color="none",fill="none")+theme(legend.position = "none"),
                          p2_2+guides(color="none",fill="none"),
                          p2_3+guides(color="none",fill="none")+theme(legend.position = "none"),
                          p2_4+guides(color="none",fill="none"),
                          ncol=4,labels = c(letters[3],"",letters[4],""),common.legend = T,legend = "bottom")

Figure_3_tot=ggarrange(Figure_3_top,Figure_3_bottom,nrow=2,heights = c(3,1))

ggsave("../Figures/Final_figs/Figure_3.pdf",Figure_3_tot,width = 9,height = 9)


#same but with patterns
Figure_3_top=ggarrange(p1_pattern+guides(pattern=F)+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                                axis.line.x = element_blank(),
                                axis.ticks.x = element_blank(),axis.title.y.right = element_text(size = 11))+
                         scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                                                name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
                         ggtitle(TeX(r'($\psi_1 = 1)')),
                       p2_pattern+guides(pattern=F)+theme(strip.background.x = element_blank(),strip.text.x = element_blank(),axis.title.y.right = element_text(size = 11))+
                         scale_y_continuous(sec.axis = sec_axis(trans = ~ .x ,
                                                                name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
                         ggtitle(TeX(r'($\psi_2 = 0)')),common.legend = T,legend = "bottom",nrow=2,
                       labels = letters[1:2])

Figure_3_tot=ggarrange(Figure_3_top,Figure_3_bottom,nrow=2,heights = c(3,1))

ggsave("../Figures/Final_figs/Figure_3_pattern.pdf",Figure_3_tot,width = 9,height = 9)


# Figure 4: Position of tipping points ----

c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)

f=.9;disp=.1
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
    
    
    d2=filter(d,Psi1==Psi_sp1)
    
    
    
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

p_tot=ggarrange(p1+theme(legend.position = "none",strip.text.x = element_text(size=13))+ggtitle(TeX("$\\psi_1 = 1$")),
                p2+ggtitle(TeX("$\\psi_2 = 0$"))+theme(strip.background.x = element_blank(),
                                                       strip.text.x = element_blank()),nrow=2,labels = letters[1:2])

ggsave("../Figures/Final_figs/Figure_4.pdf",p_tot,width = 7,height = 7)

# Figure 5: Niche consequences ----

#Describing how to compute the niche


d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                     .2,"_branch_","Degradation",
                     "_Psi1_",1,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";")%>%
  filter(., Psi2==unique(.$Psi2)[80])

p1=ggplot(d2%>%
         melt(., measure.vars=c("Stress_tolerant","Competitive") ))+
  geom_line(aes(x=Stress,y=value,color=variable))+
  the_theme+labs(x="Stress, S",y="Cover")+
  scale_color_manual(values = as.character(rev(color_rho[c(2, 4)])),labels=c("Competitive","Stress-tolerant"))+
  theme(legend.position = "none",plot.title = element_text(hjust = 0.5))+
  geom_line(data=tibble(x=c(0,d2$Stress[min(which(d2$Competitive==0))-2],0,0,d2$Stress[min(which(d2$Competitive==0))-2],d2$Stress[min(which(d2$Competitive==0))-2]),
                        y=c(.03,.03,.01,.05,.01,.05),
                        id=c(1,1,2,2,3,3)),
            aes(x=x,y=y,group=id),lwd=1)+
  geom_text(data=tibble(x=(d2$Stress[min(which(d2$Competitive==0))-2])/2,y=.1,text="Niche width (Nw)"),aes(x=x,y=y,label=text),size=4)+
  ggtitle(TeX("$\\Delta Nw_2=\\frac{Nw_{1+2}-Nw_2}{Nw_2}$"))

  
  




d[,1:2][d[,1:2] < 10^-4] = 0



#Niche for fixed traits
d_niche=read.table("../Table/2_species/PA/Niche_expansion_PA.csv",sep=";")
alpha_seq=seq(0,.3,length.out=10)
d_niche$alpha_0=alpha_seq
scale_facil="global_C_local_F"
p2=ggplot(d_niche%>%
            melt(.,measure.vars=c("Delta_niche_2"))%>%
            filter(., Scena==scale_facil,Local_disp==.1))+
  geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
  the_theme+
  scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
  labs(x=TeX(r'(Facilitation \ $f_0)'),y=TeX(r'(Competition strength \ $\alpha_e)'),fill="% of competitive \n species niche change")


#niche varying traits
d_niche = read.table("../Table/2_species/PA/Niche_expansion_varying_traits.csv",sep=";")
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
             aes(x=Psi2,y=value,shape=Facilitation,color=alpha_0),fill="white",size=3)+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,8))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  labs(y = "% of species 2 niche change", x = TeX(r'(Trait species 1, \ $\psi_1)'),
       color=TeX(r'(Interspecific competition \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f_0)'))+ggtitle(TeX("$\\psi_2 = 0$"))+
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
             aes(x=Psi2,y=value,shape=Facilitation,color=alpha_0),fill="white",size=3)+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,8))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  labs(y = "% of species 2 niche change", x = TeX(r'(Trait species 2, \ $\psi_2)'),
       color=TeX(r'(Interspecific competition \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f_0)'))+ggtitle(TeX("$\\psi_1 = 1$"))+
  the_theme
p3=ggarrange(p3_2,p3_1,labels = letters[3:4],common.legend = T,legend = "bottom")

Figure_5 = ggarrange(ggarrange(ggarrange(p1,ggplot()+theme_void(),heights = c(1,.2),nrow = 2),p2,widths = c(1,1.2),ncol=2,labels=letters[1:2]),
                     p3,nrow=2,labels=c("",""))


ggsave("../Figures/Final_figs/Figure_5.pdf",Figure_5,width = 8,height = 7)


# Figure 6: Species diversity ----

d_tipping=read.table("../Table/N_species/MF/Multistability_tipping.csv",sep=";")
d_richness2=read.table("../Table/N_species/MF/Multistability_richness2.csv",sep=";")
d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")


p1=ggplot(d_richness%>%filter(., Branch=="Degradation"))+
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

p2=ggplot(d_tot%>%filter(.,Branch==1,`Community richness (N)` %in% c(5,35)))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
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
p3_1=ggplot(filter(d_ASS_sp, `Community richness (N)` %in% c(5),Branch=="Degradation"), 
            aes(x=Trait,y=Tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .1)+
  facet_grid(.~`Community richness (N)`)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Frequency of abrupt shift",fill="")+
  scale_fill_gradientn(colours = color_Nsp(100))+
  scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))+
  theme(strip.text.x =element_blank(),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())

p3_2=ggplot(filter(d_ASS_sp, `Community richness (N)` %in% c(35),Branch=="Degradation"), 
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

Fig6=ggarrange(p1,ggarrange(p2,p3,nrow=2,labels = c(letters[2],""),hjust=-2.3,
                            common.legend = T,legend = "bottom",heights = c(1.2,1)),
              ncol=2,labels = c(letters[1],"",""),widths = c(1.3,1))
ggsave("../Figures/Final_figs/Figure_6.pdf",Fig6,width=10,height=5)


