source('./N_species_Analysis_functions.R')
source('./2_species_Analysis_functions.R')
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
                             "Desert"=  "#696969"
  ),labels=c("Stress_tolerant"="Stress-tolerant","Stress_tolerant/Desert"="Stress-tolerant/Desert","Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant"))+
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




















## Clustering of species 1 along dispersal gradient ----
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




ggsave("../Figures/Final_figs/SI/Clustering_11_dispersal.pdf",p,width = 7,height = 4)


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
                                              "Species 2/Coexistence"="Coexistence/Species 2",
                                              "Desert/Species 2"="Species 2/Desert",
                                              "Desert/Coexistence"="Coexistence/Desert",
                                              "Desert/Stress_tolerant"="Species 1/Desert",
                                              "Stress_tolerant"= "Species 1",
                                              "Stress_tolerant/Coexistence"="Coexistence/Species 1",
                                              "Stress_tolerant/Species 2"="Species 2/Species 1",
                                              "Species 2/Coexistence"="Coexistence/Species 2"))) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),              
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = alpha("black",.7),
                        pattern_colour  = alpha("black",.7))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = "Stress (S)", y = TeX(r'(Trait species 2, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Species 2" = "#ACD87B",
                                 "Coexistence/Species 2" = "#DDEFCA",
                                 "Species 1" = "#7BD8D3",
                                 "Species 1/Desert" ="#0F8E87",
                                 "Coexistence/Species 1"="#9BBBB9",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Desert"=  "#696969",
                                 "Species 2/Species 1" = "#C998CE"),
                        labels=c("Coexistence","Species 2","Coexistence/Species 2","Species 1","Species 1/Desert",
                                 "Coexistence/Species 1",
                                 "Coexistence/Desert","Desert","Species 2/Species 1"))+
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
                                              "Species 2/Coexistence"="Coexistence/Species 1",
                                              "Desert/Species 2"="Species 1/Desert",
                                              "Desert/Coexistence"="Coexistence/Desert",
                                              "Desert/Competitive"="Species 2/Desert",
                                              "Competitive" = "Species 2",
                                              "Species 2" = "Species 1",
                                              "Coexistence/Competitive" = "Coexistence/Species 2"
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
      labs(x = "Stress (S)", y = TeX(r'(Trait species 1, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Species 1" = "#7BD8D3",
                                 "Species 2" = "#ACD87B",
                                 "Coexistence/Species 1" = "#9BBBB9",
                                 "Species 1/Desert" ="#0F8E87",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Coexistence/Species 2" = "#DDEFCA",
                                 "Coexistence/Species 1"="#9BBBB9",
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
p3=ggarrange(p3_1,p3_2,labels = letters[2:3],common.legend = T,legend = "bottom")

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


## Position of tipping points ----

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
      labs(x=TeX(r'(Trait species 2, \ $\psi_2)'),y="Tipping point",color="")+
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
      labs(x=TeX(r'(Trait species 1, \ $\psi_1)'),y="Tipping point",color="")+
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

ggsave("../Figures/Final_figs/SI/Tipping_points_traits.pdf",width = 7,height = 7)

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
trait=rev(seq(0,1,length.out=15))
comp=comp2=c()
for (i in 1:15){
  comp=c(comp,sum(sapply(1:15,function(x){
    exp(-abs(trait[i]-trait[x]))
  })))
  comp2=c(comp2,sum(sapply(1:15,function(x){
    trait[i]*exp(-abs(trait[i]-trait[x]))
  })))
}


p1=ggplot(tibble(Com=comp,Trait=trait,Comp2=comp2))+
  geom_line(aes(x=Trait,y=Com))+
  geom_point(aes(x=Trait,y=Com,color=Trait),size=3)+
  the_theme+
  labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $\sum_j e^{-|\psi_i-\psi_j|})'))+
  scale_color_gradientn(colours = color_Nsp(15))
p2=ggplot(tibble(Com=comp,Trait=trait,Comp2=comp2))+
  geom_line(aes(x=Trait,y=Comp2))+
  geom_point(aes(x=Trait,y=Comp2,color=Trait),size=3)+
  the_theme+
  labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $\psi_i \sum_j e^{-|\psi_i-\psi_j|})'))+
  scale_color_gradientn(colours = color_Nsp(15))

ggsave("../Figures/Final_figs/SI/Competition_experienced.pdf",
       ggarrange(p1,p2,ncol=2,common.legend = T,legend = "bottom",labels = letters[1:2]),
       width = 9,height = 4)



# N-species ----
## Size and species specific tipping points ----


d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=40


#Species-specific tipping points size & frequency
p2=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/(N_replicate*3),fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f[0] == .(Facilitation)))+
  the_theme+
  labs(x=TeX(r'(Species trait \ $\psi)'),y="Fraction of tipping points",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))


d_ASS_sp=d_ASS_sp%>%
  group_by(.,Facilitation,Competition,Branch,Trait)%>%
  summarise(.,Size_shift=mean(Size_shift),.groups = "keep")

p3=ggplot(d_ASS_sp, aes(x=Trait,y=Size_shift,fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f[0] == .(Facilitation)))+
  the_theme+
  labs(x=TeX(r'(Species trait \ $\psi)'),y="Mean shift size (in fraction of cover)",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))


Figure_4 = ggarrange(p2+theme(strip.background.y = element_blank(),strip.text.y = element_blank())
                               ,p3,ncol=2,labels=letters[1:2],common.legend = T,legend = "bottom")
ggsave("../Figures/Final_figs/SI/Specific_specific_tipping_points.pdf",Figure_4,width = 9,height = 5)


## Restoration trajectory----
d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")
Nsp=15
p1=ggplot(d%>%filter(.,Branch=="Restoration"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.8,shape=1)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f[0]==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(Nsp),na.value = "black")+
  guides(fill="none")

ggsave("../Figures/Final_figs/SI/Restoration_trajectory.pdf",p1,width = 7,height = 5)

## Varying the threshold for shifts ----

d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")

# At species level : mean number of shifts/replicate
treshold_seq=c(.03,.05,.07,.1)
N_replicate=40;Nsp=15;traits=rev(seq(0,1,length.out=15));treshold_com=.1

for (i in 1:length(treshold_seq)){

  treshold_sp=treshold_seq[i]


  d_ASS_sp=d_ASS_com=d_undetected=tibble()
  for (f in unique(d$Facilitation)){
    for (a0 in unique(d$Competition)){
      for (branch in unique(d$Branch)){
        for (ini in unique(d$Random_ini)){
          
          #First species level
          
          d_fil=filter(d,Facilitation==f,Competition==a0,Branch==branch,Random_ini==ini)%>%
            melt(.,  measure.vars=paste0("Sp_",1:15))
          
          d_fil$value[d_fil$value<10^(-4)]=0
          
          # ggplot(d_fil)+
          #   geom_line(aes(x=Stress,y=value,color=variable))
          # ggplot(d_fil)+
          #   geom_line(aes(x=Stress,y=Rho_p))
          
          # d_shifts=tibble() #To compare when does community versus species shift respectively
          
          for (sp in 1:Nsp){ #we look for each species whether it shifted abruptly or not
            
            d_sp=filter(d_fil,variable==paste0("Sp_",sp)) 
            
            
            if (any(d_sp$value>0)){# we only keep species that are present in the system
              
              tipping=ifelse(any(abs(diff(d_sp$value))>treshold_sp),1,0) #Is there a tipping point in the system ?
              
              #for the size of the shift, we are interested in the shift were species disappear. 
              
              if (sum(abs(diff(d_sp$value))>treshold_sp)==1){ #only one shift
                size_shift=max(abs(diff(d_sp$value)))
                # d_shifts=rbind(d_shifts,tibble(Stress=d_sp$Stress[which(abs(diff(d_sp$value)) == max(abs(diff(d_sp$value))))], Species=sp))
                
              }else if (sum(abs(diff(d_sp$value))>treshold_sp)==0){ #0 shifts
                size_shift=0
                
              } else { #multiple shifts
                
                if (branch=="Degradation"){ # we take the last shift
                  size_shift=abs(diff(d_sp$value))[which(abs(diff(d_sp$value))>treshold_sp)[length(which(abs(diff(d_sp$value))>treshold_sp))]]
                  # d_shifts=rbind(d_shifts,tibble(Stress=d_sp$Stress[which(abs(diff(d_sp$value)) == abs(diff(d_sp$value))[which(abs(diff(d_sp$value))>treshold_sp)[length(which(abs(diff(d_sp$value))>treshold_sp))]])], Species=sp))
                } else{
                  size_shift=abs(diff(d_sp$value))[which(abs(diff(d_sp$value))>treshold_sp)][1]
                  # d_shifts=rbind(d_shifts,tibble(Stress=d_sp$Stress[which(abs(diff(d_sp$value)) == abs(diff(d_sp$value))[which(abs(diff(d_sp$value))>treshold_sp)][1])], Species=sp))
                  
                }
              }
              #Even if there was no shifts, we save the tipping event
              
              d_ASS_sp=rbind(d_ASS_sp,tibble(Tipping=tipping,Species=sp,Trait=traits[sp],Facilitation=f,Competition=a0,Branch=branch,Random_ini=ini,Size_shift=size_shift))
              
            }
          }
          
          #Secondly community level
          
          # d_fil=filter(d,Facilitation==f,Competition==a0,Branch==branch,Random_ini==ini)%>%
          #   melt(.,  measure.vars="Rho_p")
          # 
          # d_fil$value[d_fil$value<10^(-4)]=0
          # 
          # 
          # tipping=ifelse(any(abs(diff(d_fil$value))>treshold_com),1,0)
          # 
          # if (tipping==1){ #there is shift one or multiple
          #   for (nb_shift in 1:length(which(abs(diff(d_fil$value))>treshold_com))){ #for each community shift, we get its size, and save into the DF
          #     size_shift = abs(diff(d_fil$value))[which(abs(diff(d_fil$value))>treshold_com)[nb_shift]]
          #     d_ASS_com=rbind(d_ASS_com,tibble(Tipping=tipping,Facilitation=f,Competition=a0,Branch=branch,Random_ini=ini,Size_shift=size_shift))
          #   }
          #   # 
          #   # #We track species that are not detected at the community level
          #   # stress_tipping_com=d_fil$Stress[which(abs(diff(d_fil$value))>treshold_com)]
          #   # 
          #   # if (length(stress_tipping_com) < nrow(d_shifts) |  length(unique(d_shifts$Stress))!=1 & unique(d_shifts$Stress)!=stress_tipping_com) { #i.e. we detect less shifts at community scale compared to species scale or all species shifted at the same threshold
          #   #   (1:nrow(d_shifts))[-which(d_shifts$Stress %in% stress_tipping_com)]
          #   # }
          #   # 
          #   
          #   
          # } else{ #no shifts
          #   d_ASS_com=rbind(d_ASS_com,tibble(Tipping=tipping,Facilitation=f,Competition=a0,Branch=branch,Random_ini=ini,Size_shift=NA))
          # }
          
          
          
        }
      }
    }
  }
  
  #Species-specific tipping points size & frequency
  assign(paste0("p_",i),
         ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/(N_replicate*3),fill=as.factor(Competition))) +
           geom_bar(position="stack", stat="identity")+
           facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f == .(Facilitation)))+
           the_theme+
           labs(x=TeX("$\\psi$"),y="Fraction of tipping points",fill=TeX("$\\alpha_e$ \ \ \ "))+
           scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
           ggtitle(paste0("Threshold = ",treshold_seq[i]))
         
  )
  
}
p_tot=ggarrange(p_1,p_2,p_3,p_4,ncol=2,nrow=2,common.legend = T,legend = "bottom")
ggsave("../Figures/Final_figs/SI/Threshold_tipping_species.pdf",p_tot,width = 11,height = 9)

