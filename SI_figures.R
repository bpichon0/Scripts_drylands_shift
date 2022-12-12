source('./Dryland_shift_functions.R')
library(grid)






# 2-species ----
## Multistability fixed traits MF model ----

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




















## Multistability varying traits MF model ----

stress_seq=c(seq(0,.82, length.out = 100),seq(0.82,.9, length.out = 10)[-1])

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/MF/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",stress,"_facilitation_",.9,".csv"),sep=";"),
           read.table(paste0("../Table/2_species/MF/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",stress,"_facilitation_",.9,".csv"),sep=";"))
  d=rbind(d,d2)
  
  
  
}

d2=d
d2[,1:2][d2[,1:2] < 10^-3] = 0
d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Sp1")
  }
  if (d2[x, 1] == 0 & d2[x, 2] > 0) {
    return("Sp2")
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
  if (d_state$all_state[x] %in% c( "Sp1/Sp2","Coexistence/Sp1","Sp2/Coexistence",
                                   "Sp2/Sp1")){
    return(1)
  } else {return(0)}
  
})

d_state$Sp1_reg=d2$Sp1[which(d2$Branches=="Restoration")]
d_state$Sp2_reg=d2$Sp2[which(d2$Branches=="Restoration")]

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      
      if (i>j){
        d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,
                                     Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                     Frac_multi_raw=sum(d_fil$multistab),
                                     Deg=d_fil$Stress[min(which(d_fil$Sp2==0))-1],
                                     Rest=d_fil$Stress[min(which(d_fil$Sp2_reg==0))]
        ))
      }
    }
  }
}

trait1_for_bifu=c(unique(d2$Psi1)[15],unique(d2$Psi1)[49],unique(d2$Psi1)[49])
trait2_for_bifu=c(unique(d2$Psi2)[13],unique(d2$Psi2)[1],unique(d2$Psi2)[49])

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),
               color="black",alpha=.3,breaks = c(.1,.15,.2))+
  geom_point(data=tibble(x=trait1_for_bifu,y=trait2_for_bifu),aes(x=x,y=y),shape=19)+
  geom_text(data=tibble(x=c(trait1_for_bifu[1],trait1_for_bifu[2]+.05,trait1_for_bifu[3]+.05),
                        y=c(trait2_for_bifu[1]+.13,trait2_for_bifu[2]+.05,trait2_for_bifu[3]+0.05),
                        txt=c("(ii)","(i)","(iii)")),aes(x=x,y=y,label=txt),size=4.5)+
  geom_segment(aes(x = 0.05, y = -0.02, xend = .95, yend = -0.02),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = .5, y = -0.06,label = "Sensitivity to competition"),stat = "unique")+
  geom_segment(aes(x = 1.02, y = 0, xend = 1.02, yend = .95),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = 1.05, y = .45,label = "Trait similarity"),stat = "unique",angle=270)+
  the_theme+
  labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
       y=TeX(r'(Species 2 trait, $\psi_2$)'),
       fill="Fraction of community bistability \n along the stress gradient")+
  scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100))+
  guides(shape="none")


for (i in 1:3){
  
  if (i %in% c(1,3)){
    assign(paste0("p2_",i),
           ggplot(d2%>%
                    filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i])%>%
                    melt(.,measure.vars=c("Sp1","Sp2"))%>%
                    mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
             geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
             scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank(),
                   legend.position = "none")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             new_scale_color() +
             
             geom_line(data=d_state%>%
                         filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i]),
                       aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
             scale_color_manual(values=c("white","black"))+
             guides(color=F)
    )
    
  }else {
    assign(paste0("p2_",i),
           ggplot(d2%>%
                    filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i])%>%
                    melt(.,measure.vars=c("Sp1","Sp2"))%>%
                    mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
             geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
             scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             theme(legend.text = element_text(size=11))+
             new_scale_color() +
             
             geom_line(data=d_state%>%
                         filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i]),
                       aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
             scale_color_manual(values=c("white","black"))+
             guides(color=F)
           
           
    )
    
  }
}





pbifu=ggarrange(p2_2,p2_1,p2_3,
                ncol=3,
                labels = c("(i)","(ii)","(iii)"),common.legend = T,legend = "bottom",widths = c(1.2,1,1),hjust = c(-5,-1.3,-1.7),vjust=4)
p_tot=ggarrange(
  ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.15,1,.15)),
  pbifu,nrow=2,heights =  c(2,1),labels = c(letters[1:2]),hjust = c(-9,0))

ggsave(filename = "../Figures/Final_figs/SI/Multistability_varying_traits_MF.pdf",p_tot,width = 7,height = 7)






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


Fig_2_SI=ggplot(d2t%>%
                    mutate(.,all_state=recode_factor(all_state,
                                                     "Competitive/Coexistence"="Coexistence/Competitive",
                                                     "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
                                                     "Stress_tolerant"="Stress-tolerant",
                                                     "Stress_tolerant/Desert"="Stress-tolerant/Desert"))%>%
                    mutate(., Delta=recode_factor(Delta,"0.1"="0.1, local dispersal","0.9"="0.9, global dispersal"))) +
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
                             "Desert"=  "#696969"))+
  facet_grid(Scena~Delta)+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9),strip.text.y = element_text(size=12))+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))+
  theme(legend.key.size = unit(2, 'cm'))


ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_local_competition.pdf",Fig_2_SI,width = 8,height = 6)




## Multistability along competition gradient with global facilitation ----

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
                                             "Stress_tolerant/Desert"="Stress-tolerant/Desert",
                                             "Competitive/Stress_tolerant"="Competitive/Stress-tolerant"))%>%
            mutate(., Delta=recode_factor(Delta,"0.1"="0.1, local dispersal","0.9"="0.9, global dispersal"))%>%
            filter(., Scena=="Global facilitation")) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',
                    pattern_spacing = 0.05, 
                    pattern_density = 0.01, 
                    pattern_fill    = "black",
                    pattern_colour  = "black")+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of interspecific competition, \ $\alpha_e)'), fill = "") +
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive/Stress-tolerant"="#C998CE",
                             "Competitive" = "#ACD87B", 
                             "Competitive/Coexistence" = "#DDEFCA",
                             "Stress-tolerant" = "#7BD8D3",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(.~Delta,labeller = label_bquote(cols= delta==.(Delta) ))+
  the_theme+theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=11))+
  theme(legend.text = element_text(size = 9))+
  geom_text(data= tibble(Scena="Local facilitation",Delta="0.1, local dispersal",y=c(0.015,.415),x=1.05,labels=c("b","c")), 
            aes(x=x, y=y,label=labels), alpha=1, colour="white")+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))+
  guides(pattern="none")

  


ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_global_facilitation.pdf",p1,width = 9,height = 5)



## Multistability, varying the trade-off shape ----

#first: explaining trade-off


d=expand.grid(Psi=seq(0,1,length.out=20),Gamma=c(.5,1,1.5))
d$Psi2=sapply(1:nrow(d),function(x){
  return(d$Psi[x]^d$Gamma[x])
})

p=ggplot(d)+
  geom_line(aes(x=Psi,y=Psi2,color=as.factor(Gamma)))+
  geom_point(aes(x=Psi,y=Psi2,color=as.factor(Gamma)),shape=1,size=3,fill="white")+
  the_theme+
  labs(x=TeX("$\\psi$"),y=TeX("$\\psi^{\\gamma}$"),color=TeX("$\\gamma$"))+
  scale_color_viridis_d()
ggsave( "../Figures/Final_figs/SI/Trade_off_shape.pdf",p,width = 6,height = 4)



#second: multistability
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()

for (tradeoff in c(.5,1.5)){
  for (stress in stress_seq){
    
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Degradation_stress_",stress,"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Restoration_stress_",stress,"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"))
    d=rbind(d,d2)
    
    
    
  }
}

d2=d
colnames(d2)[1:2]=c("Sp1","Sp2")
d2[,1:2][d2[,1:2] < 10^-2] = 0
d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Sp1")
  }
  if (d2[x, 1] == 0 & d2[x, 2] > 0) {
    return("Sp2")
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
  if (d_state$all_state[x] %in% c( "Sp1/Sp2","Coexistence/Sp1","Coexistence/Sp2",
                                   "Sp2/Coexistence","Sp1/Coexistence",
                                   "Sp2/Sp1")){
    return(1)
  } else {return(0)}
  
})

d_state$Sp1_reg=d2$Sp1[which(d2$Branches=="Restoration")]
d_state$Sp2_reg=d2$Sp2[which(d2$Branches=="Restoration")]

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      for (tradeoff in unique(d_state$Tradeoff)){
        d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0,Tradeoff==tradeoff)
        
        if (i>j){
          d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,Tradeoff=tradeoff,
                                       Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                       Frac_multi_raw=sum(d_fil$multistab),
                                       Deg=d_fil$Stress[min(which(d_fil$Sp2==0))-1],
                                       Rest=d_fil$Stress[min(which(d_fil$Sp2_reg==0))]
          ))
        }
      }
    }
  }
}

trait1_for_bifu=c(unique(d_state$Psi1)[49],unique(d_state$Psi1)[13])
trait2_for_bifu=c(unique(d_state$Psi2)[1],unique(d_state$Psi2)[12])

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),
               color="black",alpha=.3,breaks = c(0,.2,.4))+
  geom_point(data=tibble(x=trait1_for_bifu,y=trait2_for_bifu),aes(x=x,y=y),shape=19)+
  geom_text(data=tibble(x=c(trait1_for_bifu[1],trait1_for_bifu[2]+.02),
                        y=c(trait2_for_bifu[1]+.06,trait2_for_bifu[2]+.12),
                        txt=c("(i)","(ii)")),aes(x=x,y=y,label=txt),size=4.5)+
  the_theme+
  facet_wrap(.~Tradeoff,labeller = label_bquote(cols = gamma == .(Tradeoff)))+
  labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
       y=TeX(r'(Species 2 trait, $\psi_2$)'),
       fill="Fraction of community bistability \n along the stress gradient")+
  scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100))+
  guides(shape="none")+
  theme(strip.text.x = element_text(size=13))


for (i in 1:2){
  
  assign(paste0("p2_",i),
         ggplot(d2%>%
                  filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==0.5)%>%
                  melt(.,measure.vars=c("Sp1","Sp2"))%>%
                  mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
           geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
           scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
           the_theme+
           labs(x="Stress (S)",y="Species cover",color="",linetype="")+
           scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
           theme(legend.text = element_text(size=11))+
           new_scale_color() +
           
           geom_line(data=d_state%>%
                       filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==0.5),
                     aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
           scale_color_manual(values=c("white","black"))+
           guides(color=F)
  )
  
  assign(paste0("p3_",i),
         ggplot(d2%>%
                  filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==1.5)%>%
                  melt(.,measure.vars=c("Sp1","Sp2"))%>%
                  mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
           geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
           scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
           the_theme+
           labs(x="Stress (S)",y="Species cover",color="",linetype="")+
           scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
           theme(legend.text = element_text(size=11))+
           new_scale_color() +
           
           geom_line(data=d_state%>%
                       filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==1.5),
                     aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
           scale_color_manual(values=c("white","black"))+
           guides(color=F)
  )  
    
  
}





pbifu=ggarrange(ggarrange(p2_1,p2_2+ylab(""), labels = c("(i)","(ii)"),ncol=2,widths = c(1.2,1),hjust = c(-4.5,-3.5),vjust=4),
                ggarrange(p3_1+ylab(""),p3_2+ylab(""), labels = c("(i)","(ii)"),ncol=2,widths = c(1,1),hjust = c(-4.5,-6),vjust=4),
                ncol=2,
                common.legend = T,legend = "bottom",
                widths = c(1.1,1),labels = c("b","c"))
p_tot=ggarrange(
  ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.15,1,.15)),
  pbifu,nrow=2,heights =  c(2,1),labels = c(letters[1],""),hjust = c(-9,0))

ggsave(filename = "../Figures/Final_figs/SI/Trade_off_multistability_PA.pdf",plot = p_tot,width = 9,height = 7)



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
         filter(., alpha_0==.2, S %in% c(.73, .77),
                delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
         group_by(., cintra,alpha_0,S,delta,Scena)%>%
         summarise(.,.groups ="keep",c11=mean(c11) ))+
  geom_point(aes(x=Scena,y=c11,color=delta),size=3,shape=1)+
  geom_line(aes(x=Scena,y=c11,color=delta,group=interaction(delta,S)),lwd=.5)+
  facet_grid(.~S,labeller=label_bquote(cols = Stress == .(S)),scales = "free_y")+
  the_theme+labs(x="",color=TeX(r'( Dispersal scale, $\delta)'),y=expression(paste("Stress-tolerant clustering (c"[11],")")))+
  geom_hline(yintercept = 1)+
  scale_color_viridis_b(option = "A")

ggsave("../Figures/Final_figs/SI/Clustering_11_dispersal.pdf",p,width = 7,height = 4)

p=ggplot(d_clustering%>%
         filter(., alpha_0==.2, S %in% c(0),
                delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=10)])%>%
         melt(., measure.vars=c("Rho_20","Rho_10"))%>%
         mutate(.,variable=recode_factor(variable,"Rho_10"="Pair Stress-tolerant sp./Fertile",
                                         "Rho_20"="Pair Competitive sp./Fertile")))+
  geom_point(aes(x=Scena,y=value,color=delta),size=3,shape=1)+
  geom_line(aes(x=Scena,y=value,color=delta,group=interaction(delta,variable)),lwd=.5)+
  ggtitle("Stress (S) = 0")+
  facet_wrap(.~variable,scales="free")+
  the_theme+labs(color=TeX(r'( Dispersal scale, $\delta)'),y=expression(paste("Pair cover")),x="")+
  scale_color_viridis_b(option = "A")

ggsave("../Figures/Final_figs/SI/Pair_10_20.pdf",p,width = 7,height = 4)




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

## Threshold for invasion ----

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
N_rep = 100
S_seq = seq(0,1,length.out=1000)
alpha_seq = c(.2)
f_seq=.9
delta_seq=seq(0,1,length.out=10)
cintra_seq=c(.3)


name_scena=c("global_C_global_F","global_C_local_F")

d_invade=tibble() #initializing the tibble

for (scena_ID in 1:2){ #for each scenario of species pairs
  
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
          
          
          state=Get_PA_initial_state(ini =c(.05,.05,.49))
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
          
          
          
          d_invade=rbind(d_invade,tibble(
            Thresh_invasion=max(d2$S[which((d2$rho_1)>0)]),
            alpha_0 = alpha0,
            f=f,delta=disp,Scena=scena_ID,
            cintra=aii
          ))
          
          
        } #end competition loop
        
      } #end facilitation loop
      
    } #end h loop
    
  } #end dispersal loop
  
} #end scenario loop

p=ggplot(d_invade%>%
           mutate(., Scena=recode_factor(Scena,'1'="Global facilitation",'2'="Local facilitation")))+
  geom_point(aes(x=Scena,y=Thresh_invasion,color=delta),size=3,shape=1)+
  geom_line(aes(x=Scena,y=Thresh_invasion,color=delta,group=interaction(delta)),lwd=.5)+
  labs(x="",y="Threshold of stress for invasion \n of a desert state",shape="",color=TeX("$\\delta$"))+
  scale_shape_manual(values = c(0,8))+the_theme+
  scale_color_viridis_b(option = "A")+
  theme(legend.text = element_text(size=12))


ggsave("../Figures/Final_figs/SI/Threshold_invasion.pdf",width = 7,height = 4)



## Restoration point trait & competition species----


stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"),
           read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"))
  d=rbind(d,d2)
  
  
  
}

d2=d
d2[,1:2][d2[,1:2] < 10^-3] = 0
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
  if (d_state$all_state[x] %in% c( "Species 2/Stress_tolerant","Coexistence/Stress_tolerant","Species 2/Coexistence",
                                   "Stress_tolerant/Species 2")){
    return(1)
  } else {return(0)}
  
})

d_state$Stress_tolerant_reg=d2$Stress_tolerant[which(d2$Branches=="Restoration")]
d_state$Competitive_reg=d2$Competitive[which(d2$Branches=="Restoration")]

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      
      if (i>j){
        d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,
                                     Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                     Frac_multi_raw=sum(d_fil$multistab),
                                     Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
                                     Deg2=d_fil$Stress[max(which(d_fil$Stress_tolerant>0))+1],
                                     Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))],
                                     Niche_1 = length(which(d_fil$Stress_tolerant>0))/length(which(d_fil$state!="Desert"))  ,
                                     Niche_2 = length(which(d_fil$Competitive>0))/length(which(d_fil$state!="Desert"))
        ))
      }
    }
  }
}



# Niche least competitive species
p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Niche_1))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_2$)'),
                 fill="Fraction of vegetation states \n where species 1 has a positive cover")+
  scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100))+
  guides(shape=F)


ggsave("../Figures/Final_figs/SI/Niche_least_competitive_species.pdf",p1,width = 6,height = 5)


#Degradation and restoration points

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Deg))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_1$)'),fill="Species 2 extinction point")+
  scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100))+
  guides(shape=F)



#Restoration point traits
p2=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Rest))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_2$)'),fill="Species 2 restoration point")+
  scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100),
                       breaks=c(0,.1,.2))+
  geom_line(data=tibble(x=c(0,1,1,1),y=c(0,0,0,1),group=c(1,1,2,2)),aes(x=x,y=y,group=group),lwd=1)+
  geom_text(data=tibble(x=c(.5,1.05),
                        y=c(0.05,.55),
                        txt=c("(c)","(d)")),aes(x=x,y=y,label=txt),size=4.5)+
  guides(shape=F)



p3=ggplot(d_final%>%filter(Psi2==0))+
  geom_smooth(aes(x=Psi1,y=Rest),se = F,color="black")+ggtitle(TeX("$\\psi_2 = 0$"))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y="Species 2 restoration point")

p4=ggplot(d_final%>%filter(Psi1==unique(d_final$Psi1)[49]))+
  geom_smooth(aes(x=Psi2,y=Rest),se = F,color="black")+ggtitle(TeX("$\\psi_1 = 1$"))+
  the_theme+labs(x=TeX(r'(Species 2 trait, $\psi_2$)'),y="Species 2 restoration point")

p_bottom=ggarrange(p3,p4,labels = letters[3:4])
p_tot=ggarrange(ggarrange(p1,p2,ncol=2,labels = letters[1:2]),
                p_bottom,nrow=2,heights = c(1.5,1))

ggsave("../Figures/Final_figs/SI/Restoration_degration.pdf",p_tot,width = 9,height = 7)








# N-species ----


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




## Number of species gradient stress ----


d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")


d_richness_mean=d_richness%>%
  group_by(., Nsp,Branch,Competition)%>%
  summarise(., mean_richness=mean(Richness),.groups = "keep")

p1=ggplot(d_richness)+
  geom_bar(aes(x=Competition,fill=as.factor(Richness),color=as.factor(Richness),y=Richness),
           color="transparent",position = "fill",stat="identity")+
  geom_point(data=d_richness_mean,aes(x=Competition,y=mean_richness/7),
             shape=1,size=3,color="white")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),
       y="Fraction of replicates",fill="Number of species \n with positive cover")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  scale_y_continuous(sec.axis=sec_axis(~.*7, name="Mean total number of species"))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))


ggsave("../Figures/Final_figs/SI/Nb_shift_diversity.pdf",p1,width=9,height=5)
