source('./Dryland_shift_functions.R')
library(grid)


# Figure 2: Dynamics and landscapes ----

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


for (i in 1:2){
  assign(paste0("p2_",i),
         ggplot(d2%>%melt(., measure.vars=c("CSI"))%>%
                  filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,.4)[i])%>%
                  mutate(.,variable=recode_factor(variable,"CSI"="Global cover")))+
           
           geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
           the_theme+
           scale_color_manual(values="black")+
           labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+guides(color=F)+
           scale_x_continuous(breaks = c(0,.25,.5,.75,1))    )
  
  
  assign(paste0("p1_",i),
         ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                  filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,.4)[i])%>%
                  mutate(.,variable=recode_factor(variable,"Stress_tolerant"="Stress-tolerant")))+
           
           geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
           the_theme+theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())+
           scale_color_manual(values=c(as.character(color_rho[c(4)]),as.character(color_rho[c(2)])))+
           labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+
           scale_x_continuous(breaks = c(0,.25,.5,.75,1))    )

}



list_land=list.files("../Table/2_species/CA/Illustration/",pattern = "Landscape")
list_dim=list.files("../Table/2_species/CA/Illustration/",pattern = "Dynamics")
color_rho2 = c("Fertile" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress-tolerant" = "#7BD8D3")

for (i in 1:length(list_dim)){
  d3=read.table(paste0("../Table/2_species/CA/Illustration/",list_dim[i]),sep=",")
  colnames(d3) = c("Time","Stress-tolerant", "Competitive", "Fertile", "Desert")
  
  assign(paste0("p3_",i),ggplot(d3 %>% melt(., id.vars = "Time"))+
           geom_line(aes(x = Time, y = value, color = variable), lwd = 1) +
           the_theme +
           scale_color_manual(values = color_rho2) +
           labs(x = "Time", y = "Densities", color = "") +
           theme(legend.text = element_text(size = 11), legend.position = "bottom")+
           xlim(0,5000))
  
  landscape=read.table(paste0("../Table/2_species/CA/Illustration/",list_land[i]),sep=",")
  assign(paste0("p4_",i),Plot_landscape(landscape,Nsp=2))
}

p4_2=p4_2+scale_fill_gradientn(colours = c("white","white",color_rho2[4]))

p_top=ggarrange(p1_1+ylab("Species cover")+ggtitle(TeX("$\\alpha_e=0$")),
                p1_2+geom_vline(xintercept = 0.2)+ggtitle(TeX("$\\alpha_e=0.4$"))+geom_text(aes(x=.24,y=.8,label="c"),state="unique"),
                ncol=2,labels = c(letters[1],""),common.legend = T,legend = "bottom")

p_middle = ggarrange(p2_1+ylab("Community index"),p2_2,ncol=2,labels = c(letters[2],""),common.legend = T,legend = "none",hjust=-2,vjust=0)
p_bottom = ggarrange(ggarrange(p3_1+ylab("Cover"),
                               ggarrange(ggplot()+theme_void(),p4_1+guides(fill="none"),ggplot()+theme_void(),nrow=3,heights = c(.1,.9,.2)),
                               ncol=2,legend = "none",widths = c(1,.5)),
                     ggarrange(p3_2+ylab(""),
                               ggarrange(ggplot()+theme_void(),p4_2+guides(fill="none"),ggplot()+theme_void(),nrow=3,heights = c(.1,.9,.2)),
                               ncol=2,legend = "none",widths = c(1,.5)),
                     ncol=2,labels = c(letters[3],""),common.legend = T,legend = "none")

p_tot=ggarrange(p_top,p_middle,p_bottom,nrow=3,heights = c(1.1,1,1.1))

ggsave("../Figures/Final_figs/Figure_2.pdf",p_tot,width = 9,height = 7)





# Figure 3: Scale facilitation, dispersal, multistability ----

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
            filter(., Scena=="Local facilitation")) +
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
  geom_hline(data=tibble(Scena="Local facilitation",Delta="0.1, local dispersal",y=c(0,.4)),aes(yintercept=y))+
  geom_text(data= tibble(Scena="Local facilitation",Delta="0.1, local dispersal",y=c(0.015,.415),x=1.05,labels=c("b","c")), 
            aes(x=x, y=y,label=labels), alpha=1, colour="black")+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))


#XXXX mechanisms

Fig_2=ggarrange(p1_pattern+guides(pattern=F),p_bottom,nrow=2,heights = c(3,1.3),labels = c(letters[1],""))

ggsave("../Figures/Final_figs/Figure_3_pattern.pdf",Fig_2,width = 9,height = 10)


# Figure 4: Multistability trait ----


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
                                     Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))]
                                     ))
      }
    }
  }
}

trait1_for_bifu=c(unique(d2$Psi1)[15],unique(d2$Psi1)[49],unique(d2$Psi1)[49])
trait2_for_bifu=c(unique(d2$Psi2)[13],unique(d2$Psi2)[1],unique(d2$Psi2)[49])

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),color="black",alpha=.3,breaks = c(.2,.3,.4,.5))+
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
                    melt(.,measure.vars=c("Stress_tolerant","Competitive"))%>%
                    mutate(., variable=recode_factor(variable,"Stress_tolerant"="Sp. 1","Competitive"="Sp. 2")))+
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
                    melt(.,measure.vars=c("Stress_tolerant","Competitive"))%>%
                    mutate(., variable=recode_factor(variable,"Stress_tolerant"="Sp. 1","Competitive"="Sp. 2")))+
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

ggsave(filename = "../Figures/Final_figs/Figure_4.pdf",plot = p_tot,width = 7,height = 7)







# Figure 5: N-species spatially explicit ----

list_nsp_sim=list.files("../Table/N_species/CA/",pattern = "Sim")[-grep(pattern = "0.4",x = list.files("../Table/N_species/CA/",pattern = "Sim"))]
list_nsp_landscape=list.files("../Table/N_species/CA/",pattern = "Landscape")[-grep(pattern = "0.4",x = list.files("../Table/N_species/CA/",pattern = "Landscape"))]
list_nsp_sim=list_nsp_sim[-grep(pattern = "0.35",x = list_nsp_sim)]
list_nsp_landscape=list_nsp_landscape[-grep(pattern = "0.35",x = list_nsp_landscape)]
for (k in 1:length(list_nsp_sim)){
  
  
  d_sim=read.table(paste0("../Table/N_species/CA/",list_nsp_sim[k]),sep=",")
  d_landscape=read.table(paste0("../Table/N_species/CA/",list_nsp_landscape[k]),sep=",")
  
  assign(paste0("p_1_",k),
         plot_dynamics(d_sim)+labs(y="Cover",x=TeX(r'(Time, x $10^5)'))+
           scale_x_continuous(breaks = seq(0,50000,by=10000),
                              labels = c('0','1',"2",'3','4','5'))
  )
  
  assign(paste0("p_2_",k),
         Plot_landscape(d_landscape,Nsp=5)+theme(legend.position = "none")
  )
}

#low competition
p_top_sim=ggarrange(p_1_1+ggtitle("Stress (S) = 0"),p_1_2+ggtitle("Stress (S) = 0.7"),
                    ncol=2,common.legend = T,legend = "bottom")
p_top_landscape=ggarrange(ggplot()+theme_void(),p_2_1+ggtitle(""),ggplot()+theme_void(),p_2_2+ggtitle(""),ggplot()+theme_void(),
                          ncol=5,common.legend = T,legend = "none",widths = c(.05,1,.1,1,.05))
p_top=ggarrange(p_top_sim,p_top_landscape,nrow=2,heights = c(1.5,1))

#higher competition
p_bottom_sim=ggarrange(p_1_3+ggtitle("Stress (S) = 0"),p_1_4+ggtitle("Stress (S) = 0.7"),
                       ncol=2,common.legend = T,legend = "bottom")
p_bottom_landscape=ggarrange(ggplot()+theme_void(),p_2_3+ggtitle(""),ggplot()+theme_void(),p_2_4+ggtitle(""),ggplot()+theme_void(),
                             ncol=5,common.legend = T,legend = "none",widths = c(.05,1,.1,1,.05))
p_bottom=ggarrange(p_bottom_sim,p_bottom_landscape,nrow=2,heights = c(1.5,1))


p_tot=ggarrange(p_top,p_bottom,nrow=2,labels = LETTERS[1:2])

ggsave("../Figures/Final_figs/Figure_5.pdf",p_tot,width=6,height=8)



# Figure 6: Alternative community states ----
d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

p1=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.35)))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.1,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())

p2=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.35)))+
  geom_point(aes(x=Stress,y=Rho_plus,color=Psi_normalized,fill=Psi_normalized),size=.3,shape=21)+
  the_theme+labs(y="Vegetation cover",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())


d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")


d_richness_mean=d_richness%>%
  group_by(., Nsp,Branch,Competition)%>%
  summarise(., mean_richness=mean(Richness),.groups = "keep")

p3=ggplot(d_richness)+
  geom_bar(aes(x=Competition,fill=as.factor(Richness),color=as.factor(Richness),y=Richness),
           color="transparent",position = "fill",stat="identity")+
  geom_point(data=d_richness_mean,aes(x=Competition,y=mean_richness/7),
             shape=1,size=3,color="white")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),
       y="Fraction of initial conditions",fill="Number of species \n with positive cover")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  scale_y_continuous(sec.axis=sec_axis(~.*7, name=TeX(r'(Mean community \n trait, \ $\bar \psi)')))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025),labels = c("0.225","0.25",'0.275','0.3','0.325','0.35'))


Fig_6=ggarrange(ggarrange(p2,p1,common.legend = T,legend = "bottom",ncol=2,labels = letters[1:2]),
                ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),ncol=3,widths = c(.05,1.5,.05),labels = c("",letters[3],"")),nrow=2,heights = c(1,1.2))

ggsave("../Figures/Final_figs/Figure_6.pdf",Fig_6,width = 9,height = 8)









# old ----

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

# Figure 3: Trait variation and multistability 

c_inter_seq=c(.1,.2, .3)
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



d[,1:2][d[,1:2] < 10^-3] = 0


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
                         ggtitle(TeX(r'(Stress-tolerant sp., $\psi_1 = 1$)')),
                       p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank(),axis.title.y.right = element_text(size = 11))+
                         scale_y_continuous(sec.axis = sec_axis(trans = ~ .x ,
                                                                name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
                         ggtitle(TeX(r'(Competitive sp., $\psi_2 = 0$)')),
                       common.legend = T,legend = "bottom",nrow=2,labels = letters[1:2])

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
                         ggtitle(TeX(r'(Stress-tolerant sp., $\psi_1 = 1$)')),
                       p2_pattern+guides(pattern=F)+theme(strip.background.x = element_blank(),strip.text.x = element_blank(),axis.title.y.right = element_text(size = 11))+
                         scale_y_continuous(sec.axis = sec_axis(trans = ~ .x ,
                                                                name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
                         ggtitle(TeX(r'(Competitive sp., $\psi_2 = 0$)')),
                       common.legend = T,legend = "bottom",nrow=2,labels = letters[1:2])

Figure_3_tot=ggarrange(Figure_3_top,Figure_3_bottom,nrow=2,heights = c(3,1))

ggsave("../Figures/Final_figs/Figure_3_pattern.pdf",Figure_3_tot,width = 9,height = 9)


# Figure 4: Position of tipping points 

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







# Figure 3& 4 bis proposition 

c_inter_seq=c(.1,.2, .3)
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



d[,1:2][d[,1:2] < 10^-3] = 0


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
                                     "Species 2/Coexistence"="Coexistence/Species 2",
                                     "Desert/Species 2"="Species 2/Desert",
                                     "Desert/Coexistence"="Coexistence/Desert",
                                     "Stress_tolerant"="Stress-tolerant",
                                     "Desert/Stress_tolerant"="Stress-tolerant/Desert",
                                     "Stress_tolerant/Coexistence"="Coexistence/Stress-tolerant",
                                     "Stress_tolerant/Species 2"="Species 2/Stress-tolerant",
                                     "Species 2/Coexistence"="Coexistence/Species 2",
                                     "Coexistence/Species 2"="Coexistence/Species 2",
                                     "Species 2" = "Species 2"))
    
    color_multistability=c("Coexistence" = "#D8CC7B",
                           "Species 2" = "#ACD87B",
                           "Coexistence/Species 2" = "#DDEFCA",
                           "Stress-tolerant" = "#7BD8D3",
                           "Stress-tolerant/Desert" ="#0F8E87",
                           "Coexistence/Stress-tolerant"="#9BBBB9",
                           "Coexistence/Desert"="#C19E5E",
                           "Desert"=  "#696969",
                           "Species 2/Stress-tolerant" = "#C998CE")
    #State multistability
    
    p1=ggplot(d_state) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = "gray20",
                        pattern_colour  = "gray20")+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait species 2, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=color_multistability)+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      scale_pattern_manual(values=rev(c("none" ,"stripe")))+
      guides(pattern=F)
    
    
    
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
                                     "Species 2/Coexistence"="Coexistence/Species 2",
                                     "Desert/Species 2"="Species 2/Desert",
                                     "Desert/Coexistence"="Coexistence/Desert",
                                     "Desert/Competitive"="Competitive/Desert",
                                     "Competitive/Coexistence" = "Coexistence/Competitive"))
    
    
    color_multistability=c("Coexistence" = "#D8CC7B",
                           "Species 2" = "#7BD8D3",
                           "Coexistence/Species 2" = "#9BBBB9",
                           "Species 2/Desert" ="#0F8E87",
                           "Coexistence/Desert"="#C19E5E",
                           "Desert"=  "#696969")
    
    p2=ggplot(d_state%>%
                mutate(Psi2=round(Psi2,5),Stress=round(Stress,5))) +
      geom_tile_pattern(aes(x=Stress,y=as.numeric(Psi2),fill=all_state,pattern=Multistability),
                        colour          = 'transparent',
                        pattern_spacing = 0.05, 
                        pattern_density = 0.01, 
                        pattern_fill    = "gray20",
                        pattern_colour  = "gray20")+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait species 2, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=color_multistability)+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      scale_pattern_manual(values=rev(c("none" ,"stripe")))+
      guides(pattern=F)
    
    
    
  }
}


c_inter_seq=c(.1, .2, .3)
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
    
    p3=ggplot(d_tipping%>%filter(., Psi2 != 1))+
      geom_smooth(aes(x=Psi2,y=Tipping_C,color=Branch,group=interaction(Competition,Branch)),se = F)+
      the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
      labs(x=TeX(r'(Trait species 2, \ $\psi_2)'),y="Threshold stress",color="")+
      scale_color_manual(values=c("black","blue"))
    
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
    
    p4=ggplot(d_tipping)+
      geom_smooth(aes(x=Psi2,y=Tipping_ST,color=Branch,group=interaction(Competition,Branch)),se = F)+
      the_theme+facet_grid(.~Competition,labeller = label_bquote(cols=alpha[e]==.(Competition)))+
      labs(x=TeX(r'(Trait species 2, \ $\psi_2)'),y="Threshold stress",color="")+
      scale_color_manual(values=c("black","blue"))
  }
}


p_fig3=ggarrange(p1+scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                                           name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
                   ggtitle(TeX(r'(Stress-tolerant sp., $\psi_1 = 1$)'))+theme(legend.text = element_text(size=11)),
                 ggarrange(p3+theme(strip.background.x = element_blank(),
                                    strip.text.x = element_blank(),legend.text = element_text(size=12)),
                           ggplot()+theme_void(),widths = c(17,1)), nrow = 2,heights = c(1.5,1),
                 labels = letters[1:2],align = 'hv')


p_fig4=ggarrange(p2+scale_y_continuous(sec.axis = sec_axis(trans = ~ .x ,
                                                           name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
                   ggtitle(TeX(r'(Competitive sp., $\psi_1 = 0$)'))+theme(legend.text = element_text(size=11)),
                 ggarrange(p4+theme(strip.background.x = element_blank(),
                                    strip.text.x = element_blank(),legend.text = element_text(size=12)),
                           ggplot()+theme_void(),widths = c(17,1)), nrow = 2,heights = c(1.5,1),
                 labels = letters[1:2],align = 'hv')


ggsave("../Figures/Final_figs/Fig3_alter.pdf",p_fig3,width = 10,height = 8)
ggsave("../Figures/Final_figs/Fig4_alter.pdf",p_fig4,width = 10,height = 8)



# Figure 6: Species diversity 

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




