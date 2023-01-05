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
         ggplot(d2%>%melt(., measure.vars=c("Rho_plus"))%>%
                  filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,.4)[i]))+
           
           geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
           the_theme+
           scale_color_manual(values="black")+
           labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+
           guides(color="none")+
           scale_x_continuous(breaks = c(0,.25,.5,.75,1))    )
  
  
  assign(paste0("p1_",i),
         ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                  filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,.4)[i])%>%
                  mutate(.,variable=recode_factor(variable,"Stress_tolerant"="Stress-tolerant  ")))+
           
           geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
           the_theme+theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x = element_blank(),axis.line.x = element_blank())+
           scale_color_manual(values=c(as.character("#858BE0"),as.character(color_rho[c(2)])))+
           labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+
           scale_x_continuous(breaks = c(0,.25,.5,.75,1))+
           theme(legend.text = element_text(size=12))+
           guides(color = guide_legend(override.aes = list(size = 1.5)),
                  linetype = guide_legend(override.aes = list(size = 1))))

}



list_land=list.files("../Table/2_species/CA/Illustration/",pattern = "Landscape")
list_dim=list.files("../Table/2_species/CA/Illustration/",pattern = "Dynamics")
color_rho2 = c("Fertile" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress-tolerant" = "#858BE0")

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

p4_2=p4_2+scale_fill_gradientn(colours = c("white","white","#858BE0"))

p_top=ggarrange(p1_1+ylab("Species cover")+ggtitle(TeX("$\\alpha_e=0$")),
                p1_2+geom_vline(xintercept = 0.2)+ggtitle(TeX("$\\alpha_e=0.4$"))+geom_text(aes(x=.24,y=.8,label="c"),state="unique"),
                ncol=2,labels = c(letters[1],""),common.legend = T,legend = "bottom")

p_middle = ggarrange(p2_1+ylab("Vegetation cover"),p2_2,ncol=2,labels = c(letters[2],""),common.legend = T,legend = "none",hjust=-.5,vjust=0)
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

d_state=d_state[order(d_state$alpha_0,d_state$Scena,d_state$Delta,d_state$Stress),]

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
            filter(., Scena=="Local facilitation")%>%
            mutate(., Delta=as.character(Delta))) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',
                    color          = 'transparent',
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
                             "Stress-tolerant" = "#BABEFF",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(.~Delta,labeller = label_bquote(cols= delta==.(Delta) ))+
  the_theme+theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=11))+
  theme(legend.text = element_text(size = 9))+
  scale_pattern_manual(values=rev(c("transparent" ,"stripe")))


#invasion and extinction
d_invade=read.table("../Table/2_species/PA/Threshold_invasion.csv",sep=";")
d_extinction=read.table("../Table/2_species/PA/Threshold_extinction.csv",sep=";")
colnames(d_invade)=colnames(d_extinction)
d_threshold=rbind(d_invade%>%add_column(., Branch="Restoration"),
                  d_extinction%>%add_column(., Branch="Desertification"))

d_ribbon=tibble(x=d_invade$delta,ymax=d_invade$Thresh_extinction,
                ymin=d_extinction$Thresh_extinction,Scena=d_invade$Scena)

p2=ggplot(d_threshold%>%filter(., Scena==2))+
  geom_point(aes(x=delta,y=Thresh_extinction,shape=Branch),size=3)+
  geom_line(aes(x=delta,y=Thresh_extinction,group=Branch),lwd=.5)+
  #geom_ribbon(data=d_ribbon%>%filter(., Scena==2),aes(x=x,ymin=ymin,ymax=ymax),fill="#0F8E87", alpha=0.15)+
  labs(x=TeX('$\\delta$'),y="Threshold of stress",color="",shape="")+
  the_theme+
  scale_shape_manual(values=c(1,8))+
  guides(shape=guide_legend(ncol=1))+
  ggtitle(TeX(r'($\alpha_e = 0.2)'))





# clustering of stress-tolerant species
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


d_clust_mean=d_clustering%>%
  filter(., alpha_0==.2, S %in% c(.73),Scena=="Local facilitation",
         delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
  group_by(., cintra,alpha_0,S,delta,Scena)%>%
  summarise(.,.groups ="keep",c11=mean(c11) )


p30=ggplot(d_clust_mean)+
  geom_point(aes(x=delta,y=c11),size=3,shape=1)+
  geom_line(aes(x=delta,y=c11,group=interaction(S)),lwd=.5)+
  the_theme+labs(x=TeX('$\\delta$'),color="",y="Stress-tolerant clustering")+
  geom_hline(yintercept = 1)+
  ggtitle(TeX(r'($\S = 0.73, \alpha_e = 0.2)'))


list_landscape=list.files("../Table/2_species/CA/",pattern = "clustering")
for (i in 1:length(list_landscape)){
  
  landscape=read.table(paste0("../Table/2_species/CA/",list_landscape[i]),sep=",")
  assign(paste0("p_land_",i),Plot_landscape(landscape,Nsp = 2)+
          scale_fill_gradientn(colours = c("white","white","#858BE0"))+
           theme(legend.position = "none"))
  
}

p3=p30 + annotation_custom(grob=ggplotGrob(p_land_1),
                       ymin = 1.5, ymax=1.9, xmin=.2, xmax=.7)+
  annotation_custom(grob=ggplotGrob(p_land_2),
                    ymin = 1.2, ymax=1.6, xmin=.6, xmax=1.1)+
  geom_segment(aes(x = d_clust_mean$delta[2], y = d_clust_mean$c11[2], xend = d_clust_mean$delta[2]+.125, yend = d_clust_mean$c11[2]),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_segment(aes(x = d_clust_mean$delta[11], y = d_clust_mean$c11[11], 
                   xend = d_clust_mean$delta[11], yend = d_clust_mean$c11[11]+.15),
               arrow = arrow(length = unit(0.2, "cm")),lwd=.1)



# density of pairs

d_pairs=d_clustering%>%
  filter(., alpha_0==.2, S %in% c(0),Scena=="Local facilitation",
         delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
  melt(., measure.vars=c("Rho_20","Rho_10"))%>%
  mutate(.,variable=recode_factor(variable,"Rho_10"="Pair Stress-tol. sp./Fertile",
                                  "Rho_20"="Pair Competitive sp./Fertile"))
p40=ggplot(d_pairs)+
  geom_point(aes(x=delta,y=value,color=variable),size=3,shape=1)+
  geom_line(aes(x=delta,y=value,color=variable,group=interaction(variable)),lwd=.5)+
  ggtitle(TeX(r'($\S = 0, \alpha_e = 0.2)'))+
  the_theme+
  labs(x=TeX('$\\delta$'),y=expression(paste("Pair cover")),color="")+
  scale_color_manual(values=c("#858BE0",as.character(color_Nsp(5)[c(2)])))+
  theme(legend.box = "vertical")+
  guides(color=guide_legend(ncol=1))



list_landscape=list.files("../Table/2_species/CA/",pattern = "Pairs")
for (i in 1:length(list_landscape)){
  landscape=read.table(paste0("../Table/2_species/CA/",list_landscape[i]),sep=",")
  assign(paste0("p_land_",i),Plot_landscape(landscape,Nsp = 2)+
           theme(legend.position = "none"))
}

p4=p40 + annotation_custom(grob=ggplotGrob(p_land_1),
                           ymin = .052, ymax=.092, xmin=.25, xmax=.7)+
  annotation_custom(grob=ggplotGrob(p_land_2),
                    ymin = .052, ymax=.092, xmin=.65, xmax=1.1)+
  
  geom_segment(aes(x = d_pairs$delta[2], y = 0.065, 
                   xend = d_pairs$delta[2], yend =.072),lwd=.1)+
  geom_segment(aes(x = d_pairs$delta[2], y = 0.08, 
                   xend = d_pairs$delta[2], yend =.072),lwd=.1)+
  
  geom_segment(aes(x = d_pairs$delta[2], y = .072, 
                   xend = d_pairs$delta[2]+.175, yend = .072),
               arrow = arrow(length = unit(0.2, "cm")),lwd=.1)+
  
  geom_segment(aes(x = d_pairs$delta[11], y = 0.04, 
                   xend = d_pairs$delta[11], yend = 0.0525),
               arrow = arrow(length = unit(0.2, "cm")),lwd=.1)






p_bottom=ggarrange(p4,
                   ggarrange(p3,ggplot()+theme_void(),nrow=2,heights = c(1,.25)),
                   p2,ncol=3,labels = letters[2:4],align = "h")
Fig_3=ggarrange(p1+guides(pattern=F),p_bottom,nrow=2,heights = c(3,2.3),labels = c(letters[1],""))

ggsave("../Figures/Final_figs/Figure_3.pdf",Fig_3,width = 9,height = 9)


# Figure 4: Multistability trait ----


stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",round(stress,4),"_delta_0.1_facilitation_0.9_scalefacilitation_local.csv"),sep=";"),
           read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",round(stress,4),"_delta_0.1_facilitation_0.9_scalefacilitation_local.csv"),sep=";"))
  d=rbind(d,d2)
  
  
  
}



d2=d
d2$Stress=round(d2$Stress,4)
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
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),
               color="black",alpha=.3,breaks = c(.2,.3,.4,.5))+
  geom_point(data=tibble(x=trait1_for_bifu,y=trait2_for_bifu),aes(x=x,y=y),shape=19)+
  geom_text(data=tibble(x=c(trait1_for_bifu[1],trait1_for_bifu[2]+.05,trait1_for_bifu[3]+.05),
                        y=c(trait2_for_bifu[1]+.13,trait2_for_bifu[2]+.05,trait2_for_bifu[3]+0.05),
                        txt=c("(ii)","(i)","(iii)")),aes(x=x,y=y,label=txt),size=4.5)+
  geom_segment(aes(x = 0.05, y = -0.02, xend = .95, yend = -0.02),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = .5, y = -0.06,label = "Sensitivity to competition"),stat = "unique")+
  # geom_segment(aes(x = 1.02, y = 0, xend = 1.02, yend = .95),arrow = arrow(length = unit(0.3, "cm")))+
  # geom_text(aes(x = 1.05, y = .45,label = "Trait similarity"),stat = "unique",angle=270)+
  geom_segment(aes(x = 1, y = 0, xend = .52, yend = .52),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = .875, y = .3,label = "Trait similarity"),stat = "unique")+
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
             scale_color_manual(values=c(as.character(color_Nsp(5)[c(2)]),"#858BE0"))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank(),
                   legend.position = "none")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             guides(color = guide_legend(override.aes = list(size = 1.5)))+
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
             scale_color_manual(values=c(as.character(color_Nsp(5)[c(2)]),"#858BE0"))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             theme(legend.text = element_text(size=11))+
             guides(color = guide_legend(override.aes = list(size = 1.5)))+
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

d_pa_equal=read.table("../Table/N_species/PA/PA_equal_ini.csv",sep=",")
d_pa_equal=d_pa_equal[,c(1:15,(ncol(d_pa_equal)-2):ncol(d_pa_equal))]
colnames(d_pa_equal)=c(paste0("Sp_",1:15),"alpha_0", "Branches", "Stress")
d_melt=d_pa_equal%>%
  melt(., measure.vars=paste0("Sp_",1:15))

d_melt$variable=as.character(d_melt$variable)

d_melt$Trait=sapply(1:nrow(d_melt),function(x){
  seq(0,1,length.out=15)[as.numeric(strsplit(d_melt$variable[x],split = "_")[[1]][2])]
})
  
p1=ggplot(d_melt%>%
            mutate(., Branches=recode_factor(Branches,"1"="Degradation","2"="Restoration"))) +
  geom_line(aes(x=Stress,y=value,color=Trait,group=interaction(Trait,alpha_0,Branches),linetype=as.factor(Branches)))+
  facet_wrap(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)),scales = "free")+
  the_theme+
  labs(x="Stress, (S)", y="Species cover",linetype="",color="")+
  scale_color_gradientn(colours = rev(color_Nsp(15)))+
  theme(strip.text.x = element_text(size=11))



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
                              labels = c('0','1',"2",'3','4','5'))+theme(legend.position = "none")
  )
  
  assign(paste0("p_2_",k),
         Plot_landscape(d_landscape,Nsp=5)+theme(legend.position = "none")
  )
}

#low competition
p_left_sim=ggarrange(p_1_1+ggtitle(TeX("$\\alpha_e = 0.1$")),p_1_3+ggtitle(TeX("$\\alpha_e = 0.25$")),
                    ncol=2,common.legend = T,legend = "none")
p_left_landscape=ggarrange(ggplot()+theme_void(),p_2_1+ggtitle(""),ggplot()+theme_void(),p_2_3+ggtitle(""),ggplot()+theme_void(),
                          ncol=5,common.legend = T,legend = "none",widths = c(.05,1,.1,1,.05))

#higher competition
p_right_sim=ggarrange(p_1_2+ggtitle(TeX("$\\alpha_e = 0.1$")),p_1_4+ggtitle(TeX("$\\alpha_e = 0.25$")),
                       ncol=2,legend="none")
p_right_landscape=ggarrange(ggplot()+theme_void(),p_2_2+ggtitle(""),ggplot()+theme_void(),p_2_4+ggtitle(""),ggplot()+theme_void(),
                             ncol=5,widths = c(.05,1,.1,1,.05),legend="none")
p_top=ggarrange(p_left_sim,p_right_sim,ncol=2,labels=letters[2:3])

p_bottom=ggarrange(p_left_landscape,p_right_landscape,ncol=2,legend="none")


p_spatial=ggarrange(p_top,p_bottom,nrow=2,heights = c(1.5,1),legend="none")
p_tot=ggarrange(ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.05,.7,.05)),
                p_spatial,nrow=2,heights = c(1.5,2),labels = c(letters[1],""))

ggsave("../Figures/Final_figs/Figure_5.pdf",p_tot,width=9,height=7)



# Figure 6: Alternative community states ----


d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")



p1=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.35),Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.1,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())

p2=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.35),Branch==1))+
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

p3=ggplot(d_richness%>%filter(., Branch=="Degradation"))+
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
  scale_y_continuous(sec.axis=sec_axis(~.*7, name="Mean number of species"))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025),labels = c("0.225","0.25",'0.275','0.3','0.325','0.35'))


Fig_6=ggarrange(ggarrange(p2,p1,common.legend = T,legend = "bottom",ncol=2,labels = letters[1:2]),
                ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),ncol=3,widths = c(.05,1.5,.05),labels = c("",letters[3],"")),nrow=2,heights = c(1,1.2))

ggsave("../Figures/Final_figs/Figure_6.pdf",Fig_6,width = 9,height = 8)












# Alternative Fig 6 ----
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


Fig_to_plot=tibble(Nsp=c(5,15,15,15,25),Number_plot=c(3,169,214,737,397))


for (i in 1:nrow(Fig_to_plot)){
  
  Nsp=Fig_to_plot$Nsp[i]
  list_csv=list.files(paste0('../Table/N_species/MF/',Nsp,'_sp/'))
  d2=read.table(paste0("../Table/N_species/MF/",Nsp,"_sp/",list_csv[Fig_to_plot$Number_plot[i]+1]),sep=",")
  colnames(d2)=c(paste0(1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
  
  d2[d2<10^(-4)]=0
  
  
  
  assign(paste0("p_",i),
         ggplot(d2%>%
                  melt(., measure.vars=paste0(1:Nsp))%>%
                  mutate(., Trait=sapply(1:nrow(.),function(x){
                    return(seq(1,0,length.out=Nsp)[as.numeric(.$variable[x])])
                  })))+
           geom_line(aes(x=Stress,y=value,color=Trait,group=variable),size=.8)+
           the_theme+labs(x="Stress, S",y="Species cover",color=TeX("$\\psi$   "))+
           scale_color_gradientn(colors=color_Nsp(Nsp))
  )
}

p=ggarrange(p_1+ggtitle(TeX(r'(N species = 5, \ $\alpha_e = 0.3)'))+theme(plot.title = element_text(size=10)),
            p_2+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.25)'))+
              theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank())+theme(plot.title = element_text(size=10)),
            p_3+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.325)'))+theme(plot.title = element_text(size=10)),
            p_5+ggtitle(TeX(r'(N species = 25, \ $\alpha_e = 0.25)'))+
              theme(axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank())+theme(plot.title = element_text(size=10)),
            nrow=1,ncol=4,common.legend = T,legend="bottom")


Fig_6=ggarrange(ggarrange(p2,p1,common.legend = T,legend = "bottom",ncol=2,labels = letters[1:2]),
                ggarrange(p,labels = letters[3]),
                ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),ncol=3,widths = c(.05,1.5,.05),labels = c("",letters[4],"")),
                nrow=3,heights = c(1.1,.75,1.3))

ggsave("../Figures/Final_figs/Possibility_fig_6.pdf",Fig_6,width = 9,height = 10)


