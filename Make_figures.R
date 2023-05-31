source('./Dryland_shift_functions.R')
library(grid)



#---------------------   Main figures   --------------------------

#*****************************************************************

## >> Figure 2: Scale facilitation, dispersal, multistability ----
"
Code to replicate figure 3. To generate the data needed for the figure, 
run chunks Step 2.1 to Step 2.3 of the Dryland_shift_main.R file.
In addition, the landscapes needed to illustrate the change in metrics along the 
dispersal gradient are generated using the region 2 and 4 of the julia file
Dryland_shift_Nspecies_main.jl
"

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
                             #"Competitive/Stress-tolerant"="#C998CE",
                             "Competitive" = "#ACD87B", 
                             #"Competitive/Coexistence" = "#DDEFCA",
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
  scale_color_manual(values=c("#858BE0",as.character(color_Nsp(5)[c(2)])),
                     labels=c(expression(paste("Pair Stress???tol. sp./Fertile ", rho["0,+"[1]])),
                              expression(paste("Pair Competitive sp./Fertile ", rho["0,+"[2]]))))+
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
Fig_2=ggarrange(p1+guides(pattern=F),p_bottom,nrow=2,heights = c(3,2.3),labels = c(letters[1],""))

ggsave("../Figures/Figure_2.pdf",Fig_3,width = 9,height = 9)


## >> Figure 3: N-species spatially explicit ----

"
Code to replicate figure 5: scaling_up. To generate the data needed for the figure, 
run regions 6 and 8 of the julia script Dryland_shift_Nspecies_main.R file
"

d_pa_equal=read.table("../Table/N_species/PA/PA_equal_ini.csv",sep=",")
d_pa_equal=d_pa_equal[,c(1:5,(ncol(d_pa_equal)-2):ncol(d_pa_equal))]
colnames(d_pa_equal)=c(paste0("Sp_",1:5),"alpha_0", "Branches", "Stress")
p1=d_pa_equal%>%
  melt(., measure.vars=paste0("Sp_",1:5))%>%
  mutate(., variable=as.character(variable))%>%
  add_column(., Trait=sapply(1:nrow(.),function(x){
    seq(1,0,length.out=5)[as.numeric(strsplit(.$variable[x],split = "_")[[1]][2])]
  }))%>%
  mutate(., Branches=recode_factor(Branches,"1"="Degradation","2"="Restoration"))%>%
  ggplot(.)+
  geom_line(aes(x=Stress,y=value,color=Trait,group=interaction(Trait,alpha_0,Branches),linetype=as.factor(Branches)))+
  facet_wrap(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)),scales = "free")+
  the_theme+
  labs(x="Stress, (S)", y="Species cover",linetype="",color="")+
  scale_color_gradientn(colours = (color_Nsp(5)))+
  theme(strip.text.x = element_text(size=11))



list_nsp_sim=list.files("../Table/N_species/CA/",pattern = "Sim")
list_nsp_landscape=list.files("../Table/N_species/CA/",pattern = "Landscape")
for (k in 1:length(list_nsp_sim)){
  
  
  d_sim=read.table(paste0("../Table/N_species/CA/",list_nsp_sim[k]),sep=",")
  d_landscape=read.table(paste0("../Table/N_species/CA/",list_nsp_landscape[k]),sep=",")
  
  assign(paste0("p_1_",k),
         plot_dynamics(d_sim)+labs(y="Cover",x=TeX(r'(Time, x $10^4)'))+
           scale_x_continuous(breaks = seq(0,50000,by=10000),
                              labels = c('0','1',"2",'3','4','5'))+theme(legend.position = "none")
  )
  
  assign(paste0("p_2_",k),
         Plot_landscape(d_landscape,Nsp=5)+theme(legend.position = "none")
  )
}

#low competition
p_left_sim=ggarrange(p_1_1+ggtitle(TeX("$\\alpha_e = 0.15, S = 0")),p_1_3+ggtitle(TeX("$\\alpha_e = 0.3, S = 0$")),
                     nrow=2,common.legend = T,legend = "none")
p_left_landscape=ggarrange(ggplot()+theme_void(),p_2_1+ggtitle(""),ggplot()+theme_void(),p_2_3+ggtitle(""),ggplot()+theme_void(),
                           nrow=5,heights = c(.15,1,.1,1,.25),legend="none")

p_left=ggarrange(p_left_sim,p_left_landscape,ncol=2,labels=c(letters[2],""),widths = c(1.3,1))


#higher competition
p_right_sim=ggarrange(p_1_2+ggtitle(TeX("$\\alpha_e = 0.15, S = 0.6$")),p_1_4+ggtitle(TeX("$\\alpha_e = 0.3, S = 0.6$")),
                      nrow=2,legend="none")
p_right_landscape=ggarrange(ggplot()+theme_void(),
                            p_2_2,
                            ggplot()+theme_void(),
                            p_2_4,
                            ggplot()+theme_void(),
                            nrow=5,heights = c(.15,1,.1,1,.25),legend="none")
p_right=ggarrange(p_right_sim,p_right_landscape,ncol=2,labels=c(letters[3],""),widths = c(1.3,1))

p_bottom=ggarrange(p_left,p_right,ncol=2,legend="none")


p_tot=ggarrange(ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.05,.7,.05)),
                p_bottom,nrow=2,heights = c(1.5,2),labels = c(letters[1],""))

ggsave("../Figures/Figure_3.pdf",p_tot,width=9,height=7)





## >> Figure 4: Alternative community states ----

"
Code to replicate figure 4. To generate the data needed for the figure, 
run chunk Step 3.1) of the Dryland_shift_main.R file
"

d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")

p1=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.30),Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.1,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="# species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())

p2=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.30),Branch==1))+
  geom_point(aes(x=Stress,y=Rho_plus,color=Psi_normalized,fill=Psi_normalized),size=.3,shape=21)+
  the_theme+labs(y="Vegetation cover",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="# species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())


d=read.table("../Table/N_species/PA/post_proc_sim1.csv",sep=";")

p3=ggplot(d%>%filter(., Competition>.2))+
  geom_smooth(aes(x=Stress,y=N_ASS,color=as.factor(Competition),group=Competition),
              se = F)+
  the_theme+
  facet_grid(Nsp~.,labeller = label_bquote(rows= "# species" ==.(Nsp)))+
  labs(x="Stress (S)",y="# of alternative states",color=TeX("$\\alpha_e \ \ $"))+
  scale_color_viridis_d()+
  theme(strip.text.y = element_text(size=13),panel.background = element_blank(),
        strip.background.y = element_blank(),legend.title = element_text(size=14))


Fig_4=ggarrange(ggarrange(p2,p1,common.legend = T,legend = "bottom",nrow=2,labels = letters[1:2]),
                ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),nrow=3,heights = c(.05,1.5,.05),labels = c("",letters[3],"")),ncol=2,widths =  c(1.2,1))

ggsave("../Figures/Figure_4.pdf",Fig_4,width = 9,height = 8)




## >> Figure 5: Type of multistability ----
"
Code to replicate figure 5. To generate the data needed for the figure, 
run chunk Step 3.1) of the Dryland_shift_main.R file
"

d=read.table("../Table/N_species/PA/post_proc_sim1.csv",sep=";")
type_bistab=read.table("../Table/N_species/PA/post_proc_sim2.csv",sep=";")
type_bistab$Type[type_bistab$Stress<.2]="Cliques"

d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")


#2D param space with mean # of sp coexisting


assign(paste0("p1_1"),
         type_bistab%>%filter(., Nsp==25,Competition>.15)%>%
         mutate(., Competition=as.factor(Competition))%>%
         ggplot(.)+
           geom_tile(aes(x=Stress,y=Competition,fill=Type))+
           the_theme+
           labs(fill="",x="Stress (S)",y=TeX(r'(Interspecific competition, \ $\alpha_e)'))+
           scale_y_discrete(breaks = unique(type_bistab$Competition),labels = paste0(unique(type_bistab$Competition)))+
           theme(strip.background.x = element_blank())+
           scale_fill_manual(values=c("Degraded"="#000000","Env"="#0F8E87","Mutual exclusion"="#C8A4D0",
                                      "No bistab"="gray50","Cliques"="#EFE8BC"),
                             labels=c("Degraded","Environmental","Mutual exclusion",
                                      "No bistability","Cliques"))+
           theme(strip.text.x = element_text(size=12)))




d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")%>%
  dplyr::group_by(.,Competition,Stress,Nsp)%>%
  dplyr::summarise(., .groups = "keep",mean_sp_div=mean(Nb_sp))%>%
  arrange(., Stress,Nsp,Competition)

d_tot$mean_sp_div[d_tot$mean_sp_div==0]=NA #for specific black color

#2D param space with types of multistability

assign(paste0("p1_2"),
         ggplot(d_tot%>%filter(., Nsp==25,Competition>.15)%>%
                  mutate(., Competition=as.factor(Competition)))+
           geom_tile(aes(x=Stress,y=Competition,fill=mean_sp_div),alpha=.8)+
           the_theme+
           labs(fill="# of coexisting species  ",x="Stress (S)",y=TeX(r'(Interspecific competition, \ $\alpha_e)'))+
           scale_y_discrete(breaks = unique(d_tot$Competition),labels = paste0(unique(d_tot$Competition)))+
           theme(strip.background.x = element_blank())+
           scale_fill_viridis_c(na.value = "black")+
           theme(strip.text.x = element_text(size=12)))
         


Fig_5=ggarrange(
  ggarrange(p1_2+theme(legend.position = "none"),
            p1_1+theme(legend.position = "none"),ncol=2,labels = letters[c(1,2)],align = "hv"),
  ggarrange(get_legend(p1_2+theme(legend.key.size = unit(.5, 'cm'))),
            get_legend(p1_1+guides(fill=guide_legend(nrow=2))+
                         theme(legend.spacing.x = unit(.5, 'cm'))),ncol=2),
  nrow=2,heights = c(1,.15)
  )

ggsave("../Figures/Figure_5.pdf",Fig_5,width = 9,height = 4)

#*****************************************************************

## >> Figure 6: Characteristics of cliques ----
"
Code to replicate figure 6. To generate the data needed for the figure, 
run chunk Step 3.1) of the Dryland_shift_main.R file
"


d_cover_cliques=read.table("../Table/N_species/PA/post_proc_sim3.csv",sep=";")
type_bistab=read.table("../Table/N_species/PA/post_proc_sim2.csv",sep=";")

p1=ggplot(d_cover_cliques%>%
            filter(., Competition>.15)%>%
            mutate(., Competition=as.factor(Competition))%>%
            filter(., Nb_sp>1,!is.na(Psi_normalized),Type=="Cliques",Nsp==25)%>%
            group_by(., Competition)%>%
            summarise(., .groups = "keep",
                      q1=mean(Nb_sp)-sd(Nb_sp),q3=mean(Nb_sp)+sd(Nb_sp),
                      q2=mean(Nb_sp)))+
  geom_pointrange(aes(x=(Competition),y=q2,ymax=q3,ymin=q1,fill=Competition),size=.75,shape=24,color="black")+
  the_theme+
  scale_y_continuous(breaks = c(1,3,5,7),limits = c(1,7))+
  scale_fill_manual(values=colorRampPalette(c("#1A41AF","#72B2D6","#B9D7E8","#FBEFCB","#F5CB61","#D26F3C"))(7))+
  labs(x=TeX(r'(Interspecific competition, \ $\alpha_e)'),y="mean # of species \n within the clique")+
  theme(legend.position = "none")
  
  


d_cover_cliques=read.table("../Table/N_species/PA/post_proc_sim3.csv",sep=";")

d_cover_summarized=d_cover_cliques%>%
  filter(., Nb_sp !=0,Competition>.225,Nsp==25)%>%
  add_column(., N_sp_clique=ifelse(.$Nb_sp>1,">2","1"))%>%
  dplyr::group_by(., N_sp_clique,Stress,Branch)%>%
  dplyr::summarise(., q2_cover=median(Rho_plus),q3_cover=quantile(Rho_plus,.75),
                   q1_cover=quantile(Rho_plus,.25),.groups = "keep")

#fingerprint of type of multistab on vegetation cover

p2=ggplot(d_cover_summarized)+
  geom_line(aes(x=Stress,y=q2_cover,group=as.factor(N_sp_clique),color=as.factor(N_sp_clique)),lwd=1)+
  geom_ribbon(aes(x=Stress,ymin=q1_cover,ymax=q3_cover,group=as.factor(N_sp_clique),
                  fill=as.factor(N_sp_clique)),alpha=.5)+
  
  scale_color_manual(values=(c("#EFE8BC","#C8A4D0")),labels=rev(c("1 species","# species > 1")))+
  scale_fill_manual(values=(c("#EFE8BC","#C8A4D0")),labels=rev(c("1 species","# species > 1")))+
  the_theme+
  guides(fill="none",color = guide_legend(override.aes = list(size = 3)))+
  labs(x="Stress (S)",y="Community cover \n ",color="")+
  theme(legend.text = element_text(size=12))


#Distribution of trait within the cliques

d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")%>%
  filter(., Random_ini>0,Stress %in% unique(.$Stress)[8])%>% #arbitrary level of stress
  filter(., Nsp==25,Competition==.35,Nb_sp>1)

d_trait=tibble()
list_trait=rev(seq(0,1,length.out=25))
for (k in 1:nrow(d_tot)){
  for (sp in as.numeric(as.vector(unlist(strsplit(d_tot$Name_sp[k],"_"))))){ #for each species in the clique
    d_trait=rbind(d_trait,tibble(Sp=sp,Compet=d_tot$Competition[k],
                                 Random_ini=d_tot$Random_ini[k],Trait=list_trait[sp],
                                 Stress=d_tot$Stress[k]))
  }
}

p3=ggplot(d_trait)+geom_bar(aes(Trait,y=2*(..count..)/sum(..count..)))+ #each clique = 2 sp, so we multiplyy by 2
  labs(x=TeX(r'(Species trait, \ $\psi_i)'),y="Frequency")+
  the_theme


#Example of dynamics

Fig_to_plot=tibble(Nsp=c(5,5),Number_plot=c(492,498))

d=tibble()
for (i in 1:nrow(Fig_to_plot)){
  Nsp=Fig_to_plot$Nsp[i]
  list_csv=list.files(paste0('../Table/N_species/PA/',Nsp,'_sp/'))
  d2=read.table(paste0("../Table/N_species/PA/",Nsp,"_sp/",list_csv[Fig_to_plot$Number_plot[i]]),sep=",")
  colnames(d2)=c(paste0(1:Nsp),"Random_ini","Competition","Facilitation","Branch","Stress")
  d2[d2<10^(-4)]=0
  d=rbind(d,d2%>%add_column(., ID=i))
}

p4=ggplot(d%>%filter(., Branch==1)%>%
            melt(., measure.vars=paste0(1:Nsp))%>%
            mutate(., Trait=sapply(1:nrow(.),function(x){
              return(seq(1,0,length.out=Nsp)[as.numeric(.$variable[x])])
            })))+
  geom_line(aes(x=Stress,y=value,color=Trait,group=interaction(variable,ID),linetype=as.factor(ID)),size=.8)+
  the_theme+labs(x="Stress, S",y="Species cover",color=TeX("$\\psi$   "),linetype="Two initial condition")+
  scale_color_gradientn(colors=color_Nsp(Nsp))+
  scale_linetype_manual(values=c("dotdash","solid"),labels=c("",""))+
  new_scale_color()+
  geom_line(data=type_bistab%>%filter(., Nsp==5,Competition==.35)%>%
              add_column(., Height=.8),
            aes(x=Stress,y=Height,color=as.factor(Type),group=as.factor(Type)),alpha=.8,size=2)+
  scale_color_manual(values=c("Degraded"="#000000","Env"="#0F8E87","Mutual exclusion"="#C8A4D0",
                              "No bistab"="gray50","Cliques"="#EFE8BC"),
                     labels=c("Degraded","Environmental","Mutual exclusion",
                              "No bistability","Cliques"))+
  guides(color="none",linetype="none")+
  ggtitle("Two different initial conditions")+
  theme(legend.title = element_text(size=13))

Fig_6=ggarrange(
  ggarrange(p1,p3,ncol=2,labels=letters[c(1,3)],widths = c(1,.8),align = "hv"),
  ggarrange(p2+theme(legend.position = "none"),
            p4+theme(legend.position = "none"),
            ncol=2,labels = letters[c(2,4)],widths =  c(1,.8),align = "hv"),
  ggarrange(get_legend(p2),get_legend(p4),ncol=2),heights = c(1,1,.2),nrow=3)

ggsave("../Figures/Figure_6.pdf",Fig_6,width = 8,height = 7)


#---------------------   SI figures 2-species  --------------------------
## >> Figure 2: Dynamics and landscapes ----
"
Code to replicate figure SI with 2 species. To generate the data needed for the figure, 
run chunk Step 2.1) of the Dryland_shift_main.R file
"


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
           the_theme+
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

p_tot=ggarrange(p_top,p_middle,p_bottom,nrow=3,heights = c(1.35,1,1.1))

ggsave("../Figures/SI/Example_2_species.pdf",p_tot,width = 9,height = 7)






## >> Advantage of CSI compared to Rho+----
"
Code to illustrate the advantage of using the community index compared to the vegetation cover
to distinguish between different alternative stable states
"

list_f=list.files("../Table/2_species/PA/Multistability_PA/Frac_gradient/",pattern = "delta_0.1_facilitation_0.9_scalefacilitation_local")

d=tibble()
for (f in list_f){
  d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/",f),sep=";")
  d=rbind(d,d2)
}

d[,1:2][d[,1:2] < 10^-4] = 0
set.seed(123)
u=runif(2)
d$CSI = sapply(1:nrow(d),function(x){
  return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
})


d2_fil=filter(d,Psi1==1,Psi2==unique(d$Psi2)[41])
#Bifurcation diagram

p=ggplot(d2_fil%>%melt(., measure.vars=c("Rho_plus","CSI"))%>%
           mutate(Psi2=round(abs(Psi2),2)))+
  geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
  labs(x = TeX(r'(Stress, \ $S)'), y = "", color = "",linetype="") +
  the_theme +
  scale_color_manual(values = c("black",alpha("red",.5)),labels=c("Global cover","Community index")) +
  scale_linetype_manual(values=c(1,9))+
  theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))

ggsave("../Figures/SI/Using_CSI.pdf",p,width = 6,height = 4)

## >> Competition experienced ----

"
Code to illustrate the asymetry in the competition function
"

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

ggsave("../Figures/SI/Competition_experienced.pdf",
       ggarrange(
         ggarrange(p1,
                   ggarrange(p_1,p_2,p_3,nrow = 3,labels = c('i','ii','iii'),heights = c(1,1,1.4)),labels = c("a","","",""),widths = c(2,1)),
         p2,ncol=2,common.legend = T,legend = "none",labels = c("",letters[2]),widths = c(1.5,1)),
       width = 11,height = 4)



## >> Comparing CA & PA ----
"
Code to illustrate the convergence between the pair approximation and the cellular automata
To run the figure, the chunk 2 of the julia file Dryland_shift_Nspecies_main.jl must be ran
"

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
                             "Stress-tolerant" = "#BABEFF",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9))
ggsave("../Figures/SI/Comparizon_CA_PA.pdf",p,width = 7,height = 6)


## >> Multistability varying traits PA model ----
"
Code to replicate changing species composition. To generate the data needed for the figure, 
run chunk Step 2.4) of the Dryland_shift_main.R file
"


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
trait2_for_bifu=c(unique(d2$Psi2)[13],unique(d2$Psi2)[1],unique(d2$Psi2)[49]-.01)

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),
               color="white",alpha=.3,breaks = c(.2,.3,.4,.5))+
  # color="black",alpha=.3,breaks = c(.2,.3,.4,.5))+
  geom_text(data=tibble(x=c(trait1_for_bifu[1],trait1_for_bifu[2]+.05,trait1_for_bifu[3]+.05),
                        y=c(trait2_for_bifu[1]+.13,trait2_for_bifu[2]+.05,trait2_for_bifu[3]+0.05),
                        txt=c("(ii)","(i)","(iii)")),aes(x=x,y=y,label=txt),size=4.5)+
  geom_segment(aes(x = 0.05, y = -0.02, xend = .95, yend = -0.02),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = .5, y = -0.06,label = "Sensitivity to competition"),stat = "unique")+
  geom_segment(aes(x = 1, y = 0, xend = .52, yend = .52),arrow = arrow(length = unit(0.3, "cm")),color="white")+
  geom_text(aes(x = .875, y = .3,label = "Trait similarity"),stat = "unique",color="white")+
  # geom_segment(aes(x = 1, y = 0, xend = .52, yend = .52),arrow = arrow(length = unit(0.3, "cm")))+
  # geom_text(aes(x = .875, y = .3,label = "Trait similarity"),stat = "unique")+
  geom_point(data=tibble(x=trait1_for_bifu,y=trait2_for_bifu),aes(x=x,y=y),shape=19,color="white",size=2.5)+
  the_theme+
  labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
       y=TeX(r'(Species 2 trait, $\psi_2$)'),
       fill="Fraction of community bistability \n along the stress gradient")+
  # scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100))+
  scale_fill_gradientn(colors=colorRampPalette(c("#EEF0F5","#AEAFB1","#77787D","black"))(100))+
  guides(shape="none")

trait1_for_bifu=c(unique(d2$Psi1)[15],unique(d2$Psi1)[49],unique(d2$Psi1)[49])
trait2_for_bifu=c(unique(d2$Psi2)[13],unique(d2$Psi2)[1],unique(d2$Psi2)[49])

for (i in 1:3){
  
  if (i %in% c(1,3)){
    assign(paste0("p2_",i),
           ggplot(d2%>%
                    filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i])%>%
                    melt(.,measure.vars=c("Stress_tolerant","Competitive"))%>%
                    mutate(., variable=recode_factor(variable,"Stress_tolerant"="Sp. 1","Competitive"="Sp. 2")))+
             geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
             scale_color_manual(values=rev(c(as.character(color_Nsp(5)[c(2)]),"#858BE0")))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank(),
                   legend.position = "none")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             scale_linetype_manual(values=c(1,9),labels=c("High initial cover","Low initial cover"))+
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
             scale_color_manual(values=rev(c(as.character(color_Nsp(5)[c(2)]),"#858BE0")))+
             scale_linetype_manual(values=c(1,2),labels=c("High initial cover   ","Low initial cover"))+
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

ggsave(filename = "../Figures/SI/Multistab_varying_trait_PA.pdf",plot = p_tot,width = 7,height = 7)








## >> Multistability fixed traits MF model ----

"
Same as figure 3 but with the mean field model.
Please run Step 1.1 of Dryland_shift_main.R prior to generate the data needed
"

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
                    color          = 'transparent',
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
                             "Stress_tolerant" = "#BABEFF",
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
  scale_pattern_manual(values=rev(c("transparent" ,"stripe")))


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


ggsave("../Figures/SI/Multistability_fixed_traits_MF.pdf",ggarrange(p+guides(pattern=F),p2,ncol = 2,widths = c(2,1),labels=c(letters[1],"")),
       width = 10,height = 6)




















## >> Multistability varying traits MF model ----

"
Same as figure 4 but with the mean field model.
Please run Step 1.2 of Dryland_shift_main.R prior to generate the data needed
"


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
  the_theme+
  labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
       y=TeX(r'(Species 2 trait, $\psi_2$)'),
       fill="Fraction of community bistability \n along the stress gradient")+
  scale_fill_gradientn(colors=colorRampPalette(c("#EEF0F5","#AEAFB1","#77787D","black"))(100))+
  guides(shape="none")


for (i in 1:3){
  
  if (i %in% c(1,3)){
    assign(paste0("p2_",i),
           ggplot(d2%>%
                    filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i])%>%
                    melt(.,measure.vars=c("Sp1","Sp2"))%>%
                    mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
             geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
             scale_color_manual(values=rev(as.character(color_Nsp(5)[c(2,5)])))+
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
             scale_color_manual(values=rev(as.character(color_Nsp(5)[c(2,5)])))+
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

ggsave(filename = "../Figures/SI/Multistability_varying_traits_MF.pdf",p_tot,width = 7,height = 7)






## >> Multistability along competition gradient with local competition ----

"
Same as figure 2 but with local competition instead of global one
Please run Step 2.1 of Dryland_shift_main.R prior to generate the data needed
"


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
                  mutate(., Delta=recode_factor(Delta,"0.1"="0.1, local dispersal","0.9"="0.9, global dispersal"))%>%
                  mutate(., Delta=as.character(Delta))) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',   
                    color          = 'transparent',
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
                             "Stress-tolerant" = "#BABEFF",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"))+
  facet_grid(Scena~Delta,labeller = label_bquote(cols = delta==.(Delta)))+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9),strip.text.y = element_text(size=12))+
  scale_pattern_manual(values=rev(c("transparent" ,"stripe")))+
  guides(pattern="none")

ggsave("../Figures/SI/Multistability_fixed_traits_local_competition.pdf",Fig_2_SI,width = 8,height = 6)





## >> Multistability along competition gradient with global facilitation ----

"
Same as figure 2 but with global facilitation instead of local one
Please run Step 2.1 of Dryland_shift_main.R prior to generate the data needed
"


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
                             "Stress-tolerant" = "#BABEFF",
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




ggsave("../Figures/SI/Multistability_fixed_traits_global_facilitation.pdf",p1,width = 9,height = 5)


## >> Multistability varying traits with global facilitation ----

"
Same as figure 4 but with global facilitation instead of local one
Please run Step 2.2 of Dryland_shift_main.R prior to generate the data needed
"


param_space=expand.grid(Scale_facil=c("local","global"),Disp=c(.1,.9))

for (num in 1:nrow(param_space)){
  
  list_f=list.files("../Table/2_species/PA/Multistability_PA/Frac_gradient",
                    pattern =paste0("delta_",param_space$Disp[num],"_facilitation_0.9_scalefacilitation_",param_space$Scale_facil[num]) )
  d=tibble()
  for (f in list_f){
    d=rbind(d,read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/",f),sep=";")%>%
              add_column(., Scale_facil=param_space$Scale_facil[num],Disp=param_space$Disp[num]))
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
  d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1,d2$Disp,d2$Scale_facil),]
  
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
          d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,Disp=param_space$Disp[num],
                                       Scale_facil=param_space$Scale_facil[num],
                                       Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                       Frac_multi_raw=sum(d_fil$multistab),
                                       Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
                                       Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))]
          ))
        }
        
        
      }
    }
  }
  
  assign(paste0("p_",num),
         ggplot(d_final)+
           geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
           the_theme+
           labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
                y=TeX(r'(Species 2 trait, $\psi_2$)'),
                fill="Fraction of community bistability \n along the stress gradient")+
           scale_fill_gradientn(colors=colorRampPalette(c("#EEF0F5","#AEAFB1","#77787D","black"))(100))+
           guides(shape="none")+
           ggtitle(paste0(ifelse(param_space$Scale_facil[num]=="global","Global","Local")," facilitation, ",
                          ifelse(param_space$Disp[num]==.1,"local","global")," dispersal"))
         
         
  )
  
  
  
}

p=ggarrange(p_1+labs(fill="")+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()),
            p_2+labs(fill="")+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank(),
                    axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            p_3+labs(fill="")+scale_fill_gradientn(colors=colorRampPalette(c("#EEF0F5","#AEAFB1","#77787D","black"))(100),breaks=c(0,.1,.2)),
            p_4+labs(fill="")+scale_fill_gradientn(colors=colorRampPalette(c("#EEF0F5","#AEAFB1","#77787D","black"))(100),breaks=c(0,.1,.2))+
              theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            nrow=2,ncol=2,heights = c(1,1.1))

ggsave("../Figures/SI/Multistability_traits_scale_facil_disp.pdf",p,width = 7,height = 7)

## >> Mechanisms for global facilitation: clustering, pairs, thresholds ----

"
Same as figure 3a-d but with global facilitation instead of local one
Please run Step 2.2 to 2.4 of Dryland_shift_main.R prior to generate the data needed
"


#invasion and extinction
d_invade=read.table("../Table/2_species/PA/Threshold_invasion.csv",sep=";")
d_extinction=read.table("../Table/2_species/PA/Threshold_extinction.csv",sep=";")
colnames(d_invade)=colnames(d_extinction)
d_threshold=rbind(d_invade%>%add_column(., Branch="Restoration"),
                  d_extinction%>%add_column(., Branch="Desertification"))

d_ribbon=tibble(x=d_invade$delta,ymax=d_invade$Thresh_extinction,
                ymin=d_extinction$Thresh_extinction,Scena=d_invade$Scena)

p1=ggplot(d_threshold%>%filter(., Scena==1))+
  geom_point(aes(x=delta,y=Thresh_extinction,shape=Branch),size=3)+
  geom_line(aes(x=delta,y=Thresh_extinction,group=Branch),lwd=.5)+
  #geom_ribbon(data=d_ribbon%>%filter(., Scena==1),aes(x=x,ymin=ymin,ymax=ymax),fill="#0F8E87", alpha=0.15)+
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
  filter(., alpha_0==.2, S %in% c(.73),Scena=="Global facilitation",
         delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
  group_by(., cintra,alpha_0,S,delta,Scena)%>%
  summarise(.,.groups ="keep",c11=mean(c11) )


p2=ggplot(d_clust_mean)+
  geom_point(aes(x=delta,y=c11),size=3,shape=1)+
  geom_line(aes(x=delta,y=c11,group=interaction(S)),lwd=.5)+
  the_theme+labs(x=TeX('$\\delta$'),color="",y="Stress-tolerant clustering")+
  geom_hline(yintercept = 1)+
  ggtitle(TeX(r'($\S = 0.73, \alpha_e = 0.2)'))


# density of pairs
d_pairs=d_clustering%>%
  filter(., alpha_0==.2, S %in% c(0),Scena=="Global facilitation",
         delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
  melt(., measure.vars=c("Rho_20","Rho_10"))%>%
  mutate(.,variable=recode_factor(variable,"Rho_10"="Pair Stress-tol. sp./Fertile",
                                  "Rho_20"="Pair Competitive sp./Fertile"))
p3=ggplot(d_pairs)+
  geom_point(aes(x=delta,y=value,color=variable),size=3,shape=1)+
  geom_line(aes(x=delta,y=value,color=variable,group=interaction(variable)),lwd=.5)+
  ggtitle(TeX(r'($\S = 0, \alpha_e = 0.2)'))+
  the_theme+
  labs(x=TeX('$\\delta$'),y=expression(paste("Pair cover")),color="")+
  scale_color_manual(values=c("#858BE0",as.character(color_Nsp(5)[c(2)])))+
  theme(legend.box = "vertical")+
  guides(color=guide_legend(ncol=1))

p_bottom=ggarrange(p3,
                   ggarrange(p2,ggplot()+theme_void(),nrow=2,heights = c(1,.25)),
                   p1,ncol=3,labels = letters[1:3],align = "h")

ggsave(filename = "../Figures/SI/Mechanisms_global_facilitation.pdf",width = 7,height = 4)


## >> Multistability, varying the trade-off shape ----

"
Same as figure 4 but with different shapes of the trade-off
Please run Step 2.5 of Dryland_shift_main.R prior to generate the data needed
"

#first: explaining trade-off


d=expand.grid(Psi=seq(0,1,length.out=20),Gamma=c(.5,1,1.5))
d$Psi2=sapply(1:nrow(d),function(x){
  return(d$Psi[x]^d$Gamma[x])
})

p=ggplot(d)+
  geom_line(aes(x=Psi,y=Psi2,color=as.factor(Gamma)))+
  geom_point(aes(x=Psi,y=Psi2,color=as.factor(Gamma)),shape=1,size=3,fill="white")+
  geom_segment(aes(x = .5, y = .15,
                   xend = 1, yend = .15),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(aes(x = .75, y = .22,label="Stress-tolerance"),stat = "unique")+
  geom_text(aes(x = .75, y = .08,label="Facilitation"),stat = "unique")+
  geom_segment(aes(x = 0.5, y = .8,
                   xend = 0, yend = .8),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(aes(x = .25, y = .9,label="Tolerance to competition"),stat = "unique")+
  the_theme+
  labs(x=TeX("$\\psi$"),y=TeX("$\\psi^{\\gamma}$"),color=TeX("$\\gamma$"))+
  scale_color_viridis_d()
ggsave( "../Figures/SI/Trade_off_shape.pdf",p,width = 6,height = 4)



#second: multistability
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()

for (tradeoff in c(.5,1.5)){
  for (stress in stress_seq){
    
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Degradation_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Restoration_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"))
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
d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1,d2$Tradeoff),]

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
  scale_fill_gradientn(colors=colorRampPalette(c("#EEF0F5","#AEAFB1","#77787D","black"))(100))+
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
           guides(color = guide_legend(override.aes = list(size = 1.5)))+
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
           guides(color = guide_legend(override.aes = list(size = 1.5)))+
           new_scale_color() +
           
           geom_line(data=d_state%>%
                       filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==1.5),
                     aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
           scale_color_manual(values=c("white","black"))+
           guides(color=F)
  )  
  
  
}





pbifu=ggarrange(ggarrange(ggarrange(p2_1,p2_2+ylab(""), labels = c("(i)","(ii)"),
                                    ncol=2,widths = c(1.2,1),hjust = c(-4.5,-3.5),vjust=4,
                                    legend="none"),
                          ggarrange(p3_1+ylab(""),p3_2+ylab(""), labels = c("(i)","(ii)"),
                                    ncol=2,widths = c(1,1),hjust = c(-4.5,-6),vjust=4,
                                    legend="none"),
                          ncol=2,
                          legend = "none",
                          widths = c(1.1,1),labels = c("b","c")),
                ggarrange(ggplot()+theme_void(),
                          get_legend(p2_1),
                          ggplot()+theme_void(),ncol=3,widths = c(.5,2,.5)),nrow=2,heights = c(4,1))


p_tot=ggarrange(
  ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.15,1,.15)),
  pbifu,nrow=2,heights =  c(2,1),labels = c(letters[1],""),hjust = c(-9,0))

ggsave(filename = "../Figures/SI/Trade_off_multistability_PA.pdf",plot = p_tot,width = 9,height = 7)





## >> Comparing degradation and restoration points trade-off ----


#convex and concave trade-off
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()

for (tradeoff in c(.5,1.5)){
  for (stress in stress_seq){
    
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Degradation_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Restoration_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"))
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

d_final1=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      for (tradeoff in unique(d_state$Tradeoff)){
        d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0,Tradeoff==tradeoff)
        if (i>j){
          d_final1=rbind(d_final1,tibble(Psi1=i,Psi2=j,alpha_0=a0,Tradeoff=tradeoff,
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



#linear trade-off
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",round(stress,4),"_delta_",.1,"_facilitation_0.9_scalefacilitation_local.csv"),sep=";"),
           read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",round(stress,4),"_delta_",.1,"_facilitation_0.9_scalefacilitation_local.csv"),sep=";"))
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

d_final2=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      if (i>j){
        d_final2=rbind(d_final2,tibble(Psi1=i,Psi2=j,alpha_0=a0,
                                       Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                       Frac_multi_raw=sum(d_fil$multistab),
                                       Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
                                       Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))]
        ))
      }
    }
  }
}

#now comparing degradation and restoration points
#making sure that they are all sorted in the same way

d_final2=d_final2[order(d_final2$Psi1,d_final2$Psi2,d_final2$alpha_0),]
d_convex=filter(d_final1,Tradeoff==.5)
d_concave=filter(d_final1,Tradeoff==1.5)

d_final2=d_final2[order(d_final2$Psi1,d_final2$Psi2,d_final2$alpha_0),]


#now making the difference between non-linear trade-off and linear trade-off

d_convex$Deg=d_convex$Deg-d_final2$Deg
d_concave$Deg=d_concave$Deg-d_final2$Deg
d_convex$Rest=d_convex$Rest-d_final2$Rest
d_concave$Rest=d_concave$Rest-d_final2$Rest
d_convex$Frac_multi_raw=d_convex$Frac_multi_raw-d_final2$Frac_multi_raw
d_concave$Frac_multi_raw=d_concave$Frac_multi_raw-d_final2$Frac_multi_raw

d_tot=rbind(d_convex,d_concave)

#now plotting

for (i in 1:2){
  assign(paste0("p1_",i),
         
         ggplot(d_tot%>%filter(., Tradeoff==c(.5,1.5)[i]))+
           geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
           the_theme+
           facet_grid(.~Tradeoff,labeller = label_bquote(cols = gamma == .(Tradeoff)))+
           labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
                y=TeX(r'(Species 2 trait, $\psi_2$)'),
                fill="")+
           scale_fill_gradient2(low = "red",mid="white",high="blue",midpoint = 0)+
           guides(shape="none")+
           theme(strip.text.x = element_text(size=13),strip.background.x = element_blank())
  )
}





param_space=expand.grid(point=c("Degradation point","Restoration point"),
                        shape=c(.5,1.5))

for (i in 1:4){
  
  assign(paste0("p2_",i),
         ggplot(d_tot%>%
                  melt(., measure.vars=c("Rest","Deg"))%>%
                  mutate(., variable=recode_factor(variable,
                                                   "Rest" = "Restoration point",
                                                   "Deg" = "Degradation point"))%>%
                  filter(., variable==param_space$point[i],Tradeoff==param_space$shape[i]))+
           geom_tile(aes(x=Psi1,y=Psi2,fill=value))+
           the_theme+
           facet_grid(.~Tradeoff,labeller = label_bquote(cols = gamma == .(Tradeoff)))+
           labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
                y=TeX(r'(Species 2 trait, $\psi_2$)'),
                fill="")+
           scale_fill_gradient2(low = "red",mid=alpha("white",.5),high="blue")+
           guides(shape="none")+
           theme(strip.text.x = element_text(size=13),strip.background.x = element_blank())
  )
  
}


p_tot=ggarrange(
  ggarrange(p1_1+ggtitle("Community bistability"),p1_2,labels = c(letters[1],""),align = "hv"),
  ggarrange(p2_1+ggtitle("Extinction point"),p2_3,labels = c(letters[2],""),align = "hv"),
  ggarrange(p2_2+ggtitle("Restoration point"),p2_4,labels = c(letters[3],""),align = "hv"),
  nrow=3
)


ggsave("../Figures/SI/Comparizon_linear_non_linear_tradeoff.pdf",width = 7,height = 10)





# 
## >> Restoration & extinction points trait space----
# 
# 
# 
# stress_seq=seq(0,.82, length.out = 100)
# 
# d=tibble()  
# d_bistab=tibble()
# 
# 
# for (stress in stress_seq){
#   
#   
#   
#   d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
#                              .3,"_branch_Degradation_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"),
#            read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
#                              .3,"_branch_Restoration_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"))
#   d=rbind(d,d2)
#   
#   
#   
# }
# 
# d2=d
# d2[,1:2][d2[,1:2] < 10^-3] = 0
# d2$state = sapply(1:nrow(d2), function(x) {
#   if (d2[x, 1] > 0 & d2[x, 2] > 0) {
#     return("Coexistence")
#   }
#   if (d2[x, 1] > 0 & d2[x, 2] == 0) {
#     return("Stress_tolerant")
#   }
#   if (d2[x, 1] == 0 & d2[x, 2] > 0) {
#     return("Species 2")
#   }
#   if (d2[x, 1] == 0 & d2[x, 2] == 0) {
#     return("Desert")
#   }
# })
# d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1),]
# 
# all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
#   if (d2$state[x] != d2$state[x+1]){
#     return(paste0(d2$state[x],"/", d2$state[x+1]))
#   }
#   else {return(d2$state[x])}
# })
# 
# d_state=d2%>%
#   filter(., Branches=="Degradation")%>%
#   select(.,-Branches)
# d_state$all_state=all_state
# 
# d_state$multistab=sapply(1:nrow(d_state),function(x){
#   if (d_state$all_state[x] %in% c( "Species 2/Stress_tolerant","Coexistence/Stress_tolerant","Species 2/Coexistence",
#                                    "Stress_tolerant/Species 2")){
#     return(1)
#   } else {return(0)}
#   
# })
# 
# d_state$Stress_tolerant_reg=d2$Stress_tolerant[which(d2$Branches=="Restoration")]
# d_state$Competitive_reg=d2$Competitive[which(d2$Branches=="Restoration")]
# 
# d_final=tibble()
# for (i in unique(d_state$Psi1)){
#   for (j in unique(d_state$Psi2)){
#     for (a0 in unique(d_state$alpha_0)){
#       d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
#       
#       if (i>j){
#         d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,
#                                      Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
#                                      Frac_multi_raw=sum(d_fil$multistab),
#                                      Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
#                                      Deg2=d_fil$Stress[max(which(d_fil$Stress_tolerant>0))+1],
#                                      Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))],
#                                      Niche_1 = length(which(d_fil$Stress_tolerant>0))/length(which(d_fil$state!="Desert"))  ,
#                                      Niche_2 = length(which(d_fil$Competitive>0))/length(which(d_fil$state!="Desert"))
#         ))
#       }
#     }
#   }
# }
# 
# 
# 
# # Niche least competitive species
# p1=ggplot(d_final)+
#   geom_tile(aes(x=Psi1,y=Psi2,fill=Niche_1))+
#   the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_2$)'),
#                  fill="Fraction of vegetation states \n where species 1 has a positive cover")+
#   scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100))+
#   guides(shape=F)
# 
# 
# ggsave("../Figures/SI/Niche_least_competitive_species.pdf",p1,width = 6,height = 5)
# 
# 
# #Degradation and restoration points
# 
# p1=ggplot(d_final)+
#   geom_tile(aes(x=Psi1,y=Psi2,fill=Deg))+
#   the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_1$)'),fill="Species 2 extinction point")+
#   scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100))+
#   guides(shape=F)
# 
# 
# 
# #Restoration point traits
# p2=ggplot(d_final)+
#   geom_tile(aes(x=Psi1,y=Psi2,fill=Rest))+
#   the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_2$)'),fill="Species 2 restoration point")+
#   scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100),
#                        breaks=c(0,.1,.2))+
#   geom_line(data=tibble(x=c(0,1,1,1),y=c(0,0,0,1),group=c(1,1,2,2)),aes(x=x,y=y,group=group),lwd=1)+
#   geom_text(data=tibble(x=c(.5,1.05),
#                         y=c(0.05,.55),
#                         txt=c("(c)","(d)")),aes(x=x,y=y,label=txt),size=4.5)+
#   guides(shape=F)
# 
# 
# 
# p3=ggplot(d_final%>%filter(Psi2==0))+
#   geom_smooth(aes(x=Psi1,y=Rest),se = F,color="black")+ggtitle(TeX("$\\psi_2 = 0$"))+
#   the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y="Species 2 restoration point")
# 
# p4=ggplot(d_final%>%filter(Psi1==unique(d_final$Psi1)[49]))+
#   geom_smooth(aes(x=Psi2,y=Rest),se = F,color="black")+ggtitle(TeX("$\\psi_1 = 1$"))+
#   the_theme+labs(x=TeX(r'(Species 2 trait, $\psi_2$)'),y="Species 2 restoration point")
# 
# p_bottom=ggarrange(p3,p4,labels = letters[3:4])
# p_tot=ggarrange(ggarrange(p1,p2,ncol=2,labels = letters[1:2]),
#                 p_bottom,nrow=2,heights = c(1.5,1))
# 
# ggsave("../Figures/SI/Restoration_degration.pdf",p_tot,width = 9,height = 7)
# 
# 
# 







## >> Landscapes along dispersal gradient ----

"
Illustration of how dispersal modulates competitive hierarchy between competitive and stress-tolerant 
species.
Please run chunk 5 of the julia file Dryland_shift_Nspecies_function.jl prior to generate the data needed
"

list_landscape=list.files("../Table/2_species/CA/",pattern = "Dispersal")
for (i in 1:length(list_landscape)){
  landscape=read.table(paste0("../Table/2_species/CA/",list_landscape[i]),sep=",")
  assign(paste0("p_land_",i),Plot_landscape(landscape,Nsp = 2)+theme(legend.position = "none",
                                                                     plot.title = element_text(size=12)))
}
p=ggarrange(p_land_1+ggtitle(TeX("$\\delta = 0$")),p_land_2+ggtitle(TeX("$\\delta = 0.1$")),
            p_land_3+ggtitle(TeX("$\\delta = 0.2$")),p_land_4+ggtitle(TeX("$\\delta = 0.3$")),
            p_land_5+ggtitle(TeX("$\\delta = 0.7$")),p_land_6+ggtitle(TeX("$\\delta = 1$")),nrow=2,ncol=3)


ggsave("../Figures/SI/Dynamics_landscape_dispersal.pdf",p,width = 7,height = 4)

## >> Bifurcation diagrams effect of competition ----

"
Code to illustrate the change in the bifurcation diagram with the interspecific competition
"


julia_setup()
de = diffeq_setup()

tspan = c(0, 5000) #to avoid long transient
t = seq(0, 5000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
# c_seq=c(0,.4)
c_seq=c(.25,.3,.35)
name_scena=c("global_C_local_F")
delta_seq=c(.1,.9)
branches=c("Degradation","Restoration")
type_ini_cond=c("equal","low_ST","low_comp")[1]

for (branch in branches){
  for (type_ini in type_ini_cond){
    
    if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
      if (type_ini=="equal") state =Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
      if (type_ini=="low_ST") state =Get_PA_initial_state(Get_MF_initial_state(c(.01,.79,.1)))
      if (type_ini=="low_comp") state =Get_PA_initial_state(Get_MF_initial_state(c(.79,.01,.1)))
      S_seq=seq(0,1, length.out = N_sim)
    }else {
      state =Get_PA_initial_state(Get_MF_initial_state(c(.0001,.0099,.49)))
      
      
      if (type_ini=="equal") state =Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
      if (type_ini=="low_ST") state =Get_PA_initial_state(Get_MF_initial_state(c(.0001,.0099,.49)))
      if (type_ini=="low_comp") state =Get_PA_initial_state(Get_MF_initial_state(c(.0099,.0001,.49)))
      
      S_seq=rev(seq(0,1, length.out = N_sim))
    }
    
    for (disp in delta_seq) {
      
      
      for (ccomp in c_seq){
        
        for (S in S_seq) { #varying dispersal scale
          
          julia_assign("state", state)
          param=Get_PA_parameters()
          param["cintra"]=.3
          param["S"] = S
          param["alpha_0"] = ccomp
          param["delta"]=disp
          julia_assign("p", param)
          
          
          julia_assign("p", param)
          prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = ccomp,type_ini,
                                                     Delta=disp,Branch=branch))
        }
      }
    }
  }
}      
d2[d2 < 10^-4] = 0


p1=ggplot(d2%>%
            filter(., type_ini=="equal",Branch=="Degradation",Delta==.1)%>%
            melt(., measure.vars=c("rho_1","rho_2"))%>%
            mutate(., variable=recode_factor(variable,"rho_1"="Stress-tolerant","rho_2"="Competitive")))+
  labs(x="Stress, S",y="Species cover",color="",linetype="")+
  ggtitle("Degradation: high initial cover")+
  geom_line(aes(x=S,y=value,color=variable))+
  scale_color_manual(values=as.character(color_rho[c(4,2)]))+
  the_theme+theme(legend.position = "none")+
  facet_grid(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)))

p2=ggplot(d2%>%
            filter(., type_ini=="equal",Branch=="Degradation",Delta==.9)%>%
            melt(., measure.vars=c("rho_1","rho_2"))%>%
            mutate(., variable=recode_factor(variable,"rho_1"="Stress-tolerant","rho_2"="Competitive")))+
  labs(x="Stress, S",y="Species cover",color="",linetype="")+
  ggtitle("Degradation: high initial cover")+
  geom_line(aes(x=S,y=value,color=variable))+
  scale_color_manual(values=as.character(color_rho[c(4,2)]))+
  the_theme+theme(legend.position = "none")+
  facet_grid(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)))



p3=ggplot(d2%>%
            filter(., type_ini=="equal",Branch=="Restoration",Delta==.1)%>%
            melt(., measure.vars=c("rho_1","rho_2"))%>%
            mutate(., variable=recode_factor(variable,"rho_1"="Stress-tolerant","rho_2"="Competitive")))+
  ggtitle("Restoration: low initial cover")+
  labs(x="Stress, S",y="Species cover",color="",linetype="")+
  geom_line(aes(x=S,y=value,color=variable))+
  scale_color_manual(values=as.character(color_rho[c(4,2)]))+
  the_theme+theme(legend.position = "none")+
  facet_grid(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)))

p4=ggplot(d2%>%
            filter(., type_ini=="equal",Branch=="Restoration",Delta==.9)%>%
            melt(., measure.vars=c("rho_1","rho_2"))%>%
            mutate(., variable=recode_factor(variable,"rho_1"="Stress-tolerant","rho_2"="Competitive")))+
  ggtitle("Restoration: low initial cover")+
  labs(x="Stress, S",y="Species cover",color="",linetype="")+
  geom_line(aes(x=S,y=value,color=variable))+
  scale_color_manual(values=as.character(color_rho[c(4,2)]))+
  the_theme+theme(legend.position = "none")+
  facet_grid(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)))

p_tot = ggarrange(p1+theme(strip.text.x = element_text(size = 10)), 
                  p3+theme(strip.text.x = element_text(size = 10)), 
                  p2+theme(strip.text.x = element_text(size = 10)),
                  p4+theme(strip.text.x = element_text(size = 10)),
                  labels = letters[1:4],heights = c(1.2,1.2,1,1),
                  ncol=2,nrow=2,common.legend = T,legend="none")

p_tot_w_legend=ggarrange(p_tot,get_legend(p4),nrow=2,heights = c(10,1))
ggsave("../Figures/SI/Mechanism_priority.pdf",p_tot_w_legend,width=10,height = 5)





## >> Net effect between and within species before the shift ----

"
Compute the net effect and within species effect before the community turnover happens
"

tspan = c(0, 30000)
t = seq(0, 30000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
julia_assign("state", state)

N_sim = 60;epsilon=10^(-8)
C_for_analyse = c(.25,.3,.35)
S_seq = seq(0, 1, length.out = N_sim)
type_seq=c("self","mutual")

d2 = tibble()
param=Get_PA_parameters()  
param=c(param[1:3],"beta1"=1,"beta2"=1,param[5:12])

for (comp in C_for_analyse) {
  
  for (type in type_seq){
    
    for (S in S_seq) {
      
      for (sp in 1:2) { #for both species
        
        param["S"] = S
        param["alpha_0"] = comp # update parameters
        param["delta"] = .1
        julia_assign("p", param)
        
        
        # first without press perturbation
        
        prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
        
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        
        if (type=="mutual"){
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-3000):nrow(d),c(1,2)[-sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],
                                  Disp=.1,Direction=type)) 
          
        } else {
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-3000):nrow(d),c(1,2)[sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "control",Species=c(1,2)[sp],
                                  Disp=.1,Direction=type)) 
          
        }
        
        # with press perturbation on growth rate proxy
        
        param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
        julia_assign("p", param)
        
        
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F_press, state, tspan, p)")
        
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        if (type=="mutual"){
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-3000):nrow(d),c(1,2)[-sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],
                                  Disp=.1,Direction=type))
          
        } else {
          d2 = rbind(d2, as_tibble(mean(d[(nrow(d)-3000):nrow(d),c(1,2)[sp]])) %>%
                       add_column(S = S, alpha_0 = comp, Type = "press",Species=c(1,2)[sp],
                                  Disp=.1,Direction=type))
          
        }
        param[paste0("beta",sp)] = 1
      }
    }
  }
}


d2[d2 < epsilon] = 0
colnames(d2) = c("Eq", "S", "alpha_0", "Type","Species","Disp","Direction")


net_effect =sapply(seq(1, nrow(d2) , by = 2),function(x){
  return((d2$Eq[x+1] - d2$Eq[x]) / epsilon)
})

d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)
d_net$value=net_effect

index=1
for (i in c(.25,.3,.35)){
  for (direction in c("mutual","self")){
    assign(paste0("p_",index),
           ggplot(NULL)+
             geom_line(data=d_net%>%filter(., Direction==direction,alpha_0==i)%>%
                         mutate(., Species=recode_factor(Species,'1'=ifelse(direction=="mutual",
                                                                            "Effect of ST sp. on Comp. sp.",
                                                                            "Effect of ST sp. on ST. sp."),
                                                         '2'=ifelse(direction=="mutual",
                                                                    "Effect of Comp. sp. on ST sp.",
                                                                    "Effect of Comp. sp. on Comp. sp."))),
                       aes(x=S,y=value,color=as.factor(Species)))+
             geom_point(data=d_net%>%filter(., Direction==direction,alpha_0==i)%>%
                          mutate(., Species=recode_factor(Species,'1'=ifelse(direction=="mutual",
                                                                             "Effect of ST sp. on Comp. sp.",
                                                                             "Effect of ST sp. on ST. sp."),
                                                          '2'=ifelse(direction=="mutual",
                                                                     "Effect of Comp. sp. on ST sp.",
                                                                     "Effect of Comp. sp. on Comp. sp.")))%>%
                          filter(., S %in% seq(0,1,length.out=N_sim)[seq(1,N_sim,length.out=12)]),
                        aes(x=S,y=value,color=as.factor(Species)),size=2,shape=1)+
             scale_color_manual(values=rev(c("#BAE886", "#84BAF3")))+
             labs(x="Stress (S)",y="Net effect",color="")+
             the_theme+
             geom_hline(yintercept = 0,linetype=9)
           
           
           
    )
    index=index+1
  }
}



#balance inter/intra

d_net=d2%>%
  filter(., Type=="control")%>%
  select(.,-Type)

d_net$value=net_effect
d_net$Species=sapply(1:nrow(d_net),function(x){
  
  if (d_net$Species[x]==1 & d_net$Direction[x]=="mutual"){
    return(2)
  }else if( d_net$Species[x]==2 & d_net$Direction[x]=="mutual"){
    return(1)
  }else {
    return(d_net$Species[x])
  }
  
})
d_net=d_net%>%group_by(., alpha_0,Species,Disp,S)%>%
  summarise(., .groups = "keep",value=sum(value))

index=7
for (i in c(.25,.3,.35)){
  for (direction in c("self")){
    assign(paste0("p_",index),
           ggplot(NULL)+
             geom_line(data=d_net%>%filter(., alpha_0==i)%>%
                         mutate(., Species=recode_factor(Species,'1'=ifelse(direction=="mutual",
                                                                            "Effect of ST sp. on Comp. sp.",
                                                                            "Effect of ST sp. on ST. sp."),
                                                         '2'=ifelse(direction=="mutual",
                                                                    "Effect of Comp. sp. on ST sp.",
                                                                    "Effect of Comp. sp. on Comp. sp."))),
                       aes(x=S,y=value,color=as.factor(Species)))+
             geom_point(data=d_net%>%filter(.,alpha_0==i)%>%
                          mutate(., Species=recode_factor(Species,'1'=ifelse(direction=="mutual",
                                                                             "Effect of ST sp. on Comp. sp.",
                                                                             "Effect of ST sp. on ST. sp."),
                                                          '2'=ifelse(direction=="mutual",
                                                                     "Effect of Comp. sp. on ST sp.",
                                                                     "Effect of Comp. sp. on Comp. sp.")))%>%
                          filter(., S %in% seq(0,1,length.out=N_sim)[seq(1,N_sim,length.out=12)]),
                        aes(x=S,y=value,color=as.factor(Species)),size=2,shape=1)+
             scale_color_manual(values=rev(c("#BAE886", "#84BAF3")))+
             labs(x="Stress (S)",y="Intra-inter",color="")+
             the_theme+
             geom_hline(yintercept = 0,linetype=9)
           
           
           
    )
    index=index+1
  }
}

p_tot=ggarrange(ggarrange(p_1+ggtitle(TeX("$\\alpha_e = 0.25$")),p_3+ggtitle(TeX("$\\alpha_e = 0.3$")),
                p_5+ggtitle(TeX("$\\alpha_e = 0.35$")),common.legend = T,legend = "bottom",ncol=3,labels = c(letters[1],"","")),
                ggarrange(p_2+ggtitle(TeX("$\\alpha_e = 0.25$")),
                p_4+ggtitle(TeX("$\\alpha_e = 0.3$")),p_6+ggtitle(TeX("$\\alpha_e = 0.35$")),
                ncol=3,common.legend = T,legend="bottom",labels = c(letters[2],"","")),
                ggarrange(p_7+ggtitle(TeX("$\\alpha_e = 0.25$")),
                          p_8+ggtitle(TeX("$\\alpha_e = 0.3$")),
                          p_9+ggtitle(TeX("$\\alpha_e = 0.35$")),
                          ncol=3,common.legend = T,legend="bottom",labels = c(letters[3],"","")),nrow=3)



ggsave("../Figures/SI/Net_effet_prior_shift.pdf",p_tot,width = 8,height = 7.5)




## >> State diagram lower degradation, higher restoration ----

"
Same as figure 3 but with a higher degradation rate compare to te restoration rate
Please run Step 2.6 of Dryland_shift_main.R prior to generate the data needed
"


d2=read.table(paste0("../Table/2_species/PA/Higher_restor_lower_deg.csv"),sep=";")
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

p=ggplot(d2t%>%
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


ggsave("../Figures/SI/Lower_degradation_higher_restoration.pdf",p,width = 9,height = 5)


#*****************************************************************


#---------------------   SI figures N-species  --------------------------
## >> Scaling-up 5sp, equal initial conditions ----


d_mf_equal=read.table("../Table/N_species/MF/MF_equal_ini.csv",sep=",")
d_mf_equal=d_mf_equal[,c(1:5,(ncol(d_mf_equal)-2):ncol(d_mf_equal))]
colnames(d_mf_equal)=c(paste0("Sp_",1:5),"alpha_0", "Branches", "Stress")
d_melt=d_mf_equal%>%
  melt(., measure.vars=paste0("Sp_",1:5))

d_melt$variable=as.character(d_melt$variable)

d_melt$Trait=sapply(1:nrow(d_melt),function(x){
  seq(1,0,length.out=5)[as.numeric(strsplit(d_melt$variable[x],split = "_")[[1]][2])]
})

p1=ggplot(d_melt%>%
            mutate(., Branches=recode_factor(Branches,"1"="Degradation","2"="Restoration"))) +
  geom_line(aes(x=Stress,y=value,color=Trait,group=interaction(Trait,alpha_0,Branches),linetype=as.factor(Branches)))+
  facet_wrap(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)),scales = "free")+
  the_theme+
  labs(x="Stress, (S)", y="Species cover",linetype="",color="")+
  scale_color_gradientn(colours = (color_Nsp(5)))+
  theme(strip.text.x = element_text(size=12))

ggsave("../Figures/SI/Equal_ini_MF.pdf",p1,width=8,height=3)






## >> Colored by number of species ----

"
Same as figure 6a-b but colored with the number of species coexisting.
As for Figure 6, you need to run the Step 3.1 of the Dryland_shift_main.R R file
to generate the data needed here
"


d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")

for (i in 1:4){
  for (j in 1:2){
    assign(paste0("p_",i,"_",j),
           ggplot(d_tot%>%filter(., Nsp %in% c(5,25)[i],Competition %in% c(.225,.30)[j],Branch==1))+
             geom_point(aes(x=Stress,y=Rho_plus,color=Nb_sp,fill=Nb_sp),size=.3,shape=21)+
             the_theme+labs(y="Vegetation cover",color="")+
             scale_color_viridis_c(option = "E")+
             scale_fill_viridis_c(option = "E")+
             guides(fill="none")+
             labs(x="Stress, S",color="Number of species \n in the community  ")+
             facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="# species"==.(Nsp)))+
             theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
                   panel.background = element_blank(),strip.background.y = element_blank())
           )
  }
}
p=ggarrange(p_1_1+labs(x=""),p_1_2,p_2_1+labs(x="",y=""),p_2_2+labs(y=""),nrow=2,ncol=2)

ggsave("../Figures/SI/Nsp_CSI_colored_nbspecies.pdf",p,width = 7,height =6 )


## >> Main fig N-species with mean-field model ----

"
Same as figure 6 but with the mean-field model
As for Figure 6, you need to run the Step 3.2 of the Dryland_shift_main.R R file
to generate the data needed here
"


d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")
d_tot=d_tot[-which(d_tot$Rho_plus==0 & d_tot$Stress<.7),]

p1=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.30),Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.1,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())

p2=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.30),Branch==1))+
  geom_point(aes(x=Stress,y=Rho_plus,color=Psi_normalized,fill=Psi_normalized),size=.3,shape=21)+
  the_theme+labs(y="Vegetation cover",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())




d=read.table("../Table/N_species/MF/post_proc_sim1.csv",sep=";")



p3=ggplot(d)+
  geom_smooth(aes(x=Stress,y=N_ASS,color=as.factor(Competition),group=Competition),
              se = F)+
  the_theme+
  facet_grid(Nsp~.,labeller = label_bquote(rows= "# species" ==.(Nsp)))+
  labs(x="Stress (S)",y="# of alternative states",color=TeX("$\\alpha_e \ \ $"))+
  scale_color_viridis_d()+
  theme(strip.text.y = element_text(size=13),panel.background = element_blank(),
        strip.background.y = element_blank(),legend.title = element_text(size=14))




Fig_6_SI=ggarrange(ggarrange(p2,p1,common.legend = T,legend = "bottom",nrow=2,labels = letters[1:2]),
                ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),nrow=3,heights = c(.05,1.5,.05),labels = c("",letters[3],"")),ncol=2,widths =  c(1.2,1))

ggsave("../Figures/SI/N_species_CSI_cover_MF.pdf",Fig_6_SI,width = 9,height = 8)


## >> Types of bistability MF ----

type_bistab=read.table("../Table/N_species/MF/post_proc_sim2.csv",sep=";")
  
for (i in 1:3){
  assign(paste0("p1_",i),
         ggplot(type_bistab%>%filter(., Nsp==unique(type_bistab$Nsp)[i],Competition>.2)%>%
                  mutate(., Competititon=as.factor(Competition),Stress=round(Stress,4)))+
           geom_tile(aes(x=Stress,y=Competition,fill=Type))+
           facet_wrap(.~Nsp,scales = "free",labeller = label_bquote(cols='# of species'==.(Nsp)))+
           the_theme+
           labs(fill="",x="Stress (S)")+
           theme(strip.background.x = element_blank())+
           scale_fill_manual(values=c("Degraded"="#000000","Env"="#0F8E87","Mutual exclusion"="#C8A4D0",
                                      "No bistab"="gray50","Cliques"="#EFE8BC"),
                             labels=c("Degraded","Environmental","Mutual exclusion",
                                      "No bistability","Cliques"))+
           theme(strip.text.x = element_text(size=12)))
}

pA=ggarrange(p1_3+theme(legend.spacing.x = unit(.5, 'cm')),
             p1_1+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()),
             p1_2+theme(axis.text.y = element_blank(),axis.ticks.y = element_blank(),axis.title.y = element_blank()),
             ncol=3,common.legend = T,legend="bottom",widths = c(1.3,1,1))

d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")%>%
  dplyr::group_by(.,Competition,Stress,Nsp)%>%
  dplyr::summarise(., .groups = "keep",mean_sp_div=mean(Nb_sp))%>%
  arrange(., Stress,Nsp,Competition)

d_tot$mean_sp_div[d_tot$mean_sp_div==0]=NA #for specific black color

#2D param space with types of multistability

pB=ggplot(d_tot%>%filter(., Competition>.15)%>%
            mutate(., Competition=as.factor(Competition),Stress=round(Stress,4)))+
  geom_tile(aes(x=Stress,y=Competition,fill=mean_sp_div),alpha=.8)+
  the_theme+
  facet_wrap(.~Nsp,labeller=label_bquote(cols="# of species "==.(Nsp)))+
  labs(fill="# of coexisting species  ",x="Stress (S)",y=TeX(r'(Interspecific competition, \ $\alpha_e)'))+
  theme(strip.background.x = element_blank())+
  scale_fill_viridis_c(na.value = "black")+
  theme(strip.text.x = element_text(size=12))



#Cover cliques vs 1 species
d_cover_cliques=read.table("../Table/N_species/MF/post_proc_sim3.csv",sep=";")

pC=ggplot(d_cover_cliques%>%
            filter(., Competition>.15)%>%
            mutate(., Competition=as.factor(Competition))%>%
            filter(., Nb_sp>1,!is.na(Psi_normalized),Type=="Cliques")%>%
            group_by(., Competition,Nsp)%>%
            summarise(., .groups = "keep",
                      q1=mean(Nb_sp)-sd(Nb_sp),q3=mean(Nb_sp)+sd(Nb_sp),
                      q2=mean(Nb_sp)))+
  geom_pointrange(aes(x=(Competition),y=q2,ymax=q3,ymin=q1,fill=Competition),size=.75,shape=24,color="black")+
  the_theme+
  facet_wrap(.~Nsp,labeller=label_bquote(cols="# of species "==.(Nsp)))+
  scale_y_continuous(breaks = c(1,3,5,7),limits = c(1,7))+
  scale_fill_manual(values=colorRampPalette(c("#1A41AF","#72B2D6","#B9D7E8","#FBEFCB","#F5CB61","#D26F3C"))(7))+
  labs(x=TeX(r'(Interspecific competition, \ $\alpha_e)'),y="mean # of species \n within the clique")+
  theme(legend.position = "none",strip.background.x = element_blank(),strip.text.x = element_text(size=12))



p_tot=ggarrange(pA,pB,pC,nrow=3,labels=letters[1:3],hjust=-5)

ggsave("../Figures/SI/Type_bistability_MF.pdf",p_tot,width = 8,height = 10)

## >> Example of species succession gradient stress N-species MF ----


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

p=ggarrange(p_1+ggtitle(TeX(r'(N species = 5, \ $\alpha_e = 0.3)'))+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()),
            p_2+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.25)'))+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank(),
                    axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            p_3+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.325)')),
            p_5+ggtitle(TeX(r'(N species = 25, \ $\alpha_e = 0.25)'))+
              theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            nrow=2,ncol=2,heights = c(1,1.1),common.legend = T,legend="bottom")

ggsave("../Figures/SI/Example_bifu_Nspecies.pdf",width = 7,height = 5)


## >> Example of species succession gradient stress N-species PA ----


Fig_to_plot=tibble(Nsp=c(5,5,15,25),Number_plot=c(492,498,36,5))


for (i in 1:nrow(Fig_to_plot)){
  
  Nsp=Fig_to_plot$Nsp[i]
  list_csv=list.files(paste0('../Table/N_species/PA/',Nsp,'_sp/'))
  d2=read.table(paste0("../Table/N_species/PA/",Nsp,"_sp/",list_csv[Fig_to_plot$Number_plot[i]]),sep=",")
  colnames(d2)=c(paste0(1:Nsp),"Random_ini","Competition","Facilitation","Branch","Stress")
  
  d2[d2<10^(-4)]=0
  
  
  
  
  assign(paste0("p_",i),
         ggplot(d2%>%
                  melt(., measure.vars=paste0(1:Nsp))%>%
                  mutate(., Trait=sapply(1:nrow(.),function(x){
                    return(seq(1,0,length.out=Nsp)[as.numeric(.$variable[x])])
                  }))%>%
                  mutate(., Branch=recode_factor(Branch, "1"="High initial cover","2"= "Low initial cover")))+
           geom_line(aes(x=Stress,y=value,color=Trait,group=interaction(variable,Branch),linetype=as.factor(Branch)),size=.8)+
           the_theme+labs(x="Stress, S",y="Species cover",color=TeX("$\\psi$   "),linetype="")+
           scale_color_gradientn(colors=color_Nsp(Nsp))+
           theme(legend.box = "vertical")
  )
}

p=ggarrange(p_1+ggtitle(TeX(r'(N species = 5, \ $\alpha_e = 0.325)'))+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()),
            p_2+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.35)'))+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank(),
                    axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            p_3+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.35)')),
            p_4+ggtitle(TeX(r'(N species = 25, \ $\alpha_e = 0.325)'))+
              theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            nrow=2,ncol=2,heights = c(1,1.1),common.legend = T,legend="bottom")

ggsave("../Figures/SI/Example_bifu_Nspecies_PA.pdf",p,width = 7,height = 5)


## >> All CSI/Rho_+ bifurcation competition/Nsp ----

# Community index

d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")

p=ggplot(d_tot%>%filter(., Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=13),strip.text.y = element_text(size=12),
        panel.background = element_blank(),strip.background.x = element_blank())

ggsave("../Figures/SI/CSI_aij_Nsp.pdf",p,width = 12,height = 7)


#community scale vegetation cover

d_tot=read.table("../Table/N_species/PA/Multistability_CSI.csv",sep=";")

p=ggplot(d_tot%>%filter(., Branch==1))+
  geom_point(aes(x=Stress,y=Rho_plus,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=13),strip.text.y = element_text(size=12),
        panel.background = element_blank(),strip.background.x = element_blank())

ggsave("../Figures/SI/Vegetation_cover_aij_Nsp.pdf",p,width = 12,height = 7)


# Community index MF

d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

p=ggplot(d_tot%>%filter(., Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=13),strip.text.y = element_text(size=12),
        panel.background = element_blank(),strip.background.x = element_blank())

ggsave("../Figures/SI/CSI_aij_Nsp_MF.pdf",p,width = 12,height = 7)

