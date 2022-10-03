source('./N_species_Analysis_functions.R')

# Step 1 : N species MF -----
## 0) First test with 4 species ----

d_4sp=as_tibble(read.table("../Table/N_species/MF/4_species.csv",sep=","))%>%
  mutate(., V9=recode_factor(V9,"1"="Low_intra","2"="Intra=inter"))
colnames(d_4sp)=c(paste0("Sp_",1:4),"Fertile","Degraded","Stress","Rela_c","Scena")  
d_4sp$Rho_p=rowSums(d_4sp[,1:4])


d_4sp=transform(d_4sp,
                Rela_c = factor(Rela_c, levels=c(1,2.5,4), labels=c("frac(alpha[ji],alpha[ij]) : 1", "frac(alpha[ji],alpha[ij]) : 2.5",
                                                              "frac(alpha[ji],alpha[ij]) : 4")),
                Scena=factor(Scena,levels=c("Low_intra","Intra=inter"),labels=c("alpha[ij] < alpha[ii]","alpha[ij] > alpha[ii]")))


p=ggplot(melt(d_4sp,measure.vars = c(paste0("Sp_",1:4),"Rho_p")))+
  geom_line(aes(x=Stress,value,color=variable),lwd=1)+
  facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+
  the_theme+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=c(color_Nsp(4),"gray50"),labels=c("Sp_1"=TeX("$\\psi_1$"),"Sp_2"=TeX("$\\psi_2$"),
                                         "Sp_3"=TeX("$\\psi_3$"),"Sp_4"=TeX("$\\psi_4$"),
                                         "Rho_p"=TeX("$\\sum_j \\psi_j$")))

ggsave("../Figures/N_species/MF/4_species_scenario.pdf",p,width = 8,height = 5)





## 1) 15 species, varying community composition ----

### a) Bifurcation diagrams ----
  
Nsp=15

for (branch in c("Degradation","Restoration")){  
  for (com_compo in paste0(c(20,50,80),".0")){
    
    trait_sp = read.table(paste0("../Table/N_species/MF/TRAITS_",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=",")
    d_Nsp=as_tibble(read.table(paste0("../Table/N_species/MF/",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=","))%>%
      mutate(., V20=recode_factor(V20,"1"="Low_intra","2"="Intra=inter"))
    colnames(d_Nsp)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Stress","Rela_c","Scena")  
  
    d_Nsp=transform(d_Nsp,
                    Rela_c = factor(Rela_c, levels=c(1,2.5,4), labels=c("R[c] : 1", "R[c] : 2.5","R[c] : 4")),
                    Scena=factor(Scena,levels=c("Low_intra","Intra=inter"),labels=c("alpha[ij] < alpha[ii]","alpha[ij] > alpha[ii]")))
    
    
    assign(paste0("p_",as.numeric(com_compo)),ggplot(melt(d_Nsp,measure.vars = c(paste0("Sp_",1:Nsp)))%>%
               mutate(., Trait_sp=rep(trait_sp$V1,each=nrow(d_Nsp))))+
      geom_line(aes(x=Stress,value,color=Trait_sp,group=variable))+
      facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+
      the_theme+
      labs(x="Stress (S)",y="Densities",color=TeX("$\\psi_i$"))+
      scale_color_gradientn(colours = color_Nsp(Nsp)))
    
    p=ggplot(melt(d_Nsp,measure.vars = c(paste0("Sp_",1:Nsp)))%>%
             mutate(., Trait_sp=rep(trait_sp$V1,each=nrow(d_Nsp))))+
      geom_line(aes(x=Stress,value,color=Trait_sp,group=variable))+
      facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+
      the_theme+
      labs(x="Stress (S)",y="Densities",color=TeX("$\\psi_i$"))+
      scale_color_gradientn(colours = color_Nsp(Nsp))
    
    ggsave(paste0("../Figures/N_species/MF/",Nsp,"_species_",branch,"_",as.numeric(com_compo),"_facilitators.pdf"),p,width = 8,height = 5)
    
  }
  p_tot=ggarrange(
    p_20+ggtitle("20 % stress-tolerant")+theme(legend.position = "none")+labs(x=""),
    p_50+ggtitle("50 % stress-tolerant")+theme(legend.position = "none")+labs(x="")+
    theme(  strip.background.x = element_blank(),strip.text.x = element_blank()),
    p_80+ggtitle("80 % stress-tolerant")+theme(  strip.background.x = element_blank(),strip.text.x = element_blank())
    ,nrow = 3)
  
  ggsave(paste0("../Figures/N_species/MF/",Nsp,"_species_",branch,"_varying_facilitators.pdf"),p_tot,width = 8,height = 12)

}

### b) Species niche ----

Nsp=15
for (branch in c("Degradation","Restoration")){  
  
  d_niche=tibble()
  
  for (com_compo in paste0(c(20,50,80),".0")){
    
    trait_sp = as.numeric(read.table(paste0("../Table/N_species/MF/TRAITS_",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=",")$V1)
    d_Nsp=as_tibble(read.table(paste0("../Table/N_species/MF/",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=","))%>%
      mutate(., V20=recode_factor(V20,"1"="Low_intra","2"="Intra=inter"))
    colnames(d_Nsp)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Stress","Rela_c","Scena")  
    
    d_niche=rbind(d_niche,Get_species_niche(d=d_Nsp,Nsp=Nsp,traits=trait_sp)%>%
                    add_column(., Community_compo=as.numeric(com_compo)))
  }
  
  
  d_niche=transform(d_niche,
                    Rela_c = factor(Rela_c, levels=c(1,2.5,4), labels=c("R[c] : 1", "R[c] : 2.5","R[c] : 4")),
                    Scena=factor(Scena,levels=c("Low_intra","Intra=inter"),labels=c("alpha[ij] < alpha[ii]","alpha[ij] > alpha[ii]")))
  
  d_niche$Niche_range[!is.finite(d_niche$Niche_range)]=0
  
  d_niche_community=d_niche%>%
    group_by(Rela_c,Scena,Community_compo)%>%
    summarise(.groups = "keep",
              Mean_niche_range_sp=mean(Niche_range),Sd_niche_range_sp=sd(Niche_range),
              Mean_niche_area_sp=mean(Niche_area),Sd_niche_area_sp=sd(Niche_area))
  
  
  
  
  p1=ggplot(NULL)+
    geom_jitter(data=d_niche,aes(x=as.numeric(Community_compo),y=Niche_range,color=Trait,fill=Trait),alpha=.5,width = 0.1,height = 0,size=2)+
    # geom_line(data=d_niche,aes(x=as.numeric(Community_compo),y=Niche_range,group=Species),color="gray40",alpha=.3)+
    geom_errorbar(data=d_niche_community,aes(x=as.numeric(Community_compo),y=Mean_niche_range_sp,
                                                 ymin=Mean_niche_range_sp-Sd_niche_range_sp,
                                                 ymax=Mean_niche_range_sp+Sd_niche_range_sp), width=0,color="gray40") +
    geom_point(data=d_niche_community,aes(x=as.numeric(Community_compo),y=Mean_niche_range_sp),color="gray40",size=3,shape=0)+
    scale_color_gradientn(colours = color_Nsp(Nsp))+
    scale_fill_gradientn(colours = color_Nsp(Nsp))+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+
    labs(y="Niche range",x="Fraction of facilitators",color=TeX("$\\psi_i$"),fill=TeX("$\\psi_i$"))
  
  
  p2=ggplot(NULL)+
    geom_jitter(data=d_niche,aes(x=as.numeric(Community_compo),y=Niche_area,color=Trait,fill=Trait),alpha=.5,width = 0.1,height = 0,size=2)+
    # geom_line(data=d_niche,aes(x=as.numeric(Community_compo),y=Niche_area,group=Species),color="gray40",alpha=.3)+
    geom_errorbar(data=d_niche_community,aes(x=as.numeric(Community_compo),y=Mean_niche_area_sp,
                                             ymin=Mean_niche_area_sp-Sd_niche_area_sp,
                                             ymax=Mean_niche_area_sp+Sd_niche_area_sp), width=0,color="gray40") +
    geom_point(data=d_niche_community,aes(x=as.numeric(Community_compo),y=Mean_niche_area_sp),color="gray40",size=3,shape=0)+
    scale_color_gradientn(colours = color_Nsp(Nsp))+
    scale_fill_gradientn(colours = color_Nsp(Nsp))+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+
    labs(y="Niche area",x="Fraction of facilitative plants",color=TeX("$\\psi_i$"),fill=TeX("$\\psi_i$"))
  
  
  p_tot=ggarrange(p1+theme( legend.position = "none")+labs(x=""),
                  p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
                  nrow=2,labels = LETTERS[1:2],align = "v")
  ggsave(paste0("../Figures/N_species/MF/Species_niche_15_sp_",branch,".pdf"),p_tot,width = 8,height = 8)
}

### c) Number of shifts ----

Nsp=15
for (branch in c("Degradation","Restoration")){  
  
  d_shift=tibble()
  
  for (com_compo in paste0(c(20,50,80),".0")){
    
    trait_sp = as.numeric(read.table(paste0("../Table/N_species/MF/TRAITS_",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=",")$V1)
    d_Nsp=as_tibble(read.table(paste0("../Table/N_species/MF/",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=","))%>%
      mutate(., V20=recode_factor(V20,"1"="Low_intra","2"="Intra=inter"))
    colnames(d_Nsp)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Stress","Rela_c","Scena")  
    
    d_shift=rbind(d_shift,Get_number_shifts(d=d_Nsp,tresh = 0.01,traits=trait_sp)%>%
                    add_column(., Community_compo=as.numeric(com_compo)))
  }
  
  
  d_shift=transform(d_shift,
                    Rela_c = factor(Rela_c, levels=c(1,2.5,4), labels=c("R[c] : 1", "R[c] : 2.5","R[c] : 4")),
                    Scena=factor(Scena,levels=c("Low_intra","Intra=inter"),labels=c("alpha[ij] < alpha[ii]","alpha[ij] > alpha[ii]")))
  
  d_shift_summed=d_shift %>%
    group_by(., Scena,Community_compo,Rela_c)%>%
    summarise(.groups = "keep",Nb_shift=sum(Shift),Mean_trait=mean(Trait*Shift)*(Nsp/sum(Shift)))
  

  p=ggplot(d_shift_summed)+
    geom_bar(aes(x=Community_compo,y=Nb_shift/Nsp,fill=Mean_trait),stat = "identity",color="gray40",width = 15)+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+scale_x_continuous(breaks = c(20,50,80))+
    labs(y="Fraction of species shifting abruptly",x="Fraction of facilitative plants",fill="Mean trait of species that shift")+
    scale_fill_gradientn(colours=color_Nsp(15),breaks=c(0,.25,.5,.75,1))
  
  ggsave(paste0("../Figures/N_species/MF/Nb_shift_",branch,".pdf"),width = 6,height = 5)
  
}


### d) Hysteresis size ----


Nsp=15
d=tibble()

for (branch in c("Degradation","Restoration")){  
  
  for (com_compo in paste0(c(20,50,80),".0")){
    
    trait_sp = as.numeric(read.table(paste0("../Table/N_species/MF/TRAITS_",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=",")$V1)
    d_Nsp=as_tibble(read.table(paste0("../Table/N_species/MF/",Nsp,"_species_",branch,"_frac_facilitator_",com_compo,".csv"),sep=","))%>%
      mutate(., V20=recode_factor(V20,"1"="Low_intra","2"="Intra=inter"))
    colnames(d_Nsp)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Stress","Rela_c","Scena")  
    d=rbind(d,d_Nsp%>%add_column(., Branch=branch,Community_compo=as.numeric(com_compo)))
  }
}

d_hysteresis=Compute_hysteresis(d,Nsp=15)

d_hysteresis=transform(d_hysteresis,
                  Rela_c = factor(Rela_c, levels=c(1,2.5,4), labels=c("R[c] : 1", "R[c] : 2.5","R[c] : 4")),
                  Scena=factor(Scena,levels=c("Low_intra","Intra=inter"),labels=c("alpha[ij] < alpha[ii]","alpha[ij] > alpha[ii]")))


ggplot(d_hysteresis)+
  geom_point(aes(Community_compo,Hysteresis_area,color=Species))+
  geom_line(aes(Community_compo,Hysteresis_area,group=Species),color="gray40")+
  the_theme+
  labs(y="Fraction of species shifting abruptly",x="Fraction of facilitative plants",fill="Mean trait of species that shift")+
  scale_color_gradientn(colours=color_Nsp(15),breaks=c(0,.25,.5,.75,1))+
  facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+scale_x_continuous(breaks = c(20,50,80))

### e) Alternative stable states ----


Nsp=15
for (com_compo in c("20.0","50.0","80.0")){
  d_ASS=tibble()
  for (s in c("low_inter","nested")){
    for (rela_c in c("1.0", "2.5","4.0")){
      trait_sp=as.vector(read.table(paste0("../Table/N_species/MF/ASS_15_species/TRAITS_",
                                        Nsp,"_species_frac_facilitator_",com_compo,"_RelaC_",rela_c,
                                        "_Scena_",s,".csv"),sep=",")$V1)
        
      d_Nsp=as_tibble(read.table(paste0("../Table/N_species/MF/ASS_15_species/",
                                        Nsp,"_species_frac_facilitator_",com_compo,"_RelaC_",rela_c,
                                        "_Scena_",s,".csv"),sep=","))%>%
        mutate(., V20=recode_factor(V20,"1"="Low_intra","2"="Intra=inter"))
      colnames(d_Nsp)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Stress","Rela_c","Scena","Nrep")  
      
      for (stress in unique(d_Nsp$Stress)){
        d_stress=filter(d_Nsp,Stress==stress)
        
        #cf Liautaud et al., ELE 2019
        CSI=sapply(1:nrow(d_stress),function(x,v=Vec_densities_ass){ 
          return(sum(v*d_stress[x,1:15]))
        })
        
        Mean_trait=sapply(1:nrow(d_stress),function(x){ 
          return(sum(as.numeric(d_stress[x,1:15] * trait_sp ))/sum(as.numeric(d_stress[x,1:15]))) #mean trait of the community
        })
        
        Global_vege=sapply(1:nrow(d_stress),function(x){ 
          return(sum(d_stress[x,1:15]))
        })
        
        
        d_ASS=rbind(d_ASS,tibble(CSI=CSI,Tot_vege=Global_vege,Mean_trait=Mean_trait,
                                 Stress=stress,Rela_c=unique(d_stress$Rela_c),Scena=unique(d_stress$Scena),Community_compo=com_compo))
      }
    }
  }

  d_ASS=transform(d_ASS,
                  Rela_c = factor(Rela_c, levels=c(1,2.5,4), labels=c("R[c] : 1", "R[c] : 2.5","R[c] : 4")),
                  Scena=factor(Scena,levels=c("Low_intra","Intra=inter"),labels=c("alpha[ij] < alpha[ii]","alpha[ij] > alpha[ii]")))
  
  p1=ggplot(d_ASS)+
    geom_point(aes(x=Stress,y=CSI),alpha=.3,size=.5)+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+
    labs(y=TeX(r'(Community index \ \ $ \hat{\rho_{\psi_i}} \bullet	\hat{v}$)'),x="Stress (S)")
  
  p2=ggplot(d_ASS)+
    geom_point(aes(x=Stress,y=Tot_vege),alpha=.3,size=.5)+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+
    labs(y=TeX(r'(Global vegetation \ \ $\sum_i \rho_{\psi_i}$)'),x="Stress (S)")
  
  
  p_tot=ggarrange(p1+theme( legend.position = "none")+labs(x=""),
                  p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
                  nrow=2,labels = LETTERS[1:2],align = "v")
  ggsave(paste0("../Figures/N_species/MF/ASS_&_CSI_frac_facilitator_",com_compo,".pdf"),p_tot,width = 8,height = 8)

  
  #colored by mean trait of species
  p1=ggplot(d_ASS)+
    geom_point(aes(x=Stress,y=CSI,color=Mean_trait,fill=Mean_trait),alpha=.3,size=.5)+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+
    scale_color_gradientn(colours = color_Nsp(Nsp))+
    scale_fill_gradientn(colours = color_Nsp(Nsp))+
    labs(y=TeX(r'(Community index \ \ $ \hat{\rho_{\psi_i}} \bullet	\hat{v}$)'),x="Stress (S)")
  
  p2=ggplot(d_ASS)+
    geom_point(aes(x=Stress,y=Tot_vege,color=Mean_trait,fill=Mean_trait),alpha=.3,size=.5)+
    facet_grid(Scena~Rela_c,scales = "free",labeller = label_parsed)+the_theme+
    scale_color_gradientn(colours = color_Nsp(Nsp))+
    scale_fill_gradientn(colours = color_Nsp(Nsp))+
    labs(y=TeX(r'(Global vegetation \ \ $\sum_i \rho_{\psi_i}$)'),x="Stress (S)")
  
  p_tot=ggarrange(p1+theme( legend.position = "none")+labs(x=""),
                  p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
                  nrow=2,labels = LETTERS[1:2],align = "v")
  ggsave(paste0("../Figures/N_species/MF/ASS_&_CSI_frac_facilitator_",com_compo,"_colored_mean_trait.pdf"),p_tot,width = 8,height = 8)
  
}
