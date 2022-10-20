source('./N_species_Analysis_functions.R')

# Step 1 : N species MF -----
## 0) First test with 15 species ----

d_15sp=as_tibble(read.table("../Table/N_species/MF/15_species.csv",sep=","))
colnames(d_15sp)=c(paste0("Sp_",1:15),"Fertile","Degraded","Stress","alpha_0","h")  
d_15sp$Rho_p=rowSums(d_15sp[,1:15])


p=ggplot(melt(d_15sp,measure.vars = c(paste0("Sp_",1:15),"Rho_p")))+
  geom_point(aes(x=Stress,value,color=variable),size=.5)+
  the_theme+facet_grid(h~alpha_0,scales = "free")+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=c(color_Nsp(15),"gray50"),labels=c(paste0("Sp",1:15),"All"))

ggsave("../Figures/N_species/MF/15_species_scenario.pdf",p,width = 8,height = 5)



trait_15sp=read.table("../Table/N_species/MF/Trait_15sp.csv",sep=",")$V1
#computing diversity, functional diversity and CWM of trait psi

d_diversity=tibble()
for (i in 1:nrow(d_15sp)){ #for each simulation along the stress gradient
  
  densi=d_15sp[i,1:15]
  densi=round(densi,4)
  if (all(densi==0) | length(which(densi==0))==14){
    d_diversity=rbind(d_diversity,tibble(FD_0 =0, FD_1=0,  FD_2 =0,  D_0=0,   D_1=0,   D_2=0,)%>%add_column(., Stress=d_15sp$Stress[i],
                                                                                                            alpha_0=d_15sp$alpha_0[i]))
  } else if (any(densi==0) & length(which(densi==0))!=14 & length(which(densi==0))!=15){
    trait_s=trait_15sp[-which(densi==0)]
    densi=densi[-which(densi==0)]
    d_diversity=rbind(d_diversity,Get_diversity_community(trait =trait_s ,densities =densi )%>%add_column(., Stress=d_15sp$Stress[i],
                                                                                                          alpha_0=d_15sp$alpha_0[i]))
  } else {
    trait_s=trait_15sp
    d_diversity=rbind(d_diversity,Get_diversity_community(trait =trait_s ,densities =densi )%>%add_column(., Stress=d_15sp$Stress[i],
                                                                                                          alpha_0=d_15sp$alpha_0[i]))
  }
}


p=ggplot(d_diversity%>%
           melt(., id.vars=c("Stress","alpha_0"))%>%
           mutate(., Type_diversity=sapply(1:nrow(.), function(x){
             if (strsplit(as.character(.$variable[x]),"_")[[1]][1] == "FD"){return("Functional diversity")
              }else {return("Species diversity")}
           }))%>%
           mutate(., Hill_number=sapply(1:nrow(.), function(x){
             return(paste0("q = ",as.numeric(strsplit(as.character(.$variable[x]),"_")[[1]][2])))
           }))
         )+
  geom_line(aes(x=Stress,y=value,color=variable))+
  facet_grid(alpha_0~Type_diversity+Hill_number,scales = "free_y")+the_theme
ggsave("../Figures/N_species/MF/Test_trait_diversity.pdf",p,width = 10,height = 4)




## 1) 15 species, varying community composition ----
d_15sp=as_tibble(read.table("../Table/N_species/MF/15_species_random_ini.csv",sep=","))
colnames(d_15sp)=c(paste0("Sp_",1:15),"Fertile","Degraded","Stress","cintra","alpha_0","facilitation")  


ggplot(melt(d_15sp,measure.vars = c(paste0("Sp_",1:15))))+
  geom_point(aes(x=Stress,value,color=variable),size=.5)+
  the_theme+facet_grid(facilitation+cintra~alpha_0,labeller = label_both)+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=c(color_Nsp(15)),labels=c(paste0("Sp",1:15)))


trait_15sp=read.table("../Table/N_species/MF/Trait_15sp.csv",sep=",")$V1
#computing diversity, functional diversity and CWM of trait psi

d_diversity=tibble()
for (i in 1:nrow(d_15sp)){ #for each simulation along the stress gradient
  
  densi=d_15sp[i,1:15]
  densi=round(densi,4)
  if (all(densi==0) | length(which(densi==0))==14){
    d_diversity=rbind(d_diversity,tibble(FD_0 =0, FD_1=0,  FD_2 =0,  D_0=0,   D_1=0,   D_2=0,)%>%
                        add_column(., Stress=d_15sp$Stress[i],
                                   alpha_0=d_15sp$alpha_0[i],
                                   cintra=d_15sp$cintra[i],
                                   facilitation=d_15sp$facilitation[i]))
  } else if (any(densi==0) & length(which(densi==0))!=14 & length(which(densi==0))!=15){
    trait_s=trait_15sp[-which(densi==0)]
    densi=densi[-which(densi==0)]
    d_diversity=rbind(d_diversity,Get_diversity_community(trait =trait_s ,densities =densi )%>%
                        add_column(., Stress=d_15sp$Stress[i],
                                    alpha_0=d_15sp$alpha_0[i],
                                    cintra=d_15sp$cintra[i],
                                    facilitation=d_15sp$facilitation[i]))
  } else {
    trait_s=trait_15sp
    d_diversity=rbind(d_diversity,Get_diversity_community(trait =trait_s ,densities =densi )%>%
                        add_column(., Stress=d_15sp$Stress[i],
                                    alpha_0=d_15sp$alpha_0[i],
                                    cintra=d_15sp$cintra[i],
                                    facilitation=d_15sp$facilitation[i]))
  }
}



ggplot(d_diversity%>%
         filter(., facilitation==.9)%>%
         rename(.,f=facilitation,aii=cintra,a0=alpha_0 )%>%
           melt(., id.vars=c("Stress","a0","f","aii")))+
  geom_line(aes(x=Stress,y=value,color=variable))+
  facet_grid(a0+aii~variable,scales = "free_y",labeller = label_both)+the_theme

