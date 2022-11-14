source('./N_species_Analysis_functions.R')
library(FactoMineR) 
library(factoextra)
library(cluster)

## 1) first exploration with CSI & rho_+ ----
d=tibble()
Nsp=15
list_csv=list.files('../Table/N_species/MF/Big_sim_Nsp/Random_ini/')

for ( k in 1:length(list_csv)){
  d2=read.table(paste0("../Table/N_species/MF/Big_sim_Nsp/Random_ini/",list_csv[k]),sep=",")
  colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
  d=rbind(d,d2)
}
d$Rho_p=rowSums(d[,1:Nsp])

d[d<10^(-4)]=0

d$CSI = sapply(1:nrow(d),function(x){
  set.seed(432)
  u=runif(Nsp)
  return(sum(d[x,1:Nsp]*u))})

d$Psi_normalized = sapply(1:nrow(d),function(x){
  trait=rev(seq(0,1,length.out=Nsp))
  if (d$CSI[x]==0){return(NA)
  }else{  return(sum(d[x,1:Nsp]*trait)/rowSums(d[x,1:Nsp]))
  }
})

d$Branch=sapply(1:nrow(d),function(x){
  return(ifelse(d$Branch[x]==1,"Degradation","Restoration"))
})

write.table(d,'../Table/N_species/MF/Random_ini.csv',sep=";")


d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")
Nsp=15

pD=ggplot(d%>%filter(., Branch=="Degradation"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+ggtitle("Degradation")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")

pR=ggplot(d%>%filter(., Branch=="Restoration"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+ggtitle("Restoration")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")
ggsave("../Figures/N_species/MF/CSI_random_ini.pdf",ggarrange(pD,pR,nrow = 2),width = 7,height = 10)

pD=ggplot(d%>%filter(., Branch=="Degradation"))+
  geom_line(aes(x=Stress,y=Rho_p,color=Psi_normalized,group=Random_ini),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Species density",color=expression(paste(bar(psi),"    ")))+ggtitle("Degradation")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")

pR=ggplot(d%>%filter(., Branch=="Restoration"))+
  geom_line(aes(x=Stress,y=Rho_p,color=Psi_normalized,group=Random_ini),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Species density",color=expression(paste(bar(psi),"    ")))+ggtitle("Restoration")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")
ggsave("../Figures/N_species/MF/Global_cover_random_ini.pdf",ggarrange(pD,pR,nrow = 2),width = 7,height = 10)




## 2) Species level : occurrence, size of tipping points  ----

d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")
treshold_sp=0.1
d_ASS_sp=tibble()

for (i in c(5,10,15,20,25,30,35)){
  
  d_t=tibble()
  
  Nsp=i
  traits=rev(seq(0,1,length.out=Nsp))
  list_csv=list.files(paste0('../Table/N_species/MF/',i,'_sp/'))
  
  for ( k in 1:length(list_csv)){ #for each replicate
    
    d2=read.table(paste0("../Table/N_species/MF/",i,"_sp/",list_csv[k]),sep=",")
    colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
    
    d2[d2<10^(-4)]=0
    d2=melt(d2,measure.vars=paste0("Sp_",1:Nsp))
    
    for (sp in 1:Nsp){
      
      d_sp=filter(d2,variable==paste0("Sp_",sp)) 
      
      
      if (any(d_sp$value[1:(nrow(d_sp)/2)]>0)){# we only keep species that are present in the system
        
        tipping=ifelse(any(abs(diff(d_sp$value[1:(nrow(d_sp)/2)]))>treshold_sp),1,0) #Is there a tipping point in the system for species sp
        
        d_ASS_sp=rbind(d_ASS_sp,tibble(Tipping=tipping,Species=sp,
                                       Trait=traits[sp],Branch="Degradation",
                                       Random_ini=k,Richness=Nsp))
        
      }
      
      if (any(d_sp$value[((nrow(d_sp)/2)+1):nrow(d_sp)] >0)){# we only keep species that are present in the system
        
        tipping=ifelse(any(abs(diff(d_sp$value[((nrow(d_sp)/2)+1):nrow(d_sp)]))>treshold_sp),1,0) #Is there a tipping point in the system for species sp
        
        d_ASS_sp=rbind(d_ASS_sp,tibble(Tipping=tipping,Species=sp,
                                       Trait=traits[sp],Branch="Restoration",
                                       Random_ini=k,Richness=Nsp))
        
      }
    }
    
  }
  
}

d_ASS_sp=d_ASS_sp[-which(d_ASS_sp$Random_ini>100),]
write.table(d_ASS_sp,"../Table/N_species/MF/d_ASS_sp.csv",sep=";")

d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=100


#Species-specific tipping points size & frequency
p=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .1)+
  facet_grid(Branch~Richness)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Fraction of tipping points")+
  scale_fill_gradientn(colours = color_Nsp(100))

ggsave("../Figures/N_species/MF/Fraction_tipping_points.pdf",p,width = 9,height = 5)





## 3) varying species diversity ----

threshold=.1
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

write.table(d_tipping,"../Table/N_species/MF/Multistability_tipping.csv",sep=";")
write.table(d_richness,"../Table/N_species/MF/Multistability_richness.csv",sep=";")
write.table(d_richness2,"../Table/N_species/MF/Multistability_richness2.csv",sep=";")


p1=ggplot(d_tipping%>%filter(., Branch=="Degradation"))+
  geom_bar(aes(x=Nsp,fill=as.factor(Nb_different_species),y=Nb_different_species),
           position = "fill",stat="identity",alpha=.75)+
  labs(x="Number of species",y="Fraction of random initial conditions",fill="Number of shifts")+
  scale_x_continuous(breaks = seq(5,35,by=5))+
  the_theme

p2=ggplot(d_richness%>%filter(., Branch=="Degradation"))+
  geom_bar(aes(x=Nsp,fill=as.factor(Richness),y=Richness),
           position = "fill",stat="identity",alpha=.75)+
  geom_point(data=d_richness2,aes(x=Nsp,y=Richness_D/Nsp),shape=0,size=2)+
  geom_line(data=d_richness2,aes(x=Nsp,y=Richness_D/Nsp))+
  scale_y_continuous(sec.axis=sec_axis(~., name="Fraction of species observed \n across replicates"))+
  labs(x="Number of species",y="Fraction of random initial conditions",fill="Number of species")+
  scale_x_continuous(breaks = seq(5,35,by=5))+
  the_theme


ggsave("../Figures/N_species/MF/Nb_shift_diversity.pdf",ggarrange(p1,p2,ncol = 2,labels = letters[1:2],align = "hv"),width=10,height=5)
