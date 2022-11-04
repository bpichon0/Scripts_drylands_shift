source('./N_species_Analysis_functions.R')
library(FactoMineR) 
library(factoextra)
library(cluster)

# Step 1 : N species MF -----
## 0) First test with 15 species ----

d=as_tibble(read.table("../Table/N_species/MF/15_species.csv",sep=","))
colnames(d)=c(paste0("Sp_",1:15),"Fertile","Degraded","Stress","Competition","Community_compo","Random_ini","Branch")  
d$Rho_p=rowSums(d[,1:15])
d$CSI = sapply(1:nrow(d),function(x){
  set.seed(123)
  u=runif(15)
  return(sum(d[x,1:15]*u))
}
)
d$Facilitation=.9
d$Disp=.1


p=ggplot(melt(d,measure.vars = c("Rho_p")))+
  geom_line(aes(x=Stress,value,group=interaction(variable,Branch,Random_ini),color=as.factor(Branch)),size=.5)+
  the_theme+facet_wrap(.~Branch)+
  labs(x="Stress (S)",y="Densities",color="")

ggsave("../Figures/N_species/MF/15_species_scenario_community.pdf",p,width = 8,height = 5)



p=ggplot(melt(d,measure.vars = c(paste0("Sp_",1:15)))%>%
           mutate(., Branch=recode_factor(Branch,"1"="Restoration","2"="Degradation")))+
  geom_line(aes(x=Stress,value,color=variable,group=interaction(variable,Branch,Random_ini)),size=.5)+
  the_theme+facet_grid(Branch~Competition,scales = "free")+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=c(color_Nsp(15),"gray50"),labels=c(paste0("Sp",1:15),"All"))

ggsave("../Figures/N_species/MF/15_species_scenario_species.pdf",p,width = 8,height = 5)


p=ggplot(melt(d,measure.vars = "CSI")%>%
           mutate(., Branch=recode_factor(Branch,"1"="Restoration","2"="Degradation")))+
  geom_point(aes(x=Stress,value,color=as.factor(Random_ini)),size=1,alpha=.7,shape=1)+
  #geom_line(aes(x=Stress,value,group=interaction(Random_ini,Branch)))+
  the_theme+facet_grid(Branch~Competition,scales = "free")+
  labs(x="Stress (S)",y="Community index",color="Random initiam condition")
ggsave("../Figures/N_species/MF/CSI_15_species.pdf",p,width = 8,height = 5)

## 1) Varying facilitation and competition ----


Nsp=15
d = read.table("../Table/N_species/MF/Varying_competition_facilitation.csv",sep=",")
colnames(d) = c(paste0("Rho_",1:Nsp),"Rho_0","Rho_d","Stress","Competition","Facilitation")
d$Rho_p=rowSums(d[,1:Nsp])
d$CSI = sapply(1:nrow(d),function(x){
  set.seed(123)
  u=runif(Nsp)
  return(sum(d[x,1:Nsp]*u))
}
)

ggplot(d%>%melt(., measure.vars=c("CSI")))+
  geom_point(aes(x=Stress,y=value))+
  the_theme+
  facet_grid(Competition~Facilitation)

ggplot(d%>%melt(., measure.vars=paste0("Rho_",1:Nsp)))+
  geom_line(aes(x=Stress,y=value,group=interaction(Facilitation,Competition,variable)))+
  the_theme+
  facet_grid(Competition~Facilitation)

## 2) Number tipping points, hysteresis----
Nsp=15
d=read.table("../Table/N_species/MF/Number_tipping.csv",sep=",")
colnames(d)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Stress","Competition","Facilitation","Random_ini")

d$Rho_p=rowSums(d[,1:Nsp])
d$CSI = sapply(1:nrow(d),function(x){
  set.seed(123)
  u=runif(Nsp)
  return(sum(d[x,1:Nsp]*u))})
pdf("../Figures/N_species/MF/Number_shifts.pdf",width = 8,height = 5)

print(ggplot(d%>%melt(., measure.vars=paste0("Sp_",1:15)))+
        geom_line(aes(x=Stress,y=value,group=interaction(Facilitation,Competition,Random_ini,variable)))+
        the_theme+
        facet_grid(Competition~Facilitation))

print(ggplot(d)+
        geom_line(aes(x=Stress,y=Rho_p,group=interaction(Facilitation,Competition,Random_ini)))+
        the_theme+
        facet_grid(Competition~Facilitation))

print(ggplot(d)+
        geom_point(aes(x=Stress,y=CSI,group=interaction(Facilitation,Competition,Random_ini)),size=.6,alpha=.3)+
        the_theme+
        facet_grid(Competition~Facilitation))
dev.off()

#For testing 
d=tibble()
for (i in 1:5) {
  Nsp=15
  d2=read.table(paste0("../Table/N_species/MF/Big_sim_NSp/Sim_Nrandom_",i,"_a0_0.3_delta_0.1_f_0.9.csv"),sep=",")
  colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
  d2$Branch=sapply(1:nrow(d2),function(x){
    return(ifelse(d2$Branch[x]==1,"Degradation","Restoration"))
  })
  d2$Rho_p=rowSums(d2[,1:Nsp])
  d2$CSI = sapply(1:nrow(d2),function(x){
    set.seed(123)
    u=runif(Nsp)
    return(sum(d2[x,1:Nsp]*u))})
  
  d=rbind(d,d2)
  
}



## 3) Hysteresis size ----

d=tibble()
Nsp=15
list_csv=list.files('../Table/N_species/MF/Big_sim_Nsp/Constant_ini/')
for ( k in 1:length(list_csv)){
  d2=read.table(paste0("../Table/N_species/MF/Big_sim_Nsp/Constant_ini/",list_csv[k]),sep=",")
  colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
  d=rbind(d,d2)
}
d$Rho_p=rowSums(d[,1:Nsp])
d$CSI = sapply(1:nrow(d),function(x){
  set.seed(123)
  u=runif(Nsp)
  return(sum(d[x,1:Nsp]*u))})
d$Branch=sapply(1:nrow(d),function(x){
  return(ifelse(d$Branch[x]==1,"Degradation","Restoration"))
})


p=ggplot(d%>%melt(., measure.vars=paste0("Sp_",1:15))%>%
           mutate(., variable=rep(seq(0,1,length.out=Nsp),each=1800)))+
  geom_line(aes(x=Stress,y=value,group=interaction(Facilitation,Competition,Random_ini,variable,Branch),
                color=variable,linetype=Branch))+
  the_theme+labs(x="Stress, S",y="Cover",color=TeX("$\\psi$ \ \ \ "),linetype="")+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  scale_color_gradientn(colours = color_Nsp(15))
ggsave("../Figures/N_species/MF/Equal_ini.pdf",p,width = 7,height = 6)









## 4) Varying initial conditions ----
### a) first exploration with CSI & rho_+ ----
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
  set.seed(123)
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





p=ggplot(d)+
  geom_line(aes(x=Stress,y=CSI,color=Branch,group=interaction(Branch,Random_ini)),size=1)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+ggtitle("Restoration")+
  scale_color_manual(values = c("red","blue"),na.value = "black")

ggsave("../Figures/N_species/MF/Global_cover_random_ini.pdf",ggarrange(pD,pR,nrow = 2),width = 7,height = 10)
















### b) Species level : occurrence, size of tipping points and hysteresis size ----
d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")

# At species level : mean number of shifts/replicate
N_replicate=40;Nsp=15;treshold_sp=0.05;traits=rev(seq(0,1,length.out=15));treshold_com=.1


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
        
        d_fil=filter(d,Facilitation==f,Competition==a0,Branch==branch,Random_ini==ini)%>%
          melt(.,  measure.vars="Rho_p")
        
        d_fil$value[d_fil$value<10^(-4)]=0
        
        
        tipping=ifelse(any(abs(diff(d_fil$value))>treshold_com),1,0)
        
        if (tipping==1){ #there is shift one or multiple
          for (nb_shift in 1:length(which(abs(diff(d_fil$value))>treshold_com))){ #for each community shift, we get its size, and save into the DF
            size_shift = abs(diff(d_fil$value))[which(abs(diff(d_fil$value))>treshold_com)[nb_shift]]
            d_ASS_com=rbind(d_ASS_com,tibble(Tipping=tipping,Facilitation=f,Competition=a0,Branch=branch,Random_ini=ini,Size_shift=size_shift))
          }
          # 
          # #We track species that are not detected at the community level
          # stress_tipping_com=d_fil$Stress[which(abs(diff(d_fil$value))>treshold_com)]
          # 
          # if (length(stress_tipping_com) < nrow(d_shifts) |  length(unique(d_shifts$Stress))!=1 & unique(d_shifts$Stress)!=stress_tipping_com) { #i.e. we detect less shifts at community scale compared to species scale or all species shifted at the same threshold
          #   (1:nrow(d_shifts))[-which(d_shifts$Stress %in% stress_tipping_com)]
          # }
          # 
          
          
        } else{ #no shifts
          d_ASS_com=rbind(d_ASS_com,tibble(Tipping=tipping,Facilitation=f,Competition=a0,Branch=branch,Random_ini=ini,Size_shift=NA))
        }
        
        
        
      }
    }
  }
}
write.table(d_ASS_sp,"../Table/N_species/MF/d_ASS_sp.csv",sep=";")
write.table(d_ASS_com,"../Table/N_species/MF/d_ASS_com.csv",sep=";")

d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=40


#Species-specific tipping points size & frequency
p=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/(N_replicate*3),fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f == .(Facilitation)))+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Fraction of tipping points",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))

ggsave("../Figures/N_species/MF/Species_level/Fraction_tipping_points.pdf",p,width = 6,height = 5)


d_ASS_sp=d_ASS_sp%>%
  group_by(.,Facilitation,Competition,Branch,Trait)%>%
  summarise(.,Size_shift=mean(Size_shift),.groups = "keep")

p=ggplot(d_ASS_sp, aes(x=Trait,y=Size_shift,fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f == .(Facilitation)))+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Mean shift size (in fraction of cover)",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))

ggsave("../Figures/N_species/MF/Species_level/Size_tipping_points.pdf",p,width = 6,height = 5)


# Finally, we compute the species specific hysteresis size.
d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")

d_hysteresis=tibble()
Nsp=15

for (f in unique(d$Facilitation)){
  for (a0 in unique(d$Competition)){
    for (ini in unique(d$Random_ini)){
      for (sp in 1:Nsp){
        
        d_h=filter(d%>%melt(., measure.vars=paste0("Sp_",1:15)),Facilitation==f,Competition==a0,Random_ini==ini,variable==paste0("Sp_",sp))
        
        if (any(d_h$value>0)){ #if species is present
          
          d_hysteresis=rbind(d_hysteresis,tibble(Facilitation=f,
                                                 Species=sp,
                                                 Competition=a0,
                                                 Random_ini=ini,
                                                 #We scale the hysteresis size with the mean cover to account for the effect of competition on global cover.
                                                 Hysteresis_scaled=abs(sum(d_h$value[1:50]-rev(d_h$value[51:100])))/((sum(d_h$value[1:50]+rev(d_h$value[51:100])))/2),
                                                 Hysteresis_not_scaled=abs(sum(d_h$value[1:50]-rev(d_h$value[51:100]))),
                                                 
                                                 
                                                 ))
          
        }
        # plot(d_h$Stress[1:50],d_h$Rho_p[1:50],col="red",type="b",main=paste0("f=",f,", a0=",a0,", ini=",ini))
        # points(rev(d_h$Stress[51:100]),rev(d_h$Rho_p[51:100]),col="blue",type="b")
      }
    }
  }
}

d_hysteresis$Trait=sapply(1:nrow(d_hysteresis),function(x){
  return(rev(seq(0,1,length.out=Nsp))[d_hysteresis$Species[x]])
})

write.table(d_hysteresis,"../Table/N_species/MF/Hysteresis_species.csv",sep=";")



d_hysteresis=read.table("../Table/N_species/MF/Hysteresis_species.csv",sep=";")




p=ggplot(NULL)+
  geom_point(data=d_hysteresis,aes(x=Trait,y=Hysteresis_not_scaled,color=as.factor(Competition)),size=1,alpha=.5)+
  geom_smooth(data=d_hysteresis,aes(x=Trait,y=Hysteresis_not_scaled,color=as.factor(Competition)),size=1,alpha=.5,se = F)+
  the_theme+labs(y="Hysteresis size",color=TeX("$\\alpha_e$"),x=TeX("$\\psi$"),fill=TeX("$\\alpha_e$"))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  ylim(0,.75)+ #to avoid outliers
  facet_wrap(.~Facilitation,labeller = label_bquote(cols=f==.(Facilitation)))+
  scale_x_continuous(breaks = c(0,.25,.5,.75,1))


ggsave("../Figures/N_species/MF/Species_level/Hystersis_size_not_scaled.pdf",p,width = 6,height = 3)


p=ggplot(NULL)+
  geom_point(data=d_hysteresis,aes(x=Trait,y=Hysteresis_scaled,color=as.factor(Competition)),size=1,width = .1,alpha=.5)+
  geom_smooth(data=d_hysteresis,aes(x=Trait,y=Hysteresis_scaled,color=as.factor(Competition)),size=1,alpha=.5,se = F)+
  the_theme+labs(y="Hysteresis size",color=TeX("$\\alpha_e$"),x=TeX("$\\psi$"),fill=TeX("$\\alpha_e$"))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  ylim(0,.75)+ #to avoid outliers
  facet_wrap(.~Facilitation,labeller = label_bquote(cols=f==.(Facilitation)))+
  scale_x_continuous(breaks = c(0,.25,.5,.75,1))


ggsave("../Figures/N_species/MF/Species_level/Hystersis_size_scaled.pdf",p,width = 6,height = 3)

### c) Number of shifts species versus community----

d_ASS_com=read.table("../Table/N_species/MF/d_ASS_com.csv",sep=";")
d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")

d_sp=d_ASS_sp%>%
  group_by(.,Facilitation,Competition,Branch,Random_ini)%>%
  summarise(.,Tipping=sum(Tipping),.groups = "keep")

d_com=d_ASS_com%>%
  group_by(.,Facilitation,Competition,Branch,Random_ini)%>%
  summarise(.,Tipping=sum(Tipping),.groups = "keep")

p1=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/(N_replicate*3),fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f == .(Facilitation)))+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Fraction of tipping points",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))





d_all=d_com%>%
  rename(., Tipping_com=Tipping)%>%
  add_column(., Tipping_sp=d_sp$Tipping)

p2=ggplot(d_all)+
  geom_jitter(aes(x=Tipping_sp,y=Tipping_com,fill=as.factor(Competition),color=as.factor(Competition),shape=as.factor(Facilitation)),
              alpha=.7,size=.5,width=0.05, height=0.05)+
  geom_abline(slope = 1,intercept = 0,linetype=9)+
  facet_wrap(.~Branch,scales = "free")+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_shape_manual(values=c(21,22,24))+
  the_theme+
  labs(fill=TeX("$\\alpha_e$"),color=TeX("$\\alpha_e$"),shape="f",
       x="# Species scale tipping point",
       y="# Community scale tipping point")

ggsave("../Figures/N_species/MF/Nb_tipping_community_versus_species.pdf",p,width = 8,height = 4)






### d) Community level : Hysteresis size ----


d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")
Nsp=15

#The problem with this definition, is that we are not estimating well systems in which there is multiple shifts...
# pdf("../Figures/N_species/MF/Community_level/Hystersis_patterns.pdf",width = 5,height = 4)
d_hysteresis=tibble()
for (f in unique(d$Facilitation)){
  for (a0 in unique(d$Competition)){
    for (ini in unique(d$Random_ini)){
      
      d_h=filter(d,Facilitation==f,Competition==a0,Random_ini==ini)
      # plot(d_h$Stress[1:50],d_h$Rho_p[1:50],col="red",type="b",main=paste0("f=",f,", a0=",a0,", ini=",ini))
      # points(rev(d_h$Stress[51:100]),rev(d_h$Rho_p[51:100]),col="blue",type="b")
      
      d_hysteresis=rbind(d_hysteresis,tibble(Facilitation=f,
                                             Competition=a0,
                                             Random_ini=ini,
                                             #We scale the hysteresis size with the mean cover to account for the effect of competition on global cover.
                                             Hysteresis_scaled=abs(sum(d_h$Rho_p[1:50]-rev(d_h$Rho_p[51:100])))/((sum(d_h$Rho_p[1:50]+rev(d_h$Rho_p[51:100])))/2),
                                             Hysteresis_not_scaled=abs(sum(d_h$Rho_p[1:50]-rev(d_h$Rho_p[51:100]))))
      
      )
      
    }
  }
}

# dev.off()


d_hysteresis_summary = d_hysteresis%>%
  group_by(., Competition,Facilitation)%>%
  summarise(.,.groups = "keep",Mean_hystersis=mean(Hysteresis_not_scaled),Sd_hysteresis=sd(Hysteresis_not_scaled))

p=ggplot(NULL)+
  geom_jitter(data=d_hysteresis,aes(x=as.factor(Facilitation),y=Hysteresis_not_scaled,color=as.factor(Competition)),size=1,width = .1,alpha=.5)+
  geom_pointrange(data=d_hysteresis_summary,aes(x=as.factor(Facilitation),
                                                y=Mean_hystersis,
                                                ymin=Mean_hystersis-Sd_hysteresis,
                                                ymax=Mean_hystersis+Sd_hysteresis,
                                                fill=as.factor(Competition)),size=.5,color="black",shape=21)+
  the_theme+labs(y="Hysteresis size",color=TeX("$\\alpha_e$"),x="Facilitation, f",fill=TeX("$\\alpha_e$"))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))

ggsave("../Figures/N_species/MF/Community_level/Hysteresis_size_community_not_scaled.pdf",p,width = 6,height =4)



d_hysteresis_summary = d_hysteresis%>%
  group_by(., Competition,Facilitation)%>%
  summarise(.,.groups = "keep",Mean_hystersis=mean(Hysteresis_scaled),Sd_hysteresis=sd(Hysteresis_scaled))

p=ggplot(NULL)+
  geom_jitter(data=d_hysteresis,aes(x=as.factor(Facilitation),y=Hysteresis_scaled,color=as.factor(Competition)),size=1,width = .1,alpha=.5)+
  geom_pointrange(data=d_hysteresis_summary,aes(x=as.factor(Facilitation),
                                                y=Mean_hystersis,
                                                ymin=Mean_hystersis-Sd_hysteresis,
                                                ymax=Mean_hystersis+Sd_hysteresis,
                                                fill=as.factor(Competition)),size=.5,color="black",shape=21)+
  the_theme+labs(y="Hysteresis size",color=TeX("$\\alpha_e$"),x="Facilitation, f",fill=TeX("$\\alpha_e$"))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))

ggsave("../Figures/N_species/MF/Community_level/Hysteresis_size_community_scaled.pdf",p,width = 6,height =4)


### e) Idea 01/11 : at community level push the CSI distribution on PCA and identify ASS -> to discuss ----
d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")

pdf("../Figures/N_species/MF/Community_level/Detecting_ASS.pdf",width = 6,height = 4)

N_replicate=40

d_ASS=tibble()
for (f in unique(d$Facilitation)){
  for (a0 in unique(d$Competition)){
    for (branch in unique(d$Branch)){
      
      d_fil=filter(d,Facilitation==f,Competition==a0,Branch==branch)
      for_pca=cbind(tibble(Random_ini=1:N_replicate),as_tibble(matrix(d_fil$CSI,N_replicate,50,byrow = T)))
      for_pca=round(for_pca,2)
      
      res.pca=PCA(for_pca[,2:ncol(for_pca)], scale.unit = F, ncp = 3,  graph=F)
      print(ggplot(d_fil)+
              geom_point(aes(x=Stress,y=CSI,group=Random_ini))+
              the_theme+
              ggtitle(paste0("f = ",f,", alpha_e = ",a0,", branch = ",branch)))
      
      print(fviz_pca_ind(res.pca, geom.ind = "point", 
                         axes=c(1,2), 
                         pointsize = 2, 
                         label = "var",
                         repel = T,mean.point = FALSE)+
              labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
                   y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),col.ind="")+
              ggtitle("") +theme_classic()+theme(legend.position = "bottom",
                                                 legend.text = element_text(size=14)))
      
      #dendrogram
      res.hcpc=HCPC(res.pca, graph = FALSE)
      print(fviz_cluster(res.hcpc)+
              the_theme+
              ggtitle(""))
      
      
      d_ASS=rbind(d_ASS,tibble(Nb_ASS=length(res.hcpc$desc.ind$dist),Competition=a0,Facilitation=f))
    }
  }
}
dev.off()


#Same but merging the two branches

N_replicate=40

d_ASS=tibble()
for (f in unique(d$Facilitation)){
  for (a0 in unique(d$Competition)){
    
    d_fil=filter(d,Facilitation==f,Competition==a0)
    for_pca=cbind(tibble(Random_ini=1:N_replicate),as_tibble(matrix(d_fil$CSI,N_replicate,50,byrow = T)))
    for_pca=round(for_pca,2)
    
    res.pca=PCA(for_pca[,2:ncol(for_pca)], scale.unit = F, ncp = 3,  graph=F)
    print(ggplot(d_fil)+
            geom_point(aes(x=Stress,y=CSI,group=interaction(Random_ini),color=Branch))+
            the_theme+
            ggtitle(paste0("f = ",f,", alpha_e = ",a0)))
    
    print(fviz_pca_ind(res.pca, geom.ind = "point", 
                       axes=c(1,2), 
                       pointsize = 2, 
                       label = "var",
                       repel = T,mean.point = FALSE)+
            labs(x=paste0("PC 1 (",round(res.pca$eig[1,2], 1)," %)"),
                 y=paste0("PC 2 (",round(res.pca$eig[2,2], 1)," %)"),col.ind="")+
            ggtitle("") +theme_classic()+theme(legend.position = "bottom",
                                               legend.text = element_text(size=14)))
    
    #dendrogram
    res.hcpc=HCPC(res.pca, graph = FALSE)
    print(fviz_cluster(res.hcpc)+the_theme+ggtitle(""))
    
    
    d_ASS=rbind(d_ASS,tibble(Nb_ASS=length(res.hcpc$desc.ind$dist),Competition=a0,Facilitation=f))
    
  }
}



