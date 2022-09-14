# Other) Recruitment rate ----
rm(list=ls())
library(tidyverse)
d=expand_grid(emax=1,e=.5,S=seq(0,1,length.out=200),psi=c(0,.5,1))

p=ggplot(d)+geom_line(aes(x=S,y=emax*(1-S*(1-e*psi)),color=as.factor(psi),group=psi))+
scale_color_manual(values=c("blue","green","red"))+theme_classic()+
theme(legend.position = "bottom")+labs(x="S",
y=TeX("$\\epsilon_{max} (1-S (1-e \\psi))$"),color=TeX("\\psi"))



ggsave("../Figures/Recruitment_rate.pdf",p,width = 6,height =4 )



# Step 1) Two species MF model ----
x=c("tidyverse","ggpubr","latex2exp","deSolve","reshape2","JuliaCall","diffeqr")
lapply(x, require, character.only = TRUE)
julia_setup()
de = diffeq_setup()


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))


MF_two_species_julia = julia_eval("

function MF_two_species(du,u,p,t)

  r,d,f,beta,m,e,emax,cintra,cinter1,cinter2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) + rho_2 * m + rho_1 * m

end")

plot_dynamics=function(d){
  
  color_rho=c("fertile"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  if ('time' %in% colnames(d) | 'Time' %in% colnames(d)){
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)])%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive","rho_0"="fertile","rho_d"="desert")))+
             geom_line(aes(x=time,y=value,color=variable),lwd=1)+
             the_theme+scale_color_manual(values=color_rho)+labs(x="Time",y="Densities",color=""))
  }
  else {
    d$time=1:nrow(d)
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)])%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive","rho_0"="fertile","rho_d"="desert")))+
             geom_line(aes(x=time,y=value,color=variable),lwd=1)+
             the_theme+scale_color_manual(values=color_rho)+labs(x="Time",y="Densities",color=""))
  }
}



## 1) No difference of competition, influence of tolerance ----
 

param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=0.1,emax=1.2,c=0.3,S=0)
state = c(rho_1=0.4,rho_2=0.4,rho_d=0.1,rho_0=0.1)
tspan = c(0, 10000)
t = seq(0, 10000, by=1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)

MF_two_species_julia = julia_eval( "
function MF_two_species(du,u,p,t)
  r,d,f,beta,m,e,emax,c,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c*(rho_2 + rho_1))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c*(rho_2 + rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c*(rho_2 + rho_1)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c*(rho_2 + rho_1))) + rho_2 * m + rho_1 * m


end")




d2=tibble();S_seq= seq(0,1,length.out=100)
for (S in S_seq){
  for (e in c(.1,.4)){
    
    param[9]=S;param[6]=e
    julia_assign("p", param)
    prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
    sol = de$solve(prob,de$Tsit5(),saveat=t)
    d=as.data.frame(t(sapply(sol$u,identity)))
    d2=rbind(d2,d[nrow(d),]%>%add_column(e=e,S=S))
    
  }
}
colnames(d2)=c("rho_1","rho_2",'rho_d',"rho_0","e","S")
d2$rho_plus=d2$rho_1+d2$rho_2

p1=ggplot(d2%>%melt(., id.vars=c("S","e")))+geom_point(aes(x=S,y=value,color=as.factor(variable),group=variable),size=1)+
  theme_classic()+
  theme(legend.position = "bottom")+labs(x="Stress",y="Density")+
  scale_color_manual(name="",values=color_rho,
                     labels=c(TeX("$\\rho_1$"), TeX("$\\rho_2$"),TeX("$\\rho_d$"),TeX("$\\rho_0$"),TeX("$\\rho_+$")))+
  facet_wrap(.~e,labeller = label_both)+
  theme(legend.text = element_text(size=14))

ggsave("../Figures/2_species/Bifu_kefi_2_species_no_compet.pdf",p1,width = 7,height = 4)
















## 2) State diagram with competitive ability and stress for different type of competition ----

#We now try to replicate Danet figure with abiotic stress in x-axis and differential competitive ability in y-axis




#Intra + interspecific competition

MF_two_species_julia = julia_eval("

function MF_two_species(du,u,p,t)

  r,d,f,beta,m,e,emax,cintra,cinter1,cinter2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) + rho_2 * m + rho_1 * m

end")



Post_processing_MF=function(d2,name,C_seq=c_seq){
  
  write.table(d2,paste0("../Table/2_species/2_species_",name,".csv"),sep=";")
  
  
  d2$state=sapply(1:nrow(d2),function(x){
    
    if (d2[x,1]>0 & d2[x,2]>0) return("coexistence")
    if (d2[x,1]>0 & d2[x,2]==0) return("stress_tol")
    if (d2[x,1]==0 & d2[x,2]>0) return("competitive")
    if (d2[x,1]==0 & d2[x,2]==0) return('desert')
    
  })
  
  c_values_bifu=C_seq[c(round(length(C_seq)/10),round(length(C_seq)/2),length(C_seq))]
  
  
  color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  #state at equilibrium
  p1=ggplot(d2)+geom_tile(aes(x=S,y=cinter1/param["cinter2"],fill=state))+
    theme_classic()+scale_fill_manual(values=color_rho)+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="")+
    geom_hline(yintercept = c_values_bifu/param["cinter2"],lwd=.1,color="gray40")+
    theme(legend.text = element_text(size=11))
  
  
  #density of global vegetation
  density_col=colorRampPalette(c("red","white","blue"))
  p2= ggplot(d2)+geom_tile(aes(x=S,y=cinter1/param["cinter2"],fill=rho_plus))+
    theme_classic()+scale_fill_gradientn(colours=density_col(100))+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="Global vegetation density")
  
  
  #some bifurcation diagrams
  
  c_values_bifu=c_seq[c(1,round(length(c_seq)/2),length(c_seq))]
  d_bifu=filter(d2,cinter1 %in% c_values_bifu)
  for (c_bifu in 1:3){
    assign(paste0("p3_",c_bifu),
           ggplot(d_bifu %>% filter(.,cinter1==c_values_bifu[c_bifu]) %>% melt(.,measure.vars=c("rho_1","rho_2"))%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive")))+
             geom_point(aes(x=S,y=value,color=variable))+labs(x="Stress (S)",y="Density")+
             the_theme+ scale_color_manual(values=color_rho[c(2,4)])+  theme(legend.text = element_text(size=11))
    )
  }
  p3=ggarrange(p3_3,p3_2,p3_1,nrow=3,common.legend = T,legend = "bottom")
  
  
  p_tile=ggarrange(p1,p2,ncol=2)
  p_tot=ggarrange(p_tile,p3,ncol=2,widths = c(3,1))
  ggsave(paste0("../Figures/2_species/2_species_",name,".pdf"),width = 14,height = 6)
  
  
  
  
}


#preparing the 3 types of simulations


state = c(rho_1=0.4,rho_2=0.4,rho_d=0.1,rho_0=0.1)
tspan = c(0, 5000)
t = seq(0, 5000, by=1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)


for_sim=tibble(scena=c("inter=0","intra>inter","inter=intra","intra=0"),
               cinter2=c(0,.05,.1,.15),
               cintra=c(.25,.2,.1,0))


for (i in 1:nrow(for_sim)){
  
  param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=.1,emax=1.2,cintra=for_sim$cintra[i],
          cinter1=for_sim$cinter2[i],cinter2=for_sim$cinter2[i],S=0)
  
  
  d2=tibble();S_seq= seq(0,1,length.out=100);c_seq=seq(param["cinter1"],4*param["cinter1"],length.out=100)
  for (S in S_seq){
    for (c1 in c_seq){
      
      param[11]=S;param[9]=c1
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob,de$Tsit5(),saveat=t)
      d=as.data.frame(t(sapply(sol$u,identity)))
      colnames(d)=c("rho_1","rho_2","rho_d","rho_0")
      d2=rbind(d2,d[nrow(d),]%>%add_column(S=S,cinter1=c1))
      
    }
  }
  d2[d2<10^-3]=0
  colnames(d2)=c("rho_1","rho_2",'rho_d',"rho_0","S","cinter1")
  d2$rho_plus=d2$rho_1+d2$rho_2
  
  
  Post_processing_MF(d2,for_sim$scena[i],C_seq = c_seq)
}






## 3) Hysteresis size as a function of competition coefficient ----

Run_dynamics_hysteresis=function(plot=F,N_seq_c=100,N_seq_S=100,write_table=T){
  
  color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")

  tspan = c(0, 5000)
  t = seq(0, 5000, by=1)
  julia_library("DifferentialEquations")
  julia_assign("tspan", tspan)
  
  param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=.1,emax=1.2,cintra=.1,
            cinter1=.1,cinter2=.1,S=0)
  N_seq_s=N_seq_S
  d2=tibble()
  c_seq=seq(param["cinter1"],4*param["cinter1"],length.out=N_seq_c)
  
  
  
  for (c1 in c_seq){ #for each competition coefficient 
    param[9]=c1
    
    for (branch in c("Degradation","Restoration")){ #do the bifurcation plot
    
      if (branch == "Degradation"  ){ #forward
        state = c(rho_1=0.4,rho_2=0.4,rho_d=0.1,rho_0=0.1)
        julia_assign("state", state)
        S_seq= seq(0,1,length.out=N_seq_s)
      } else { #i.e. backward
        state = c(rho_1=0.005,rho_2=0.005,rho_d=0.49,rho_0=0.5)
        julia_assign("state", state)
        S_seq= rev(seq(0,1,length.out=N_seq_s))
      }
      
      for (s in S_seq){ #each stress value
        
        param[11]=s
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
        sol = de$solve(prob,de$Tsit5(),saveat=t)
        d=as.data.frame(t(sapply(sol$u,identity)))
        colnames(d)=c("rho_1","rho_2","rho_d","rho_0")
        d2=rbind(d2,d[nrow(d),]%>%add_column(S=s,cinter1=c1,orientation=branch))
        
      }
    }
  }
  
  d2[d2<10^-3]=0
  d2$rho_plus=d2$rho_1+d2$rho_2
  d2$cinter1=d2$cinter1/param["cinter2"]
  
  if (plot){
    
    #global vegetation
    p=ggplot(d2)+geom_line(aes(x=S,y=rho_plus,linetype=orientation),lwd=.6,alpha=.8)+
      the_theme+facet_wrap(.~cinter1,labeller = label_both)+
      scale_linetype_manual(values=c(1,9))+
      labs(x="Stress (S)",y=TeX("$\\rho_{+}$"),color="")
    
    #by species
    p2=ggplot(d2%>%melt(., measure.vars=c("rho_1","rho_2"))%>%
                mutate(., variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive")))+
      geom_line(aes(x=S,y=value,color=variable,linetype=orientation),lwd=1)+
      the_theme+facet_wrap(.~cinter1,labeller = label_both)+
      scale_linetype_manual(values=c(1,9))+
      scale_color_manual(values=color_rho[c(2,4)])+
      labs(x="Stress (S)",y=TeX("$\\rho_{+}$"),color="")
    
    p_tot=ggarrange(p,p2,nrow = 2)
    ggsave("../Figures/2_species/Example_hysteresis.pdf",width = 7,height = 6)
    
  }
  
  if (write_table){
    write.table(d2,paste0("../Table/2_species/Hysteresis_size_length_C_",N_seq_c,"_length_S_",N_seq_S,".csv"),sep=";")
  }
  return(d2)
}

Run_dynamics_hysteresis(plot=T,N_seq_c=3,N_seq_S=300)

#Analysis of hysteresis size

Compute_hysteresis=function(d,n_species=2){

  nb_S=2*length(unique(d$S)) #we want both the degradation & restoration orientation
  d_hysteresis=tibble()
  
  for (i in 1:(nrow(d)/(nb_S))){ #for each bifurcation plot
    
    interval=((i-1)*nb_S+1):(i*nb_S)
    
    for (sp in 1:n_species){
      
      d2=d[interval,c(paste0("rho_",sp),"S","orientation")]
      #hysteresis size
      biomass_D=filter(d2,orientation=="Degradation")
      biomass_R=filter(d2,orientation=="Restoration")
      
      if (sp==1){ #we delete the left part of the bifurcation diagram before species become dominant
        max_S=biomass_D$S[which(biomass_D[,1]==max(biomass_D[,1]))]
        biomass_D=filter(biomass_D,S>max_S)
        biomass_R=filter(biomass_R,S>max_S)
      }
      
      hysteresis=sum(biomass_D[,1]-rev(biomass_R[,1]))
      d_hysteresis=rbind(d_hysteresis,tibble(Hysteresis=hysteresis,Species=sp))
      
      # #now tipping points
      # Tipping_D=((i-1)*n)+which(abs(diff(biomass_D))==max(abs(diff(biomass_D))))
      # size_shift_D=abs(diff(biomass_D))[which(abs(diff(biomass_D))==max(abs(diff(biomass_D))))]
      # 
      # point_5_restoration=((i-1)*n)+which(abs(diff(biomass_R))==max(abs(diff(biomass_R))))
      # size_shift_R=abs(diff(biomass_R))[which(abs(diff(biomass_R))==max(abs(diff(biomass_R))))]
      
    }
  }
  return(d_hysteresis)
}


param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=.1,emax=1.2,cintra=.1,
          cinter1=.1,cinter2=.1,S=0)
Nc=30
color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")

d=Run_dynamics_hysteresis(plot=F,N_seq_c=Nc,N_seq_S=500,write_table = T)
d_hysteresis=Compute_hysteresis(d)%>%add_column(Rel_comp=rep(seq(param["cinter2"],4*param["cinter1"],length.out=Nc),each=2))

p=ggplot(d_hysteresis%>%mutate(Species=as.character(Species))%>%
         mutate(.,Species=recode_factor(Species,"1"="stress_tol","2"="competitive")))+
  geom_line(aes(x=Rel_comp/param["cinter2"],y=Hysteresis,color=Species,group=Species),lwd=1)+
  the_theme+scale_color_manual(values=color_rho[c(2,4)])+
  labs(x="Relative competitive strength",y="Hysteresis size")

ggsave("../Figures/2_species/Evolution_hysteresis_competition.pdf",width = 7,height = 4)










## 4) Interplay between facilitation & competition  ----

# interplay between facilitation and competition on the coexistence and type of transition
# For that, we use 3 values of stress for which we do heat maps of the density of species.

state = c(rho_1=0.4,rho_2=0.4,rho_d=0.1,rho_0=0.1)
tspan = c(0, 5000)
t = seq(0, 5000, by=1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)

param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=.1,emax=1.2,cintra=.1,cinter1=.1,cinter2=.1,S=0)

c_seq=seq(param["cinter1"],4*param["cinter1"],length.out=30)
f_seq=seq(0,1.5,length.out=30);S_seq=c(0,.25,.5,.75)
d2=tibble()

for (S in S_seq){
  for (f in f_seq){
    for (c1 in c_seq){
      
      param[3]=f;param[9]=c1;param[11]=S
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob,de$Tsit5(),saveat=t)
      d=as.data.frame(t(sapply(sol$u,identity)))
      colnames(d)=c("rho_1","rho_2","rho_d","rho_0")
      d2=rbind(d2,d[nrow(d),]%>%add_column(facil=f,cinter1=c1,S=S))
      
    }
  }
}
d2[d2<10^-3]=0
d2$rho_plus=d2$rho_1+d2$rho_2


density_col=colorRampPalette(c("red","white","blue"))

p=ggplot(d2)+geom_tile(aes(x=facil,y=cinter1/param["cinter2"],fill=rho_plus))+
  theme_classic()+
  labs(x="Facilitation strength (f)", y="Relative competition strength",fill=TeX("$\\rho_{+}$"))+
  scale_fill_gradientn(colours=density_col(100))+
  facet_wrap(.~S)

ggsave("../Figures/2_species/Interplay_facilitation_competition.pdf",width = 7,height = 6)










## 5) Pair approximation between species ----




#As I get oscillations, function to get the mean densities

get_mean_densities=function(d){

  length_transient=2000
  d=d[(length_transient+1):nrow(d),]
  return(colMeans(d))
  
}


param = c(r=0.04,d=0.1,f=1,beta=0.8, m=0.1, e=0.05,emax=1.2,cintra=.02,
          cinter1=.04,cinter2=.02,S=0,delta=.1,z=4) #PA parameters from 2007 model
state = c(rho_1=0.4,rho_2=0.4,rho_m=0.1)
state_pair = c(rho_12=state[1]*state[2],rho_1m=state[1]*state[3],rho_2m=state[2]*state[3],
               rho_11=state[1]*state[1],rho_22=state[2]*state[2],rho_mm=state[3]*state[3])
state = c(state,state_pair)
names(state)= c("rho_1","rho_2","rho_m","rho_12","rho_1m","rho_2m","rho_11","rho_22","rho_mm")

tspan = c(0, 2000)
t = seq(0, 2000, by=1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)



PA_two_species_julia = julia_eval("

function PA_two_species(du,u,p,t)

    r,d,f,beta,m,e,emax,cintra,cinter1,cinter2,S,delta,z=p
    rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u
    
    #rho_1
    du[1] =  (1- rho_1-rho_2-rho_m) * beta * (delta * rho_1 + (1 - delta) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * (1 - e) ) - (cintra * rho_1 + cinter1 * rho_2)) - rho_1 * m
    
    #rho_2
    du[2] =  (1- rho_1-rho_2-rho_m) * beta * (delta * rho_2 + (1 - delta) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S) - (cintra * rho_2 + cinter2 * rho_1)) - rho_2 * m
    
    #rho_m
    du[3] =  (1- rho_1-rho_2-rho_m) * d - rho_m * (r + f * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m)))
    
    #rho_12
    du[4] = (rho_1 - rho_11 - rho_12 - rho_1m) *  (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S ) - (cintra * rho_2 + cinter2 * rho_1))  +
            (rho_2 - rho_22 - rho_12 - rho_2m) *  (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * (1 - e) ) - (cintra * rho_1 + cinter1 * rho_2))-
            2 * rho_12 * m
    
    #rho_1m
    du[5] = (rho_1 - rho_11 - rho_12 - rho_1m) * d + (rho_m - rho_mm - rho_1m - rho_2m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * (1 - e) ) - (cintra * rho_1 + cinter1 * rho_2)) -
            rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m)) ))
    
    #rho_2m
    du[6] = (rho_2 - rho_22 - rho_12 - rho_2m) * d + (rho_m - rho_mm - rho_1m - rho_2m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m)) * (emax * (1 - S * (1 - e) ) - (cintra * rho_2 + cinter2 * rho_1))) -
            rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m) )))
    
    #rho_11
    du[7] = (rho_1 - rho_11 - rho_12 - rho_1m) * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S * (1 - e) ) - (cintra * rho_1 + cinter1 * rho_2)) -
            2 * rho_11 * m
    
    #rho_22
    du[8] = (rho_2 - rho_22 - rho_12 - rho_2m) * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S ) - (cintra * rho_2 + cinter2 * rho_1)) -
            2 * rho_22 * m
    
    #rho_mm
    du[9] = 2 * (rho_m - rho_mm - rho_1m - rho_2m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m)))
  
  end
  

")

N_rep=25
d2=tibble();S_seq= seq(0,1,length.out=N_rep);c_seq=seq(param["cinter2"],20*param["cinter2"],length.out=N_rep)
for (S in S_seq){
  for (c1 in c_seq){
        
    param[11]=S;param[9]=c1
    julia_assign("p", param)
    prob = julia_eval("ODEProblem(PA_two_species, state, tspan, p)")
    sol = de$solve(prob,de$Tsit5(),saveat=t)
    d=as.data.frame(t(sapply(sol$u,identity)))
    colnames(d)= c("rho_1","rho_2","rho_m","rho_12","rho_1m","rho_2m","rho_11","rho_22","rho_mm")
    
    d2=rbind(d2,as_tibble(t(get_mean_densities(d)))%>%add_column(S=S,c1=c1))
    
  }
}
d2[d2<10^-3]=0
colnames(d2)=c("rho_1","rho_2","rho_m","rho_12","rho_1m","rho_2m","rho_11","rho_22","rho_mm","S","c1")
d2$rho_plus=d2$rho_1+d2$rho_2
#write.table(d2,"../Table/2_species/Both_competition_PA.csv",sep=";")


#postprocessing
post_processing_2species_PA=function(d2,C_seq=c_seq,S_seq=s_seq){
  

  d2$state=sapply(1:nrow(d2),function(x){
    
    if (d2$rho_1[x]>0 & d2$rho_2[x]>0) return("coexistence")
    if (d2$rho_1[x]>0 & d2$rho_2[x]==0) return("stress_tol")
    if (d2$rho_1[x]==0 & d2$rho_2[x]>0) return("competitive")
    if (d2$rho_1[x]==0 & d2$rho_2[x]==0) return('desert')
    
  })
  
  c_values_bifu=c_seq[c(1,round(length(c_seq)/2),length(c_seq))]
  
  
  color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  #state at equilibrium
  p1=ggplot(d2)+geom_tile(aes(x=S,y=c1/param["cinter2"],fill=state))+
    theme_classic()+scale_fill_manual(values=color_rho)+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="")+
    geom_hline(yintercept = c_values_bifu/param["cinter2"],lwd=.1,color="gray40")+
    theme(legend.text = element_text(size=11))
  
  
  #density of global vegetation
  density_col=colorRampPalette(c("red","white","blue"))
  p2= ggplot(d2)+geom_tile(aes(x=S,y=c1/param["cinter2"],fill=rho_plus))+
    theme_classic()+scale_fill_gradientn(colours=density_col(100))+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="Global vegetation density")
  
  #pair 12
  density_pair=colorRampPalette(c("yellow","#CE7604"))
  d2$rho_12[which(d2$rho_1==0 | d2$rho_2==0)]=NA
  p4= ggplot(d2)+geom_tile(aes(x=S,y=c1/param["cinter2"],fill=rho_12))+
    theme_classic()+scale_fill_gradientn(colours=density_pair(100),na.value = "grey")+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="Pair sp1 & sp2")
  
  #some bifurcation diagrams
  
  c_values_bifu=c_seq[c(1,round(length(c_seq)/2),length(c_seq))]
  d_bifu=filter(d2,c1 %in% c_values_bifu)
  for (c_bifu in 1:3){
    assign(paste0("p3_",c_bifu),
           ggplot(d_bifu %>% filter(.,c1==c_values_bifu[c_bifu]) %>% melt(.,measure.vars=c("rho_1","rho_2"))%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive")))+
             geom_point(aes(x=S,y=value,color=variable))+labs(x="Stress (S)",y="Density")+
             the_theme+ scale_color_manual(values=color_rho[c(2,4)])+  theme(legend.text = element_text(size=11))
    )
  }
  p3=ggarrange(p3_3,p3_2,p3_1,nrow=3,common.legend = T,legend = "bottom")
  
  
  p_tile=ggarrange(p1,p4,ncol=2)
  p_tot=ggarrange(p_tile,p3,ncol=2,widths = c(3,1))
  ggsave("../Figures/2_species/PA.pdf",width = 14,height = 6)
  
}




## 6) CA between both species ----



CA_2_species=function(landscape, param){
  
    # Variables : 1 = stress_tol, 2 = competitive, 0 = fertile, -1 = degraded   
    rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
    rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
    rho_f = length(which((landscape == 0)))   / length(landscape) #fraction fertile
    rho_d = 1-rho_1-rho_2-rho_f # fraction degraded
    
    # Neighbors :

    #using simcol package  
    neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)

    
    r=param[1];d=param[2];f=param[3];beta=param[4];m=param[5];e=param[6]
    emax=param[7];cintra=param[8];cinter2=param[9];cinter1=param[10];S=param[11];delta=param[12]
    dt=param[13];z=param[14]

    colonization_1 = beta*(delta * rho_1 + (1 - delta) * neigh_1 / z) * (emax* ( 1 - S * (1-e)) - (cintra*rho_1+cinter1*rho_2) ) * dt
    colonization_2 = beta*(delta * rho_2 + (1 - delta) * neigh_2 / z) * (emax* ( 1 - S ) - (cintra*rho_2+cinter2*rho_1) ) * dt
      
    
    # calculate regeneration, degradation & mortality rate
    death = m *dt
    regeneration = (r + f * neigh_1 / z)*dt
    degradation = d*dt 
      
    # Apply rules
    rnum = runif(length(landscape))# one random number between 0 and 1 for each cell
    landscape_update=landscape
      
    ## New stress_tol
    landscape_update[which((landscape ==0) & (rnum <= colonization_2))] = 1
    
    ## New competitive
    landscape_update[which((landscape ==0) & (rnum > (colonization_2) & rnum <= (colonization_2 + colonization_1)))] = 2
      
    ## New fertile
    landscape_update[which((landscape == 1) & (rnum <= death))] = 0
    landscape_update[which((landscape == 2) & (rnum <= death))] = 0
    landscape_update[which((landscape == -1) & (rnum <= regeneration))] = 0
      
    ## New degraded 
    landscape_update[which((landscape == 0) & (rnum > (colonization_2 + colonization_1)) & (rnum <= (colonization_2 + colonization_1 + degradation)))] = -1
  
    rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
    rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
    rho_f = length(which((landscape == 0)))    / length(landscape) #fraction fertile
    rho_d = 1-rho_1-rho_2-rho_f # fraction degraded
    
    # print(max(colonization_1))
    # print(max(colonization_2))
    # print(max(regeneration))
    
    
    if (any(c(colonization_1 + degradation,
              colonization_2 + degradation,
              colonization_2 + colonization_1 + degradation,
              death, regeneration) > 1 )) {
      warning("a set probability is exceeding 1 in run! decrease delta!!!")
    }
    if (any(c(colonization_2, colonization_1,
              degradation, regeneration,
              death) < 0)) {
      warning("a set probability falls below 0 in run balance parameters!!!")
    }
    
      
    return(list(State=c(rho_1, rho_2, rho_f, rho_d),Landscape=landscape_update))
    
}



Run_CA_2_species=function(time=seq(1,1000,1),params,ini){
  
  d=tibble(rho_1=sum(ini == 1) / length(ini),rho_2=sum(ini == 2) / length(ini),
           rho_0=sum(ini == 0) / length(ini),rho_d=sum(ini == -1) / length(ini),time=1)
  
  for (k in 2:length(time)){
    
    params["dt"]=time[k]-time[k-1]
    output_CA=CA_2_species(ini,param = params)
    ini=output_CA$Landscape
    d=rbind(d,tibble(rho_1=output_CA$State[1],rho_2=output_CA$State[2],
                     rho_0=output_CA$State[3],rho_d=output_CA$State[4],
                     time=k))

  }
  
  return(list(d=d,State=ini))
  
}


Get_initial_lattice=function(frac=c(.4,.4,.1,.1),size=25){
  return(matrix(sample(c(2,1,0,-1),replace = T,size=size*size,prob = frac),ncol=size,nrow=size))
}


Plot_landscape=function(landscape){
  color_CA=c("1"="#7BD8D3","2"="#ACD87B","0"= "#D8CC7B","-1"="#696969")
  
  ggplot(melt(landscape))+geom_tile(aes(x=Var1,y=Var2,fill=as.character(value)))+
    theme_transparent()+scale_fill_manual(values=color_CA,labels=c("Stress tol","Competitive","Fertile","Desert"))+
    theme(panel.border = element_blank())+
    theme(legend.position = "bottom")+labs(fill="")
  }


t_seq=seq(1,2000,1)

param = c(r=0.02,d=0.1,f=.9,beta=0.8, m=0.1, e=.1,emax=1,cintra=0.1,cinter1=.4,cinter2=.1,S=0,delta=.1,dt=1,z=4)

Lattice_ini=Get_initial_lattice()
test=Run_CA_2_species(time=t_seq,params = param,ini = Lattice_ini)
Plot_landscape(test$State)

d=tibble()
for (s in c(0,.25,.5,.75,1)){
  param["S"]=s
  Lattice_ini=Get_initial_lattice()
  CA=Run_CA_2_species(time = t_seq,params = param,ini = Lattice_ini)
  d=rbind(d,as_tibble(melt(CA$State))%>%add_column(.,S=s))
}

color_CA=c("1"="#7BD8D3","2"="#ACD87B","0"= "#D8CC7B","-1"="#696969")
p=ggplot(d)+geom_tile(aes(x=Var1,y=Var2,fill=as.character(value)))+
  theme_transparent()+scale_fill_manual(values=color_CA,labels=c("Stress tol","Competitive","Fertile","Desert"))+
  facet_grid(.~S)+theme(panel.border = element_blank())+
  theme(legend.position = "bottom")+labs(fill="")

ggsave("../Figures/2_species/CA_example.pdf",width = 14,height = 4)






# Step 2) Testing EWS on dynamics ----

rm(list=ls())
x=c("tidyverse","ggpubr","latex2exp","deSolve","reshape2","JuliaCall","diffeqr","tseries")
lapply(x, require, character.only = TRUE)

julia_setup()
de = diffeq_setup()


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))

plot_dynamics=function(d){
  
  color_rho=c("fertile"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  if ('time' %in% colnames(d) | 'Time' %in% colnames(d)){
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)])%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive","rho_F"="fertile","rho_D"="desert")))+
             geom_line(aes(x=time,y=value,color=variable),lwd=1)+
             theme_classic())
    

  }
  else {
    d$time=1:nrow(d)
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)])%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive","rho_F"="fertile","rho_D"="desert")))+
             geom_line(aes(x=time,y=value,color=variable),lwd=1)+
             theme_classic())
  }
}




MF_two_species_julia = julia_eval("

function MF_two_species(du,u,p,t)

  r,d,f,beta,m,e,emax,cintra,cinter1,cinter2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) + rho_2 * m + rho_1 * m

end")

Noise_MF=julia_eval("

function Noise_MF_2_species(du,u,p,t)

  # sigma = 0.01 
  
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = 0.01*u[1]
  du[2] = 0.01*u[2]
  du[3] = 0.01*u[3]
  du[4] = 0.01*u[4]
  
  
end")

param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=.1,emax=1.2,cintra=.1,cinter1=.25,cinter2=.1,S=0)
state = c(rho_1=0.4,rho_2=0.4,rho_d=0.1,rho_0=0.1)
tspan = c(0, 6000)
t = seq(0, 6000, by=1)

julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)
julia_assign("p", param)

C_seq=seq(0,1,length.out=30)
d2=tibble()

pdf("../Figures/2_species/EWS_dynamics.pdf",width = 6,height = 3)
for (S in C_seq){
  d3=tibble()
  for (rep in 1:20){
    
    param[11]=S
    julia_assign("p", param)
    
    # 
    # 
    prob = julia_eval("SDEProblem(MF_two_species,Noise_MF_2_species, state, tspan, p)")
    sol = de$solve(prob,de$Tsit5(),saveat=t)
    d=as.data.frame(t(sapply(sol$u,identity)))
    colnames(d)=c("rho_1","rho_2","rho_d","rho_0")
    d$rho_plus=d$rho_1+d$rho_2
  
    d=d[((nrow(d)-100):nrow(d)),]
  
    # print(plot_dynamics(d))
    
    #getting the variance of species abundances
    d3=rbind(d3,tibble(S=S,acf_1=acf(d$rho_1,lag=1,pl=F)$acf[2],
                       acf_2=acf(d$rho_2,lag=1,pl=F)$acf[2],
                       acf_tot=acf(d$rho_plus,lag=1,pl=F)$acf[2],
                       var_1=mean(d$rho_1)/var(d$rho_1),
                       var_2=mean(d$rho_2)/var(d$rho_2),
                       var_plus=mean(d$rho_plus)/var(d$rho_plus)))
  }
  d2=rbind(d2,colMeans(d3))
 
#Just to check if there is a shift
# 
#   prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
#   sol = de$solve(prob,de$Tsit5(),saveat=t)
#   d=as.data.frame(t(sapply(sol$u,identity)))
#   d$rho_plus=d[,1]+d[,2] # sum species
#   colnames(d)=c("rho_1","rho_2","rho_d","rho_0","rho_plus")
#   d2=rbind(d2,d[nrow(d),]%>%add_column(S=S))

  
}
dev.off()

ggplot(d2%>%add_column(., S=C_seq)%>%melt(.,id.vars="S"))+geom_line(aes(x=S,y=value,color=variable))+
  theme_classic()


      
# C) Comparing Gillespie simulations with classical IBM using 2007 model ----
rm(list=ls())
library(tidyverse)
library(simecol)
library(ggpubr)
library(latex2exp)
library(deSolve)
library(reshape2)

## 1) IBM R classic ----


Get_params=function(delta,b,c,m,d,r,f,z){
  
  return(c(
    
    delta  =  delta ,
    b      =  b     ,
    c      =  c     ,
    m      =  m     ,
    d      =  d     ,
    r      =  r     ,
    f      =  f     ,
    z      =  z
  ))
}

Get_classical_param=function(delta= 0.1,b=0.6,c= 0.3,m=.1,d= 0.2,r=0.0001,f=0.9,z=4){
  return(Get_params(     delta= delta,b=b,c= c,m=m,d= d,r=r,f=f,z=z))
}

Get_initial_lattice=function(frac=c(.8,.1,.1),size=25){
  return(matrix(sample(1:3,replace = T,size=size*size,prob = frac),ncol=size,nrow=size))
}


fourneighbors = function(landscape, state = 1, bounds = 1) {
  
  neighborhood = matrix(
    c(0, 1, 0,
      1, 0, 1,
      0, 1, 0),
    nrow = 3)
  nb_neighbors = simecol::neighbors(
    x = landscape,
    state = state,
    wdist = neighborhood,
    bounds = bounds)
  
}


#for (k in 1:length(params))assign(names(params)[k],params[[k]])

CA_Kefi = function(init, params) {
  
  # Variables : 1 = vegetation, 2 = fertile, 3 = degraded   
  landscape = init
  rho_v = sum(landscape == 1) / length(landscape)
  rho_f = sum(landscape == 2) / length(landscape)
  rho_d = 1-rho_v-rho_f
  
  # Neighbors :
  neigh_v = fourneighbors(landscape, state = 1, bounds = 1)
  neigh_f = fourneighbors(landscape, state = 2, bounds = 1)
  neigh_d = fourneighbors(landscape, state = 3, bounds = 1)
  
  
  delta =params[1];    b  =params[2];   c  =params[3];   m =params[4];    d =params[5]
  r =params[6];    f  =params[7];   z =params[8];
  colonization = (delta * rho_v + (1 - delta) * neigh_v / z) *(b - c * rho_v )*dt
  
  # calculate regeneration, degradation & mortality rate
  death = m *dt
  regeneration = (r + f * neigh_v / z)*dt
  degradation = d*dt 
  
  # Apply rules
  rnum = runif(length(landscape)) # one random number between 0 and 1 for each cell
  landscape_update = landscape
  
  ## New vegetation
  landscape_update[which(landscape == 2 & rnum <= colonization)] = 1
  
  ## New fertile
  landscape_update[which(landscape == 1 & rnum <= death)] = 2
  landscape_update[which(landscape == 3 & rnum <= regeneration)] = 2
  
  ## New degraded 
  landscape_update[which(landscape == 2 & rnum > colonization & rnum <= colonization + degradation)] = 3
  
  
  #length(which(landscape == 1 & rnum <= death))
  #length(which(landscape == 2 & rnum > colonization & rnum <= colonization + degradation))
  #length(which(landscape == 3 & rnum <= regeneration))
  #length(which(landscape == 2 & rnum <= colonization))
  
  
  return(   list(Rho_v = sum(landscape_update == 1) / length(landscape_update),
                 Rho_f = sum(landscape_update == 2) / length(landscape_update),
                 Rho_D = 1-sum(landscape_update == 1) / length(landscape_update)-sum(landscape_update == 2) / length(landscape_update),
                 Landscape=landscape_update))
}


## 2) Gillespie R ----

CA_Kefi_Gillespie = function(init, params) {
  
  
  # Variables : 1 = vegetation, 2 = fertile, 3 = degraded   
  landscape = init
  rho_v = sum(landscape == 1) / length(landscape)
  rho_f = sum(landscape == 2) / length(landscape)
  rho_d = 1-rho_v-rho_f
  
  # Neighbors :
  neigh_v = fourneighbors(landscape, state = 1, bounds = 1)
  neigh_f = fourneighbors(landscape, state = 2, bounds = 1)
  neigh_d = fourneighbors(landscape, state = 3, bounds = 1)
  
  delta =params[1];    b  =params[2];   c  =params[3];   m =params[4];    d =params[5]
  r =params[6];    f  =params[7];   z =params[8];dt=params[9]
  
  colonization = (delta * rho_v + (1 - delta) * neigh_v / z) *(b - c * rho_v )*dt
  death = m *dt
  regeneration = (r + f * neigh_v / z)*dt
  degradation = d*dt 
  
  
  
  for (i in 1:nrow(landscape)){
    for (k in 1:ncol(landscape)){
      
      #compute propensity
      if (landscape[i,j]==1){ #vegetation
        rate_tot = sum(death)
      } else if (landscape[i,j]==2){ #fertile
        rate_tot = sum(colonization[i,j],degradation) 
      }else {
        rate_tot = sum(regeneration) 
      }
      
      #random number
      rnumber=runif(1)
      time_next=-log(rnumber)/rate_tot
      
      rnumber2=runif(1)
      
      P=rnumber2*rate_tot
      
      
      
    }
  }
  
  
  
  # Apply rules
  rnum = runif(length(landscape)) # one random number between 0 and 1 for each cell
  landscape_update = landscape
  
  ## New vegetation
  landscape_update[which(landscape == 2 & rnum <= colonization)] = 1
  
  ## New fertile
  landscape_update[which(landscape == 1 & rnum <= death)] = 2
  landscape_update[which(landscape == 3 & rnum <= regeneration)] = 2
  
  ## New degraded 
  landscape_update[which(landscape == 2 & rnum > colonization & rnum <= colonization + degradation)] = 3
  
  
  #length(which(landscape == 1 & rnum <= death))
  #length(which(landscape == 2 & rnum > colonization & rnum <= colonization + degradation))
  #length(which(landscape == 3 & rnum <= regeneration))
  #length(which(landscape == 2 & rnum <= colonization))
  
  
  
  

}


params=Get_classical_param()
init=Get_initial_lattice()
