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
rm(list=ls())
x=c("tidyverse","ggpubr","latex2exp","deSolve","reshape2","JuliaCall","diffeqr")
lapply(x, require, character.only = TRUE)
source("./2_Species_analysis_functions.R")

julia_setup()
de = diffeq_setup()




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


Run_dynamics_hysteresis(plot=T,N_seq_c=3,N_seq_S=300)

#Analysis of hysteresis size
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




N_rep=50
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
write.table(d2,"../Table/2_species/PA.csv",sep=";")


#postprocessing
post_processing_2species_PA(d2,c_seq,S_seq)

## 6) CA between both species : example and clustering ----

### a) An example along interspecific competition gradient ----
t_seq=seq(1,1000,1)
param = c(r=0.02,d=0.1,f=.9,beta=0.8, m=0.05, e=.1,emax=1,cintra=0.1,cinter1=.3,cinter2=.1,S=0,delta=.1,dt=1,z=4)

Lattice_ini=Get_initial_lattice()
test=Run_CA_2_species(time=t_seq,params = param,ini = Lattice_ini)
p1=plot_dynamics(test$d)

d=tibble()
for (s in c(0,.25,.5,.75,1)){
  param["S"]=s
  Lattice_ini=Get_initial_lattice()
  CA=Run_CA_2_species(time = t_seq,params = param,ini = Lattice_ini)
  d=rbind(d,as_tibble(melt(CA$State))%>%add_column(.,S=s))
}

color_CA=c("1"="#7BD8D3","2"="#ACD87B","0"= "#D8CC7B","-1"="#696969")
p2=ggplot(d)+geom_tile(aes(x=Var1,y=Var2,fill=as.character(value)))+
  theme_transparent()+scale_fill_manual(values=color_CA,labels=c("Stress tol","Competitive","Fertile","Desert"))+
  facet_grid(.~S)+theme(panel.border = element_blank())+
  theme(legend.position = "bottom")+labs(fill="")


p_tot=ggarrange(p1,p2,ncol=2,widths = c(1,4))

ggsave("../Figures/2_species/CA_example.pdf",p_tot,width = 16,height = 4,device = "pdf")

system.time(
  for (k in 1:10){Run_CA_2_species(time = t_seq,params = param,ini = Lattice_ini)}
)

### b) Clustering metrics along inter-specific competition gradient ----

#we use the Julia Code 2_species_clustering.jl
param = c(r=0.02,d=0.1,f=.9,beta=0.8, m=0.05, e=.1,emax=1,cintra=0.1,cinter1=.3,cinter2=.1,S=0,delta=.1,dt=1,z=4)

d=read.table("../Table/2_species/Clustering_species_interspe_compet_gradient.csv",sep=",")
colnames(d)=c("rep","q12","c12","cpp","Relative_compet","S")

#q12
p1= d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    q12_m=median(q12),q12_q1=quantile(q12,0.25),q12_q3=quantile(q12,0.75),.groups = "keep")%>%
  ggplot(.)+
  geom_line(aes(x=Relative_compet/param["cinter2"],y=q12_m),color="#F1B943",lwd=1)+
  geom_ribbon(aes(x=Relative_compet/param["cinter2"],ymin=q12_q1,ymax=q12_q3),fill="#F1B943",alpha=.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~S,labeller = label_both)+labs(x="",y=TeX("$\\q_{1|2}$"))

#c12
p2= d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    c12_m=median(c12),c12_q1=quantile(c12,0.25),c12_q3=quantile(c12,0.75),.groups = "keep")%>%
  ggplot(.)+
  geom_line(aes(x=Relative_compet/param["cinter2"],y=c12_m),color="#119680")+
  geom_ribbon(aes(x=Relative_compet/param["cinter2"],ymin=c12_q1,ymax=c12_q3),fill="#119680",alpha=.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~S)+labs(x="",y=TeX("$\\c_{12}$"))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


#cpp
p3= d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    cpp_m=median(cpp),cpp_q1=quantile(cpp,0.25),cpp_q3=quantile(cpp,0.75),.groups = "keep")%>%
  ggplot(.)+
  geom_line(aes(x=Relative_compet/param["cinter2"],y=cpp_m),color="#E63535")+
  geom_ribbon(aes(x=Relative_compet/param["cinter2"],ymin=cpp_q1,ymax=cpp_q3),fill="#E63535",alpha=.5)+
  theme_classic()+theme(legend.position = "bottom")+
  facet_wrap(.~S)+labs(x="Relative competition ability",y=TeX("$\\c_{++}$"))+
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )



p_tot=ggarrange(p1,p2,p3,nrow=3)
ggsave("../Figures/2_species/Clustering_figure_interspe.pdf",width = 8,height = 6)




### c) Same but along stress gradient for c21=3*c12 ----

d=read.table("../Table/2_species/Clustering_species_stress_gradient.csv",sep=",")
colnames(d)=c("rep","q12","c12","cpp","Relative_compet","S")

#q12
p1= d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    q12_m=median(q12),q12_q1=quantile(q12,0.25),q12_q3=quantile(q12,0.75),.groups = "keep")%>%
  ggplot(.)+
  geom_line(aes(x=S,y=q12_m),color="#F1B943",lwd=1)+
  geom_ribbon(aes(x=S,ymin=q12_q1,ymax=q12_q3),fill="#F1B943",alpha=.5)+
  theme_classic()+theme(legend.position = "bottom")+
  labs(x="",y=TeX("$\\q_{1|2}$"))

#c12
p2= d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    c12_m=median(c12),c12_q1=quantile(c12,0.25),c12_q3=quantile(c12,0.75),.groups = "keep")%>%
  ggplot(.)+
  geom_line(aes(x=S,y=c12_m),color="#119680")+
  geom_ribbon(aes(x=S,ymin=c12_q1,ymax=c12_q3),fill="#119680",alpha=.5)+
  theme_classic()+theme(legend.position = "bottom")+
  labs(x="",y=TeX("$\\c_{12}$"))

#cpp
p3= d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    cpp_m=median(cpp),cpp_q1=quantile(cpp,0.25),cpp_q3=quantile(cpp,0.75),.groups = "keep")%>%
  ggplot(.)+
  geom_line(aes(x=S,y=cpp_m),color="#E63535")+
  geom_ribbon(aes(x=S,ymin=cpp_q1,ymax=cpp_q3),fill="#E63535",alpha=.5)+
  theme_classic()+theme(legend.position = "bottom")+
  labs(x="Relative competition ability",y=TeX("$\\c_{++}$"))


p_tot=ggarrange(p1,p2,p3,nrow=3)
ggsave("../Figures/2_species/Clustering_figure_stress_gradient.pdf",width = 5,height = 5)





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

No_Noise_MF=julia_eval("

function No_Noise_MF(du,u,p,t)

  # sigma = 0.01 
  
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = 0*u[1]
  du[2] = 0*u[2]
  du[3] = 0*u[3]
  du[4] = 0*u[4]
  
  
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

#pdf("../Figures/2_species/EWS_dynamics.pdf",width = 6,height = 3)
for (S in C_seq){
  d3=tibble()
  for (rep in 1:20){
    
    param[11]=S
    julia_assign("p", param)
     
    prob = julia_eval("SDEProblem(MF_two_species,Noise_MF_2_species, state, tspan, p)")
    prob = julia_eval("SDEProblem(MF_two_species,No_Noise_MF, state, tspan, p)")
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


      
# Step 3) Comparing Gillespie simulations with classical IBM  ----
rm(list=ls())
library(tidyverse)
library(simecol)
library(ggpubr)
library(latex2exp)
library(deSolve)
library(reshape2)

## 1) Comparing Gillespie and classical CA on R ----

t_seq=seq(1,1000,1)
param = c(r=0.02,d=0.1,f=.9,beta=0.8, m=0.05, e=.1,emax=1,cintra=0.1,cinter1=.1,cinter2=1,S=0.7,delta=.1,z=4,leap=0.02)
Lattice_ini=Get_initial_lattice()

t_compare=tibble()
for (k in 1:100){
  
  t1=system.time(Run_CA_2_species(t_seq,param,ini = Lattice_ini))[1]
  t2=system.time(Gillespie_tau_leaping_R(Lattice_ini,t_seq,param))[1]
  t_compare=rbind(cbind(rbind(t1,"Classic"),rbind()))
  
  
}
plot_dynamics(output$state)



## 2) Gillespie Julia ----


## Comparing computational speed ----

t_seq=seq(1,2000,1)
param = c(r=0.02,d=0.1,f=.9,beta=0.8, m=0.05, e=.1,emax=1,cintra=0.1,cinter1=.4,cinter2=.1,S=0,delta=.1,z=4,dt=1)

Lattice_ini=Get_initial_lattice()
output=Run_CA_2_species(time=t_seq,params = param,ini = Lattice_ini)
plot_dynamics(output$d)



