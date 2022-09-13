# Other) Recruitment rate ----
rm(list=ls())
library(tidyverse)
d=expand_grid(emax=1,e=.5,S=seq(0,1,length.out=200),psi=c(0,.5,1))

p=ggplot(d)+geom_line(aes(x=S,y=emax*(1-S*(1-e*psi)),color=as.factor(psi),group=psi))+
scale_color_manual(values=c("blue","green","red"))+theme_classic()+
theme(legend.position = "bottom")+labs(x="S",
y=TeX("$\\epsilon_{max} (1-S (1-e \\psi))$"),color=TeX("\\psi"))



ggsave("../Figures/Recruitment_rate.pdf",p,width = 6,height =4 )



# A) Two species MF model ----
rm(list=ls())
library(tidyverse)
library(ggpubr)
library(latex2exp)
library(deSolve)
library(reshape2)
library(JuliaCall)
julia_setup()
library(diffeqr)
de = diffeq_setup()


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))



# MF_two_species=function(t,state,param){
#   with(as.list(c(state, param)), {
#     
#     state[state < 10^(-8)] = 0
#     drho_1 = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c*(rho_2 + rho_1))) - rho_1 * m 
#     drho_2 = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c*(rho_2 + rho_1))) - rho_2 * m 
#     drho_d = d*rho_0 - rho_d*(r+f*(rho_1))
#     drho_0 = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c*(rho_2 + rho_1)))  -
#       rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c*(rho_2 + rho_1))) + rho_2 * m + rho_1 * m 
#     
#   list(c(drho_1,drho_2,drho_d,drho_0))
#   })
# }

# 1) Simple bifu along stress ----
 

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
















# 2) State diagram with competitive ability and stress for different type of competition ----

#We now try to replicate Danet figure with abiotic stress in x-axis and differential competitive ability in y-axis

param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=0.1,emax=1.2,c1=0.3,c2=.3,S=0)
state = c(rho_1=0.4,rho_2=0.4,rho_d=0.1,rho_0=0.1)
tspan = c(0, 10000)
t = seq(0, 10000, by=1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)

#Intra + interspecific competition
MF_two_species_julia = julia_eval("

function MF_two_species(du,u,p,t)

  r,d,f,beta,m,e,emax,c1,c2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c1*(rho_2 + rho_1))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c2*(rho_2 + rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c1*(rho_2 + rho_1)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c2*(rho_2 + rho_1))) + rho_2 * m + rho_1 * m

end")

#Asymetric competition

# MF_two_species_julia = julia_eval( "
# 
# function MF_two_species(du,u,p,t)
# 
#   r,d,f,beta,m,e,emax,c1,c2,S=p
#   rho_1,rho_2,rho_d,rho_0=u
# 
#   du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c1*(rho_2+rho_1))) - rho_1 * m
#   du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c2*(rho_2))) - rho_2 * m
#   du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
#   du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c1*(rho_2+rho_1)))  -
#       rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c2*(rho_2))) + rho_2 * m + rho_1 * m
# 
# end")


#Interspecific competition

# MF_two_species_julia = julia_eval( "
# 
# function MF_two_species(du,u,p,t)
# 
#   r,d,f,beta,m,e,emax,c1,c2,S=p
#   rho_1,rho_2,rho_d,rho_0=u
# 
#   du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c1*(rho_2))) - rho_1 * m
#   du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c2*(rho_1))) - rho_2 * m
#   du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
#   du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c1*(rho_2)))  -
#       rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c2*(rho_1))) + rho_2 * m + rho_1 * m
# 
# end")

d2=tibble();S_seq= seq(0,1,length.out=100);c_seq=seq(.3,.7,length.out=100)
for (S in S_seq){
  for (c1 in c_seq){
    
    param[10]=S;param[8]=c1
    julia_assign("p", param)
    prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
    sol = de$solve(prob,de$Tsit5(),saveat=t)
    d=as.data.frame(t(sapply(sol$u,identity)))
    d2=rbind(d2,d[nrow(d),]%>%add_column(S=S,c1=c1))
    
  }
}
d2[d2<10^-3]=0
colnames(d2)=c("rho_1","rho_2",'rho_d',"rho_0","S","c1")
d2$rho_plus=d2$rho_1+d2$rho_2
write.table(d2,"../Table/2_species/Both_competition.csv",sep=";")
# write.table(d2,"../Table/2_species/Asymmetric_competition.csv",sep=";")
# write.table(d2,"../Table/2_species/Interspe_competition.csv",sep=";")

#postprocessing
d2$state=sapply(1:nrow(d2),function(x){
  
  if (d2[x,1]>0 & d2[x,2]>0) return("coexistence")
  if (d2[x,1]>0 & d2[x,2]==0) return("stress_tol")
  if (d2[x,1]==0 & d2[x,2]>0) return("competitive")
  if (d2[x,1]==0 & d2[x,2]==0) return('desert')
  
})

c_values_bifu=c_seq[c(1,round(length(c_seq)/2),length(c_seq))]


color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")

#state at equilibrium
p1=ggplot(d2)+geom_tile(aes(x=S,y=c1/param["c2"],fill=state))+
  theme_classic()+scale_fill_manual(values=color_rho)+
  theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="")+
  geom_hline(yintercept = c_values_bifu/param["c2"],lwd=.1,color="gray40")+
  theme(legend.text = element_text(size=11))


#density of global vegetation
density_col=colorRampPalette(c("red","white","blue"))
p2= ggplot(d2)+geom_tile(aes(x=S,y=c1/param["c2"],fill=rho_plus))+
  theme_classic()+scale_fill_gradientn(colours=density_col(100))+
  theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="Global vegetation density")


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


p_tile=ggarrange(p1,p2,ncol=2)
p_tot=ggarrange(p_tile,p3,ncol=2,widths = c(3,1))
ggsave("../Figures/2_species/Trade_off_competitive_ability.pdf",width = 14,height = 6)
# ggsave("../Figures/2_species/Trade_off_competitive_ability_asymmetry.pdf",width = 14,height = 6)
# ggsave("../Figures/2_species/Trade_off_competitive_ability_interspe.pdf",width = 14,height = 6)












# 3) Pair approximation between species ----


plot_dynamics=function(d){
  if ('time' %in% colnames(d) | 'Time' %in% colnames(d)){
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)]))+
                    geom_line(aes(x=time,y=value,color=variable),lwd=1)+
      theme_classic())
  }
  else {
    d$time=1:nrow(d)
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)]))+
                    geom_line(aes(x=time,y=value,color=variable),lwd=1)+
      theme_classic())
  }
}


#As I get oscillations, function to get the mean densities

get_mean_densities=function(d){

  length_transient=2000
  d=d[(length_transient+1):nrow(d),]
  return(colMeans(d))
  
}


param = c(r=0.04,d=0.1,f=1,beta=0.8, m=0.1, e=0.05,emax=1.2,c1=0.02,c2=.02,S=0,delta=.1,z=4) #PA parameters from 2007 model
state = c(rho_1=0.4,rho_2=0.4,rho_m=0.1)
state_pair = c(rho_12=state[1]*state[2],rho_1m=state[1]*state[3],rho_2m=state[2]*state[3],
               rho_11=state[1]*state[1],rho_22=state[2]*state[2],rho_mm=state[3]*state[3])
state = c(state,state_pair)
names(state)= c("rho_1","rho_2","rho_m","rho_12","rho_1m","rho_2m","rho_11","rho_22","rho_mm")

tspan = c(0, 10000)
t = seq(0, 10000, by=1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)



PA_two_species_julia = julia_eval("

function PA_two_species(du,u,p,t)

    r,d,f,beta,m,e,emax,c1,c2,S,delta,z=p
    rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u
    
    #rho_1
    du[1] =  (1- rho_1-rho_2-rho_m) * beta * (delta * rho_1 + (1 - delta) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * (1 - e) ) - c1 * (rho_1 + rho_2)) - rho_1 * m
    
    #rho_2
    du[2] =  (1- rho_1-rho_2-rho_m) * beta * (delta * rho_2 + (1 - delta) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S) - c2 * (rho_2 + rho_2)) - rho_2 * m
    
    #rho_m
    du[3] =  (1- rho_1-rho_2-rho_m) * d - rho_m * (r + f * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m)))
    
    #rho_12
    du[4] = (rho_1 - rho_11 - rho_12 - rho_1m) *  (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S ) - c2 * (rho_1 + rho_2))  +
            (rho_2 - rho_22 - rho_12 - rho_2m) *  (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * (1 - e) ) - c1 * (rho_1 + rho_2))-
            2 * rho_12 * m
    
    #rho_1m
    du[5] = (rho_1 - rho_11 - rho_12 - rho_1m) * d + (rho_m - rho_mm - rho_1m - rho_2m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * (1 - e) ) - c1 * (rho_1 + rho_2)) -
            rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m)) ))
    
    #rho_2m
    du[6] = (rho_2 - rho_22 - rho_12 - rho_2m) * d + (rho_m - rho_mm - rho_1m - rho_2m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m)) * (emax * (1 - S * (1 - e) ) - c2 * (rho_1 + rho_2))) -
            rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m) )))
    
    #rho_11
    du[7] = (rho_1 - rho_11 - rho_12 - rho_1m) * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S * (1 - e) ) - c1 * (rho_1 + rho_2)) -
            2 * rho_11 * m
    
    #rho_22
    du[8] = (rho_2 - rho_22 - rho_12 - rho_2m) * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S ) - c2 * (rho_1 + rho_2)) -
            2 * rho_22 * m
    
    #rho_mm
    du[9] = 2 * (rho_m - rho_mm - rho_1m - rho_2m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m) / (1- rho_1-rho_2-rho_m)))
  
  end
  

")

N_rep=50
d2=tibble();S_seq= seq(0,1,length.out=N_rep);c_seq=seq(.02,1,length.out=N_rep)
for (S in S_seq){
  print(S)
  for (c1 in c_seq){
        
    param[10]=S;param[8]=c1
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
write.table(d2,"../Table/2_species/Both_competition_PA.csv",sep=";")


#postprocessing
post_processing_2species=function(d2,C_seq=C_seq,S_seq=S_seq){
  

  d2$state=sapply(1:nrow(d2),function(x){
    
    if (d2$rho_1[x]>0 & d2$rho_2[x]>0) return("coexistence")
    if (d2$rho_1[x]>0 & d2$rho_2[x]==0) return("stress_tol")
    if (d2$rho_1[x]==0 & d2$rho_2[x]>0) return("competitive")
    if (d2$rho_1[x]==0 & d2$rho_2[x]==0) return('desert')
    
  })
  
  c_values_bifu=c_seq[c(1,round(length(c_seq)/2),length(c_seq))]
  
  
  color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  #state at equilibrium
  p1=ggplot(d2)+geom_tile(aes(x=S,y=c1/param["c2"],fill=state))+
    theme_classic()+scale_fill_manual(values=color_rho)+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="")+
    geom_hline(yintercept = c_values_bifu/param["c2"],lwd=.1,color="gray40")+
    theme(legend.text = element_text(size=11))
  
  
  #density of global vegetation
  density_col=colorRampPalette(c("red","white","blue"))
  p2= ggplot(d2)+geom_tile(aes(x=S,y=c1/param["c2"],fill=rho_plus))+
    theme_classic()+scale_fill_gradientn(colours=density_col(100))+
    theme(legend.position = "bottom")+labs(x="Stress",y="Relative competitive ability",fill="Global vegetation density")
  
  #pair 12
  density_pair=colorRampPalette(c("yellow","#CE7604"))
  d2$rho_12[which(d2$rho_1==0 | d2$rho_2==0)]=NA
  p4= ggplot(d2)+geom_tile(aes(x=S,y=c1/param["c2"],fill=rho_12))+
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
  ggsave("../Figures/2_species/Trade_off_competitive_ability_PA.pdf",width = 14,height = 6)
  
}




# 4) CA between both species ----


param = c(r=0.05,d=0.1,f=0.9,beta=0.8, m=0.1, e=0.1,emax=1.2,c1=0.3,c2=.3,S=0,delta=.1,dt=1,z=4)

CA_2_species=function(landscape, param){
  
    # Variables : 1 = stress_tol, 2 = competitive, 0 = fertile, -1 = degraded   
    rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
    rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
    rho_f = length(which((landscape == 0)))    / length(landscape) #fraction fertile
    rho_d = 1-rho_1-rho_2-rho_f # fraction degraded
    
    # Neighbors :

    #using simcol package  
    neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)

    
    r=param[1];d=param[2];f=param[3];beta=param[4];m=param[5];e=param[6]
    emax=param[7];c1=param[8];c2=param[9];S=param[10];delta=param[11]
    dt=param[12];z=param[13]

    colonization_1 = beta*(delta * rho_1 + (1 - delta) * neigh_1 / z) * (emax* ( 1 - S * (1-e)) - c1 * (rho_1+rho_2) ) * dt
    colonization_2 = beta*(delta * rho_2 + (1 - delta) * neigh_2 / z) * (emax* ( 1 - S        ) - c2 * (rho_1+rho_2) ) * dt
      
    
    # calculate regeneration, degradation & mortality rate
    death = m *dt
    regeneration = (r + f * neigh_1 / z)*dt
    degradation = d*dt 
      
    # Apply rules
    rnum = runif(length(landscape))# one random number between 0 and 1 for each cell
    landscape_update=landscape
      
    ## New stress_tol
    landscape_update[which((landscape ==0) & (rnum <= colonization_1))] = 1
    
    ## New competitive
    landscape_update[which((landscape ==0) & (rnum <= colonization_2 + colonization_1))] = 2
      
    ## New fertile
    landscape_update[which((landscape == 1) & (rnum <= death))] = 0
    landscape_update[which((landscape == 2) & (rnum <= death))] = 0
    landscape_update[which((landscape == -1) & (rnum <= regeneration))] = 0
      
    ## New degraded 
    landscape_update[which((landscape == 2) & (rnum > colonization_2 + colonization_1) & (rnum <= (colonization_2 + colonization_1 + degradation)))] = 3
  
    rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
    rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
    rho_f = length(which((landscape == 0)))    / length(landscape) #fraction fertile
    rho_d = 1-rho_1-rho_2-rho_f # fraction degraded
      
    return(list(State=c(rho_1, rho_2, rho_f, rho_d),Landscape=landscape_update))
    
}



Run_CA_2_species=function(time=seq(1,2000,1),params,ini){
  
  d=tibble(Time=1,Rho_1=sum(ini == 1) / length(ini),Rho_2=sum(ini == 2) / length(ini),
           Rho_F=sum(ini == 0) / length(ini),Rho_D=sum(ini == -1) / length(ini))
  
  for (k in 2:length(time)){
    
    params["dt"]=time[k]-time[k-1]
    output_CA=CA_2_species(ini,param = params)
    ini=output_CA$Landscape
    d=rbind(d,tibble(Time=k,Rho_1=output_CA$State[1],Rho_2=output_CA$State[2],Rho_F=output_CA$State[3],Rho_D=output_CA$State[4]))

  }
  
  return(list(d=d,State=ini))
  
}


Get_initial_lattice=function(frac=c(.4,.4,.1,.1),size=25){
  return(matrix(sample(c(2,1,0,-1),replace = T,size=size*size,prob = frac),ncol=size,nrow=size))
}


Plot_landscape=function(landscape){
  image(landscape,xaxt = "n",yaxt ="n",col=c("1"="#7BD8D3","2"="#ACD87B","0"="#FFECBC","-1"="#696969") )
}


Plot_dynamics=function(d){
  
  d_melt=melt(d,id.vars=c("Time"))
  p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_1"="stress_tol","Rho_2"="competitive","Rho_F"="fertile","Rho_D"="desert")))+
    geom_line(aes(x=Time,y=value,color=variable))+scale_color_manual(values=c(color_rho[2:4],"fertile"="#FFECBC"))+
    labs(x="Time steps",y="Fraction",color="")+theme_classic()+theme(legend.position = "bottom")

}

Lattice_ini=Get_initial_lattice()
test=Run_CA_2_species(params = param,ini = Lattice_ini)
