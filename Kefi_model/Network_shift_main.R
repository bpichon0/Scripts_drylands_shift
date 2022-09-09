rm(list=ls())
library(simecol)
library(tidyverse);library(reshape2);library(latex2exp)
library(animation);library(magick)
library(deSolve);library(rootSolve)

setwd("C:/Users/Benoi/OneDrive/Documents/Phd/Network_shift/Scripts")
source("./Network_shift_function.R")

the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))



# Step 0 : testing the model ----

params=Get_classical_param(seed = 9,a_min = 0,a_max=4,S = 3)
Plot_network(params)
Plot_species_niche(params)

state=Get_initial_conditions(params = params)

landscape=Get_initial_lattice(params = params)
Plot_landscape(landscape,params)

mortality=Get_mortality_rates(params,landscape)

ode_MF=Compute_ode(state = state,params = params,optim_time = F,n_time_no_opt = 2000)
Plot_dynamics(ode_MF,params = params)

pdf("../Figures/Species_niche.pdf",width = 6,height = 4)
for (i in 1:20){
  params=Get_classical_param(seed = i,a_min = -4,a_max=4,S = 3)
  Plot_network(params)
  Plot_species_niche(params)
  state=Get_initial_conditions(params)
  
  
  d2=tibble()
  for (b in seq(.2,.8,length.out=20)){
    for (f in c(.9)){
      params$b=b;params$f=f
      d=Compute_ode(state,params,optim_time = F)
      d=d[nrow(d),-1]
      d2=rbind(d2,d %>% add_column(., b=b,f=f))
    }
  }
  d2$Rho_V=rowSums(d2[,1:params$S])
  
  p=ggplot(melt(d2,measure.vars = c(paste0("Rho_V",1:params$S),'Rho_V')))+geom_line(aes(x=b,color=variable,y=value),shape=1)+theme_classic()+theme(legend.position = "bottom")+
    labs(x='aridity (b)',y=TeX("$\\rho_{+}$"),color="f")
  print(p)
}
dev.off()


# Step 1 : Some bifurcations ----
S=3 #number of species

d2=tibble();params=Get_classical_param(S=S,seed = 9,a_min = 0,a_max=4);state=Get_initial_conditions(params)

for (b in seq(.2,.8,length.out=20)){
  for (f in c(.9)){
    params$b=b;params$f=f
    d=Compute_ode(state,params)
    d=d[nrow(d),-1]
    d2=rbind(d2,d %>% add_column(., b=b,f=f))
  }
}
d2$Rho_V=rowSums(d2[,1:3])

p=ggplot(melt(d2,measure.vars = c(paste0("Rho_V",1:params$S),'Rho_V')))+geom_point(aes(x=b,color=variable,y=value),shape=1)+theme_classic()+theme(legend.position = "bottom")+
  labs(x='aridity (b)',y=TeX("$\\rho_{+}$"),color="f")
ggsave(paste0("../Figures/Com_",S,"_sp.pdf"),width = 7,height = 4)


#Testing different motifs

S=3 #number of species
for (i in 1:4){
  
  d2=tibble();params=Get_classical_param(S=S,seed = 9,a_min = 0,a_max=4);state=Get_initial_conditions(params)
  
  params$A=Get_modules(i)
  
  for (b in seq(.2,.8,length.out=20)){
    for (f in c(.9)){
      params$b=b;params$f=f
      d=Compute_ode(state,params)
      d=d[nrow(d),-1]
      d2=rbind(d2,d %>% add_column(., b=b,f=f))
    }
  }
  d2$Rho_V=rowSums(d2[,1:3])
  
  p=ggplot(melt(d2,measure.vars = c(paste0("Rho_V",1:params$S),'Rho_V')))+geom_line(aes(x=b,color=variable,y=value),shape=1)+theme_classic()+theme(legend.position = "bottom")+
    labs(x='aridity (b)',y=TeX("$\\rho_{+}$"),color="f")
  ggsave(paste0("../Figures/Com_sp_module_",i,".pdf"),width = 7,height = 4)
  
}




