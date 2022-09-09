rm(list=ls())
library(simecol)
library(tidyverse);library(reshape2);library(latex2exp)
library(animation);library(magick)
library(deSolve);library(rootSolve)

setwd("C:/Users/Benoi/OneDrive/Documents/Phd/Network_shift/Scripts")


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90, face = "italic"),
                                strip.text.x = element_text(size = 10, face = "italic"),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))

# Step 1 : Functions -----


Get_params=function(delta,b,c,m,d,r,f,z){
  
  return(list(
    
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
  

  with(params, {
    
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
  })
  
}


Run_CA_kefi=function(time=seq(1,1000,1),params,ini,plot=F){
  
  
  d=tibble(Time=1,Rho_V=sum(ini == 1) / length(ini),Rho_F=sum(ini == 2) / length(ini),Rho_D=sum(ini == 3) / length(ini))
  state=list(Landscape=ini,Rho_v=d$Rho_V,Rho_f=d$Rho_F,Rho_D=d$Rho_D)
  
  for (k in 2:length(time)){
    
    params$dt=time[k]-time[k-1]
    
    state=CA_Kefi(state$Landscape,params = params)
    if (plot==T & k%%200==0){
      png(paste0("../Figures/State_",k,".png"))
      Plot_landscape(state$Landscape)
      dev.off()
    }
    
    d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
  }
  
  return(list(d=d,State=state$Landscape))
  
}

Plot_dynamics=function(d,different_sim=F,simple=F){
  
  if (different_sim==T & simple==F){  
    d_melt=melt(d,id.vars=c("Time","Type"))
    p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
      geom_line(aes(x=Time,y=value,color=variable,linetype=Type),lwd=1)+scale_color_manual(values=c("#9DD696","#D2C48F","#777777"))+
      labs(x="Time steps",y="Fraction",color="")+the_theme
    return(p)
    
  } 
  
  if (different_sim==F & simple==F){  
    
    d_melt=melt(d,id.vars=c("Time"))
    p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
      geom_line(aes(x=Time,y=value,color=variable))+scale_color_manual(values=c("#9DD696","#D2C48F","#777777"))+
      labs(x="Time steps",y="Fraction",color="")+the_theme
    return(p)
    
  }
  
  if (simple == T){
    plot(d$Time,d$Rho_V,xlab="Time steps","l",ylab="Fraction",ylim=c(0,1),col="#9DD696")
    lines(d$Time,d$Rho_F,ylim=c(0,1),col="#D2C48F")
    lines(d$Time,d$Rho_D,ylim=c(0,1),col="#777777")
  }
    
}

Plot_landscape=function(landscape){
  image(landscape,xaxt = "n",yaxt ="n",col=c("#9DD696","#D2C48F","#777777") )
}


Mean_field_kefi=function(t,state,param){
  
  names(state)=c("rho_p","rho_m")
  
  with (as.list(c(state, param)),{
    
    Drho_p=rho_p*(b-c*rho_p)*(1-rho_p-rho_m)-m*rho_p
    Drho_m=d*(1-rho_p-rho_m)-(r+f*rho_p)*rho_m
    
    
    list(c(    Rho_p  = Drho_p    ,  Rho_p  = Drho_m))
  })
  
}

Function_kefi_model=function(rho,param){
  
  with (as.list(c(param)),{
    
    Drho_p=rho[1]*(b-c*rho[1])*(1-rho[1]-rho[2])-m*rho[1]
    Drho_m=d*(1-rho[1]-rho[2])-(r+f*rho[1])*rho[2]
    
    
    c(  Drho_p    ,Drho_m)
  })
  
}


Compute_ode=function(state,params,method_ode="lsoda",optim_time=T,h=.1){
  
  
  if ( optim_time==T){
    n_time_ode=0
    count_non_eq=1
    
    while(count_non_eq!=0 & n_time_ode<5000){ #while we are not at equilibrium for all the organisms/nutrients/detritus
      
      
      n_time_ode=n_time_ode+1000
      time=seq(0,n_time_ode,h)
      count_non_eq=0
      
      dynamics = as.data.frame(ode(state,time,func=Mean_field_kefi,parms=unlist(params) ,method = method_ode))
     
      
      #check if we are at equilibrium : fit a linear model at each "patch" and 
      #see if the slope is ~ 0. If not : increase time of integration
      
      length_lm=(1/h)*(n_time_ode-300):n_time_ode 
      
      cols_to_lm=(1+which(state!=0)) #values of cols for which we perform the linear regression
      
      for (c in cols_to_lm){#we check whether the last 300 values are equal by fitting a linear model
        data=data.frame(time=dynamics[length_lm,1],value=dynamics[length_lm,c])
        slope_patch=lm(formula = value~time,data)$coefficients[2] #slope of the linear regression
        if (abs(slope_patch)>1e-5){ #Threshold for considering equilibrium state
          count_non_eq=count_non_eq+1
        }
      }
    }
    
    final_point=as.numeric(dynamics[nrow(dynamics),-1])
    dynamics[nrow(dynamics),-1]=final_point
    colnames(dynamics)=c("Time","Rho_V","Rho_D")
    
  } 
  
  dynamics$Rho_F=1-dynamics$Rho_D-dynamics$Rho_V
  dynamics=dynamics[,c(1,2,4,3)]
  return(as_tibble(dynamics))
}

Make_gif_landscape=function(path="../Figures/"){
  
  imgs=list.files(path, full.names = TRUE,pattern = ".png")
  img_list=lapply(imgs, image_read)
  
  ## join the images together
  img_joined=image_join(img_list)
  
  ## animate at 2 frames per second
  img_animated=image_animate(img_joined, fps = 2)
  
    ## save to disk
  image_write(image = img_animated,
              path = paste0(path,"/Gif/Landscape.gif"))
  
  
}



# Step 2 : CA & mean-field -----

set.seed(3)
params=Get_classical_param(r = 0,d = .1,c = .2,delta = .1,f = .9,z = 4,m = .2,b = .57)

ini_land=  Get_initial_lattice(frac = c(.8,.1,.1),size=25)

d_CA=Run_CA_kefi(time = seq(1,1000,by=1),params = params,ini=ini_land,plot = F)

Plot_dynamics(d_CA$d,simple = T)
Plot_landscape(d_CA$State)

system.time(for (k in 1:10) Run_CA_kefi(time = seq(1,1000,by=1),params = params,ini=ini_land,plot = F))


state=c(rho_p=.8,rho_m=.1)
d_meanfield=Compute_ode(state,params)

Plot_dynamics(d_meanfield)

d_meanfield$Type="Mean field"
d_CA$d$Type="CA"
d_tot=rbind(d_meanfield,d_CA$d)

Plot_dynamics(d_tot,different_sim =T)



# Diag bifu de Kefi et al., 2007 Theo. Ecol.

d2=tibble();params=Get_classical_param(r = 0.05,d = .2,c = .3,delta = .1,f = .9,z = 4,m = .1,b = .57);state=c(.9,.05)

for (b in seq(.2,.8,length.out=30)){
  for (f in c(.1,.3,.9)){
    params$b=b;params$f=f
    d=Compute_ode(state,params)
    d=d[nrow(d),-1]
    d2=rbind(d2,d %>% add_column(., b=b,f=f))
  }
}

ggplot(d2)+geom_point(aes(x=b,color=as.factor(f),y=Rho_V),shape=1)+theme_classic()+theme(legend.position = "bottom")+
  labs(x='aridity (b)',y=TeX("$\\rho_{+}$"),color="f")


params=Get_classical_param(r = 0.05,d = .2,c = .3,delta = .1,f = .9,z = 4,m = .1,b = .57);Eq=tibble()

for (b in seq(.2,.8,length.out=500)){
  for (f in c(.1,.3,.9)){
    
    params$b=b;params$f=f;eq=c()
    
    for (init in seq(0.1,1,length.out=10)){
      eq=c(eq,multiroot(Function_kefi_model,start = c(init,1-init-.1),param=params)$root[1])}
    
    eq=unique(round(eq[eq>=0 & eq<=1],2))
    if (length(eq)>1){stability=c("Unstable","Stable")}else{stability="Stable"}
    Eq=rbind(Eq,tibble(Eq=eq,b=rep(params$b,length(eq)),f=rep(params$f,length(eq)),Stability=stability))
  }
}

ggplot(Eq)+geom_point(aes(x=b,y=Eq,color=Stability,shape=as.factor(f)))+theme_classic()+scale_color_manual(values=c("black","gray70"))


# Gif of landscapes

d2=tibble();params=Get_classical_param(r = 0.05,d = .2,c = .3,delta = .1,f = .9,z = 4,m = .1,b = .57);state=c(.9,.05);k=1

for (b in rev(seq(.2,.8,length.out=10))){
  for (f in c(.9)){
    params$b=b;params$f=f
    
    d=Run_CA_kefi(params = params,size = 50,plot = F)
    
    png(paste0("../Figures/State_",k,".png"),width = 800,height = 400)
    par(mfrow=c(1,2))
    Plot_landscape(d$State)
    Plot_dynamics(d$d,simple=T)
    mtext(paste0("b = ",params$b," & f = ",params$f))
    dev.off()
    
    k=k+1
    
  }
}

Make_gif_landscape()


# Step 3 : test with Julia ----



