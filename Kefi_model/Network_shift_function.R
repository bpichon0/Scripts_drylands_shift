library(simecol)
library(tidyverse);library(reshape2);library(latex2exp)
library(animation);library(magick)
library(deSolve);library(rootSolve);library(igraph);library(rethinking)



## Parameters ----



Get_params=function(delta,b,c,d,r,f,z,S,teta,m0,seed,a_max,a_min){
  
  set.seed(seed)
  
  m0  = m0 
  A   = Get_interaction_coefficient(S,seed = seed,a_max=a_max,a_min=a_min)
  
  return(list(
    b_opt  =  runif(n=S,0.1,.9),
    Tol_sp =  runif(n=S,.01,.7),
    m0     =  m0    ,
    teta   =  teta  ,
    A      =  A     ,
    delta  =  delta ,
    b      =  b     ,
    c      =  c     ,
    m0     =  m0    ,
    d      =  d     ,
    r      =  r     ,
    f      =  f     ,
    z      =  z     ,
    S      =  S     
    
  ))
}

Get_classical_param=function(delta= 0.1,b=0.6,c= 0.3,d= 0.2,r=0.0001,f=0.9,z=4,m0=.15,S=3,teta=.5,seed=8,a_max=4,a_min=-4){
  return(Get_params(     delta= delta,b=b,c= c,d= d,r=r,f=f,z=z,teta=teta,S=S,m0=m0,seed=seed,a_max=a_max,a_min=a_min))
}

Get_interaction_coefficient=function(S,symetry=T,seed=8,a_min=-4,a_max=4){
  
  "To ADD : case of asymmetric interactions"
  
  set.seed(seed)
  nS=S
  
  alpha=runif(n=nS*(nS-1)/2,a_min,a_max) #number of pair interactions between species
  
  A = matrix(0,nrow=nS,ncol=nS)
  A[which(upper.tri(x = A,diag = F)==T)]   = alpha
  
  if (symetry==T){
    A[which(lower.tri(x = A,diag = F)==T)] = alpha
  }
  
  else {
    set.seed(seed+1)
    A[which(lower.tri(x = A,diag = F)==T)] = runif(n=nS,-4,4)
  }
  
  return (A)
}

Get_modules=function(type_module=1){
  
  # Cf Losapio et al., PNAS 2021
  
  # intransitive competition
  mat1=matrix(c(0,0,-1,-1,0,0,0,-1,0),3,3)
  
  # two facilitation, 1 competition
  mat2=t(matrix(c(0,0,0,1,0,-1,1,-1,0),3,3))
  
  # two competition, 1 facilitation
  mat3=matrix(c(0,0,0,-1,0,1,-1,1,0),3,3)
  
  # two competition, 1 facilitation bis
  mat4=matrix(c(0,-1,-1,0,0,1,0,1,0),3,3)
  
  if (type_module==1) return(mat1)
  if (type_module==2) return(mat2)
  if (type_module==3) return(mat3)
  if (type_module==4) return(mat4)
  
  
}


Get_mortality_rates=function(params,landscape,seed=8){
  
  set.seed(seed)
  
  nS       =  params$S     # number species
  b_opt    =  params$b_opt # optimal aridity
  theta    =  params$teta  # coefficient of decrease
  m0       =  params$m0    # basal mortality
  Tol_sp   =  params$Tol_sp
  
  mat=list() # neighbors for each species i
  
  for (i in 1:params$S){
    mat[[i]]=fourneighbors(landscape, state = i, bounds = 1)
  }
  
  mortality_interaction=mortality_aridity=matrix(0,nrow(landscape),ncol(landscape))
  for (nr1 in 1:nrow(landscape)){
    for (nr2 in 1:nrow(landscape)){
      
      if (landscape[nr1,nr2] %in% 1:params$S){
        
        #Effect of interactions
        
        coeff_interac =  params$A[landscape[nr1,nr2],] #coefficient of interaction the species in the patch (nr1,nr2)
        mortality_interaction[nr1,nr2] =  1  -    (sum(sapply(1:params$S,function(x,...){ #coefficient of interaction weighted by the fraction of species in the neighbouring 
          return( (coeff_interac[x] * ( mat[[x]][nr1,nr2]/params$z ))/((abs((coeff_interac[x] * ( mat[[x]][nr1,nr2]/params$z ))))+1)  )
        })))
        
        #Effect of aridity
        
        mortality_aridity[nr1,nr2]=1-params$teta*exp(-((params$b-b_opt[landscape[nr1,nr2]])/Tol_sp[landscape[nr1,nr2]])**2)
        
      } else {
        mortality_aridity[nr1,nr2]=mortality_interaction[nr1,nr2]=0
      }
      
      
    }
  }
  
  mortality=m0*mortality_aridity*mortality_interaction
  
  
  
  return(mortality)
  
}



## Lattice ----


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


Get_initial_lattice=function(frac=c(.1,.1,.8),size=25,params=params){
  
  #frac being the relative proportion of vegetation, fertile and degraded space respectively
  
  prob_each_sp=rep(frac[1]/params$S,length.out=params$S)
  
  newfrac=c(frac[1],frac[2],prob_each_sp)
  
  return(matrix(sample(c(-1,0,1:params$S),replace = T,size=size*size,prob = newfrac),ncol=size,nrow=size))
}



Get_initial_conditions=function(params,frac=c(.8,.1,.1)){
  
  #frac being the relative proportion of vegetation, fertile and degraded space respectively
  
  prob_each_sp=rep(frac[1]/params$S,length.out=params$S)
  
  newfrac=c(prob_each_sp,frac[2],frac[3])
  
  return(newfrac)
}



# CA ----

#for (i in 1:length(params)) assign(names(params[i]),params[i][[1]])






CA_network = function(landscape, params) {
  
  # Variables : 1 : S = vegetation, 0 = fertile, -1 = degraded   
  rho_v = sum(landscape %in% 1:params$S) / length(landscape)
  rho_f = sum(landscape == 0) / length(landscape)
  rho_d = 1-rho_v-rho_f
  
  # Neighbors & fraction of each species
  
  neigh=list();rho=c() # neighbors & global density for each species i
  for (i in 1:params$S){
    neigh[[i]]=fourneighbors(landscape, state = i, bounds = 1)
    rho[i]=sum(landscape ==i) / length(landscape)
  }
  neigh_f = fourneighbors(landscape, state =  0, bounds = 1)  # fertile cells
  neigh_d = fourneighbors(landscape, state = -1, bounds = 1)  # degraded cells
  neigh_v = sum(neigh[[1]],neigh[[2]],neigh[[3]])             # vegetation cells
  
  
  #mortality
  mortality=Get_mortality_rates(params,landscape);params$mortality=mortality
  
  with(params, { #dt being the time step between two events
    
    colonization=list()
    for (i in 1:S){
      colonization[[i]]=(delta * rho[i] + (1 - delta) * neigh[[i]] / z) *(b - c * rho[i] )*dt
    }
    total_colonization=colonization[[1]]
    for (sp in 2:S){
      total_colonization=total_colonization+colonization[[sp]]
    }
    
    # calculate regeneration, degradation & mortality rate
    death =  mortality*dt
    regeneration = (r + f * neigh_v / z)*dt
    degradation = d*dt 
    
    # Apply rules
    rnum = runif(length(landscape)) # one random number between 0 and 1 for each cell
    landscape_update = landscape
    
    ## New vegetation
    
    ## XXXX problem here : how to update the cells ???
    
    for(i in 1:S){
      for (k in 1:S){
        landscape_update[which(landscape == k & rnum <= colonization[[i]])] = 1
        
        
      }
      
      
      
    }
    
    ## New fertile
    landscape_update[which(landscape == -1 & rnum <= death)] = 0
    
    
    landscape_update[which(landscape == -1 & rnum <= regeneration)] = 0
    
    ## New degraded 
    landscape_update[which(landscape == 0 & rnum > total_colonization & rnum <= total_colonization + degradation)] = -1
    
    
    return(   list(
                   Rho_v = sum(landscape_update == 1) / length(landscape_update),
                   Rho_f = sum(landscape_update == 2) / length(landscape_update),
                   Rho_D = 1-sum(landscape_update == 1) / length(landscape_update)-sum(landscape_update == 2) / length(landscape_update),
                   Landscape=landscape_update))
  })
  
}

# Ploting functions ----

Run_CA=function(time=seq(1,1000,1),params,frac=c(.8,.1,.1),size=25,plot=F){
  
  ini=Get_initial_lattice(frac = frac,size=size)
  
  d=tibble(Time=1,Rho_V=sum(ini == 1) / length(ini),Rho_F=sum(ini == 2) / length(ini),Rho_D=sum(ini == 3) / length(ini))
  state=list(Landscape=ini,Rho_v=d$Rho_V,Rho_f=d$Rho_F,Rho_D=d$Rho_D)
  
  for (k in 2:length(time)){
    
    params$dt=time[k]-time[k-1]
    
    state=CA_network(state$Landscape,params = params)

    d=rbind(d,tibble(Time=k,Rho_V=state$Rho_v,Rho_F=state$Rho_f,Rho_D=state$Rho_D))
  }
  
  return(list(d=d,State=state$Landscape))
  
}

Plot_dynamics=function(d,different_sim=F,simple=F,params){
  col_pal=colorRampPalette(c("#216707","#359E0F","#68C745","#9BD088","#D2E6CB"))
  
  if (different_sim==T & simple==F){  
    d_melt=melt(d,id.vars=c("Time","Type"))
    p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_V"="Vegetation","Rho_F"="Fertile","Rho_D"="Degraded")))+
      geom_line(aes(x=Time,y=value,color=variable,linetype=Type),lwd=1)+scale_color_manual(values=rev(c(col_pal(params$S),"#D2C48F","#777777") ))+
      labs(x="Time steps",y="Fraction",color="")+the_theme
    return(p)
    
  } 
  
  if (different_sim==F & simple==F){  
    
    d_melt=melt(d,id.vars=c("Time"))
    p=ggplot(d_melt%>%mutate(.,variable=recode_factor(variable,"Rho_F"="Fertile","Rho_D"="Degraded")))+
      geom_line(aes(x=Time,y=value,color=variable))+scale_color_manual(values=rev(c(col_pal(params$S),"#D2C48F","#777777") ))+
      labs(x="Time steps",y="Fraction",color="")+the_theme
    return(p)
    
  }
  
  if (simple == T){
    plot(d$Time,d$Rho_V,xlab="Time steps","l",ylab="Fraction",ylim=c(0,1),col="#9DD696")
    lines(d$Time,d$Rho_F,ylim=c(0,1),col="#D2C48F")
    lines(d$Time,d$Rho_D,ylim=c(0,1),col="#777777")
  }
  
}

Plot_landscape=function(landscape,params){
  col_pal=colorRampPalette(c("#216707","#359E0F","#68C745","#9BD088","#D2E6CB"))
  image(landscape,xaxt = "n",yaxt ="n",col=rev(c(col_pal(params$S),"#D2C48F","#777777") ))
}


Plot_network=function(params){
  
  graph_net=graph_from_adjacency_matrix(params$A,weighted = T)
  color_link=sapply(1:length(E(graph_net)$weight),function(x){ifelse(E(graph_net)$weight[x]>0,"blue","red")})
  plot(graph_net,edge.color=color_link, vertex.color = col.alpha("forestgreen",.5),
       vertex.size=30,vertex.label="",layout=layout_in_circle,
       edge.width=3, edge.arrow.size=0.1,edge.arrow.width=1
  )
  
}


Plot_species_niche=function(params){
  d=tibble()
  for (i in 1:params$S){
    d=rbind(d,tibble(Species=as.factor(i),Niche=exp(-((seq(0,1,length.out=100)-params$b_opt[i])/params$Tol_sp[i])**2),
                     Aridity=seq(0,1,length.out=100)))
  }
  print(ggplot(d%>%melt(., measure.vars="Niche"))+geom_line(aes(x=Aridity,y=value,color=Species,group=Species),lwd=.8)+theme_classic()+
    theme(legend.position = "none")+scale_color_brewer(palette="Dark2")+geom_area(aes(x=Aridity,y=value,fill=Species),position = 'identity',alpha=.2)+
      scale_fill_brewer(palette="Dark2"))
  
}


# Mean field fonction ----

Mean_field_network=function(t,state,param){
  
  with (as.list(c(param)),{
    
    state[state<10^-8]=0 #prevent numerical problems
    
    
    rho_V=sum(state[1:(length(state)-2)]) #fraction vegetation
    rho_D=state[length(state)] #fraction degraded
    rho_S=state[(length(state)-1)] #fraction fertile
    
    dstate_dt=matrix(0,nrow=S+2)
    
    for (sp in 1:S){
      dstate_dt[sp,] = state[sp]*(rho_S*(b-c*rho_V)-
                          m0*(1- sum(A[sp,-sp]*state[1:S][-sp])/(1+abs(sum(A[sp,-sp]*state[1:S][-sp])))) * # interactions
                      (1-teta*exp(((-(b-b_opt[sp])/Tol_sp[sp])**2))) )                         # aridity
      
    }  
    

    dstate_dt[length(dstate_dt)-1,]=(r+f*rho_V)*rho_D  - #regeneration
      d*rho_S                   + #degradation
      sum(sapply(1:S,function(x){
        state[x] *(m0*(1- sum(A[x,-x]*state[1:S][-x])/(1+abs(sum(A[x,-x]*state[1:S][-x])))) * # interactions
                      (1-teta*exp(((-(b-b_opt[x])/Tol_sp[x])**2))) )                         # aridity
        
      }))  - rho_S*rho_V*(b-c*rho_V)  

    
    dstate_dt[length(dstate_dt),]=d*rho_S  - (r+f*rho_V)*rho_D  
    
    list(dstate_dt)
  })
  
    
}




Compute_ode=function(state,params,method_ode="lsoda",optim_time=T,h=.1,n_time_no_opt=2000){
  
  if ( optim_time==T){
    n_time_ode=0
    count_non_eq=1
    
    while(count_non_eq!=0 & n_time_ode<5000){ #while we are not at equilibrium for all the organisms/nutrients/detritus
      
      n_time_ode=n_time_ode+500
      time=seq(0,n_time_ode,h)
      count_non_eq=0
      dynamics = as.data.frame(ode(state,time,func=Mean_field_network,parms=params ,method = method_ode))
    
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
    
  }else{
    n_time_ode=n_time_no_opt
    time=seq(0,n_time_ode,.1)
    dynamics = as.data.frame(ode(state,time,func=Mean_field_network,parms=params ,method = method_ode))
    
  } 
  
  final_point=as.numeric(dynamics[nrow(dynamics),-1])
  dynamics[nrow(dynamics),-1]=final_point
  colnames(dynamics)=c("Time",paste0("Rho_V",1:params$S),"Rho_S","Rho_D")
  
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




