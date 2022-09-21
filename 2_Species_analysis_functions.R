x=c("tidyverse","ggpubr","latex2exp","deSolve","reshape2","JuliaCall","diffeqr","simecol","tseries","phaseR","ggquiver","scales")
lapply(x, require, character.only = TRUE)


the_theme=theme_classic()+theme(legend.position = "bottom",
                                strip.background = element_rect(fill = "#CCE8D8"),
                                strip.text.y = element_text(size = 10, angle = -90),
                                strip.text.x = element_text(size = 8),axis.text = element_text(size=11),axis.title = element_text(size=13),
                                legend.text = element_text(size = 10),text = element_text(family = "NewCenturySchoolbook"))


color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")

# A) Mean field analysis ----

Get_MF_parameters=function(){
  
  return(c(
    r = 0.05, d = 0.1, f = 0.9, beta = .8, m = 0.1, e = .1, emax = 1.2, cintra = .1,
    cinter1 = .1, cinter2 = .1, S = 0
  )
  
  )
}

Get_MF_initial_state=function(){
  state <- c(rho_1 = 0.4, rho_2 = 0.4, rho_m = 0.1,rho_0=.1)
  names(state) <- c("rho_1", "rho_2", "rho_0")
  return(state)
}



#normal function
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


#function of press perturbation
MF_two_species_julia_press = julia_eval(" 

function MF_two_species_press(du,u,p,t)

  r,d,f,beta1,beta2,m,e,emax,cintra,cinter1,cinter2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta1 * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2))) - rho_1 * m
  du[2] = rho_0 * ( beta2 * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta1 * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2)))  -
      rho_0 * ( beta2 * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) + rho_2 * m + rho_1 * m

end")




MF_two_species_julia_SGH = julia_eval(" 

function MF_two_species_SGH(du,u,p,t)

  r,d,f,beta1,beta2,m,e,emax,cintra,cinter1,cinter2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+S*f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+S*f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) + rho_2 * m + rho_1 * m

end")


MF_two_species_julia_SGH_press = julia_eval(" 

function MF_two_species_SGH_press(du,u,p,t)

  r,d,f,beta1,beta2,m,e,emax,cintra,cinter1,cinter2,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta1 * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2))) - rho_1 * m
  du[2] = rho_0 * ( beta2 * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+S*f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+S*f*(rho_1))-rho_0 * ( beta1 * rho_1 *  (emax *( 1 -S*(1-e)) - (cintra*rho_1 + cinter1*rho_2)))  -
      rho_0 * ( beta2 * rho_2 *  (emax *( 1 -S) - (cintra*rho_2 + cinter2*rho_1))) + rho_2 * m + rho_1 * m

end")


#Analysis and ploting figs

Post_processing_MF=function(d2,name,C_seq=c_seq){
  
  # write.table(d2,paste0("../Table/2_species/2_species_",name,".csv"),sep=";")
  

  d2$state=sapply(1:nrow(d2),function(x){
    
    if (d2[x,1]>0 & d2[x,2]>0) return("coexistence")
    if (d2[x,1]>0 & d2[x,2]==0) return("stress_tol")
    if (d2[x,1]==0 & d2[x,2]>0) return("competitive")
    if (d2[x,1]==0 & d2[x,2]==0) return('desert')
    
  })
  
  c_values_bifu=C_seq[c(1,round(length(C_seq)/2),length(C_seq))]
  
  
  color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  #state at equilibrium
  p1=ggplot(d2)+geom_tile(aes(x=S,y=cinter1/param["cinter2"],fill=state))+
    theme_classic()+scale_fill_manual(values=color_rho)+
    annotate("text",x=rep(1.02,3),y= (c_values_bifu/param["cinter2"])+0.1,label=c("A","B","C"),color="black")+
    theme(legend.position = "bottom")+labs(x="Stress (S)",y=TeX(r'(Relative competitive ability \ $\frac{c_{2,1}}{c_{1,2}})'),fill="")+
    geom_hline(yintercept = c_values_bifu/param["cinter2"],lwd=.1,color="gray40")+
    theme(legend.text = element_text(size=11))
  
  
  #density of global vegetation
  density_col=colorRampPalette(c("red","white","blue"))
  p2= ggplot(d2)+geom_tile(aes(x=S,y=cinter1/param["cinter2"],fill=rho_plus))+
    theme_classic()+scale_fill_gradientn(colours=density_col(100))+
    theme(legend.position = "bottom")+labs(x="Stress (S)",y=TeX(r'(Relative competitive ability \ $\frac{c_{2,1}}{c_{1,2}})'),
                                           fill=TeX(r'(Global vegetation \ $\rho_1 + \rho_2$)'))
  
  
  #some bifurcation diagrams
  
  c_values_bifu=c_seq[c(1,round(length(c_seq)/2),length(c_seq))]
  d_bifu=filter(d2,round(cinter1,4) %in% round(c_values_bifu,4))
  for (c_bifu in 1:3){
    assign(paste0("p3_",c_bifu),
           ggplot(d_bifu %>% filter(.,round(cinter1,4)==round(c_values_bifu[c_bifu],4)) %>% melt(.,measure.vars=c("rho_1","rho_2"))%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive")))+
             geom_point(aes(x=S,y=value,color=variable),size=.5)+labs(x="Stress (S)",y="Density",color="")+
             the_theme+ scale_color_manual(values=color_rho[c(2,4)])+  theme(legend.text = element_text(size=12))
    )
  }
  p3=ggarrange(p3_3,p3_2,p3_1,nrow=3,common.legend = T,legend = "bottom",labels = LETTERS[1:3])
  
  
  p_tile=ggarrange(p1,p2,ncol=2)
  p_tot=ggarrange(p_tile,p3,ncol=2,widths = c(3,1))
  ggsave(paste0("../Figures/2_species/2_species_",name,".pdf"),width = 14,height = 6)  
  
}




# Hysteresis functions 


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
      if (max(abs(diff(biomass_D[,1])))>.1){
        # #now tipping points
        Tipping_D=biomass_D$S[which(abs(diff(biomass_D[,1]))==max(abs(diff(biomass_D[,1]))))]
        Tipping_R=biomass_R$S[which(abs(diff(biomass_R[,1]))==max(abs(diff(biomass_R[,1]))))]
        
      } else{
        Tipping_D=0
        Tipping_R=0
      }
      
      d_hysteresis=rbind(d_hysteresis,tibble(Hysteresis=hysteresis,Species=sp,Tipping_D=Tipping_D,Tipping_R=Tipping_R))
      
      
    }
  }
  return(d_hysteresis)
}



plot_dynamics=function(d){
  
  color_rho=c("fertile"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  colnames(d)=c("rho_1",'rho_2',"rho_0",'rho_d')
  if ('time' %in% colnames(d) | 'Time' %in% colnames(d)){
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)])%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive","rho_0"="fertile","rho_d"="desert")))+
             geom_line(aes(x=time,y=value,color=variable),lwd=1)+
             theme_classic()+scale_color_manual(values=color_rho)+labs(x="Time",y="Densities",color="")+
             theme(legend.text = element_text(size=11),legend.position = "bottom"))
  }
  else {
    d$time=1:nrow(d)
    return(ggplot(d%>%melt(.,id.vars=colnames(d)[ncol(d)])%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive","rho_0"="fertile","rho_d"="desert")))+
             geom_line(aes(x=time,y=value,color=variable),lwd=1)+
             the_theme+scale_color_manual(values=color_rho)+labs(x="Time",y="Densities",color="")+
             theme(legend.text = element_text(size=11),legend.position = "bottom"))
  }
}


# Recruitment rate

Get_recruitment_rates=function(eq,param){
  
  recruit_rate = c(
    r1 =  (1-eq$rho_1-eq$rho_2-eq$rho_d) * param["beta"]  *  eq$rho_1  *   (param["emax"] *( 1 -param["S"]  * (1-param["e"])) - (param["cintra"]*eq$rho_1 + param["cinter1"]*eq$rho_2)) - eq$rho_1*param["m"],
    r2 =  (1-eq$rho_1-eq$rho_2-eq$rho_d) * param["beta"]  *  eq$rho_2  *   (param["emax"] *( 1 -param["S"])                   - (param["cintra"]*eq$rho_2 + param["cinter2"]*eq$rho_1)) - eq$rho_2*param["m"]  )
  
  names(recruit_rate)=c("r1","r2")
  return(recruit_rate)

}


# B) Pair approximation analysis ----

Get_PA_parameters=function(){
  
  return(c(    r = 0.05, d = 0.1, f = .9, beta = 0.8, m = 0.1, e = .9, emax = 1.2, cg = .1, alpha11 = .1,
    alpha12 = .1, alpha21 = .1, alpha22 = .1, S = 0, delta = .1, z = 4
  )
  )
}

Get_PA_initial_state=function(){
  state <- c(rho_1 = 0.4, rho_2 = 0.4, rho_m = 0.1)
  state_pair <- c(
    rho_12 = state[1] * state[2], rho_1m = state[1] * state[3], rho_2m = state[2] * state[3],
    rho_11 = state[1] * state[1], rho_22 = state[2] * state[2], rho_mm = state[3] * state[3]
  )
  state <- c(state, state_pair)
  names(state) <- c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
  return(state)
}

get_mean_densities=function(d){
  
  length_transient=3000
  d=d[(length_transient+1):nrow(d),]
  return(colMeans(d))
  
}


PA_two_species_julia = julia_eval("

function PA_two_species(du,u,p,t)

r,d,f,beta,m,e,emax,cg,alpha11,alpha12,alpha21,alpha22,S,delta,z=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

#rho_1
du[1] =  (1- rho_1-rho_2-rho_m) * beta * (delta * rho_1 + (1 - delta) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) + alpha21 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) )) - rho_1 * m

#rho_2
du[2] =  (1- rho_1-rho_2-rho_m) * beta * (delta * rho_2 + (1 - delta) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S)  - (cg *(rho_1 + rho_2) + (alpha22 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) + alpha12 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) )) - rho_2 * m

#rho_m
du[3] =  (1- rho_1-rho_2-rho_m) * d - rho_m * (r + f * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m)))

#rho_12
du[4] = ((rho_1 - rho_11 - rho_12 - rho_1m)) * beta * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S ) - (cg *(rho_1 + rho_2) + (alpha22 *((z - 1) / z)*((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) + (alpha12/z) + alpha12 *((z - 1) / z)*((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) ))  +
        ((rho_2 - rho_22 - rho_12 - rho_2m)) * beta * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) + (alpha21/z) + alpha21 * ((z - 1) / z) *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) ))-
        2 * rho_12 * m

#rho_1m
du[5] = ((rho_1 - rho_11 - rho_12 - rho_1m)) * d + (rho_m - rho_mm - rho_1m - rho_2m) * beta * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) + alpha21 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) )) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m)) ))

#rho_2m
du[6] = ((rho_2 - rho_22 - rho_12 - rho_2m)) * d + (rho_m - rho_mm - rho_1m - rho_2m) * beta * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m)) * (emax * (1 - S  ) - (cg *(rho_1 + rho_2) + (alpha22 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) + alpha12 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) ))) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m) )))

#rho_11
du[7] = 2* ((rho_1 - rho_11 - rho_12 - rho_1m)) *  beta * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) * ((z - 1) / z) + (alpha11/z) + alpha21 * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) )) -
        2 * rho_11 * m

#rho_22
du[8] = 2* ((rho_2 - rho_22 - rho_12 - rho_2m)) * beta * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S ) - (cg *(rho_1 + rho_2) + (alpha22 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) * ((z - 1) / z) + (alpha22/z) + alpha12 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) )) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_m - rho_mm - rho_1m - rho_2m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m)))
  
end
  

")

PA_two_species_julia_press = julia_eval("

function PA_two_species_press(du,u,p,t)

r,d,f,beta1,beta2,m,e,emax,cg,alpha11,alpha12,alpha21,alpha22,S,delta,z=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

#rho_1
du[1] =  (1- rho_1-rho_2-rho_m) * beta1 * (delta * rho_1 + (1 - delta) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) + alpha21 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) )) - rho_1 * m

#rho_2
du[2] =  (1- rho_1-rho_2-rho_m) * beta2 * (delta * rho_2 + (1 - delta) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S)  - (cg *(rho_1 + rho_2) + (alpha22 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) + alpha12 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) )) - rho_2 * m

#rho_m
du[3] =  (1- rho_1-rho_2-rho_m) * d - rho_m * (r + f * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m)))

#rho_12
du[4] = ((rho_1 - rho_11 - rho_12 - rho_1m)) * beta2 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S ) - (cg *(rho_1 + rho_2) + (alpha22 *((z - 1) / z)*((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) + (alpha12/z) + alpha12 *((z - 1) / z)*((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) ))  +
        ((rho_2 - rho_22 - rho_12 - rho_2m)) * beta1 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) + (alpha21/z) + alpha21 * ((z - 1) / z) *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) ))-
        2 * rho_12 * m

#rho_1m
du[5] = ((rho_1 - rho_11 - rho_12 - rho_1m)) * d + (rho_m - rho_mm - rho_1m - rho_2m) * beta1 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m))) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) + alpha21 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) )) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m)) ))

#rho_2m
du[6] = ((rho_2 - rho_22 - rho_12 - rho_2m)) * d + (rho_m - rho_mm - rho_1m - rho_2m) * beta2 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m)) * (emax * (1 - S  ) - (cg *(rho_1 + rho_2) + (alpha22 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) + alpha12 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) ))) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m) )))

#rho_11
du[7] = 2* ((rho_1 - rho_11 - rho_12 - rho_1m)) *  beta1 * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S * e ) - (cg *(rho_1 + rho_2) + (alpha11 *((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m)) * ((z - 1) / z) + (alpha11/z) + alpha21 * ((z - 1) / z) * ((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m))) )) -
        2 * rho_11 * m

#rho_22
du[8] = 2* ((rho_2 - rho_22 - rho_12 - rho_2m)) * beta2 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (((rho_2 - rho_22 - rho_12 - rho_2m)) / (1- rho_1-rho_2-rho_m) )) * (emax * (1 - S ) - (cg *(rho_1 + rho_2) + (alpha22 *((rho_2 - rho_22 - rho_12 - rho_2m)/(1- rho_1-rho_2-rho_m)) * ((z - 1) / z) + (alpha22/z) + alpha12 * ((z - 1) / z) * ((rho_1 - rho_11 - rho_12 - rho_1m)/(1- rho_1-rho_2-rho_m))) )) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_m - rho_mm - rho_1m - rho_2m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (((rho_1 - rho_11 - rho_12 - rho_1m)) / (1- rho_1-rho_2-rho_m)))
  
end
  

")


post_processing_2species_PA=function(d2,Alpha_seq=alpha_seq,S_seq=s_seq,title=F,title_name=""){
  
  
  d2$state=sapply(1:nrow(d2),function(x){
    
    if (d2$rho_1[x]>0 & d2$rho_2[x]>0) return("coexistence")
    if (d2$rho_1[x]>0 & d2$rho_2[x]==0) return("stress_tol")
    if (d2$rho_1[x]==0 & d2$rho_2[x]>0) return("competitive")
    if (d2$rho_1[x]==0 & d2$rho_2[x]==0) return('desert')
    
  })
  
  c_values_bifu=Alpha_seq[c(1,round(length(Alpha_seq)/2),length(Alpha_seq))]
  
  
  color_rho=c("coexistence"="#D8CC7B","competitive"="#ACD87B","desert"="#696969","stress_tol"="#7BD8D3")
  
  #state at equilibrium
  p1=ggplot(d2)+geom_tile(aes(x=S,y=alpha21/param["alpha12"],fill=state))+
    theme_classic()+scale_fill_manual(values=color_rho)+
    theme(legend.position = "bottom")+labs(x="Stress (S)",y=TeX(r'(Relative competitive ability \ $\frac{\alpha_{2,1}}{\alpha_{1,2}})'),fill="")+
    geom_hline(yintercept = c_values_bifu/param["alpha12"],lwd=.1,color="gray40")+
    theme(legend.text = element_text(size=11))
  
  
  # #density of global vegetation
  # density_col=colorRampPalette(c("red","white","blue"))
  # p2= ggplot(d2)+geom_tile(aes(x=S,y=alpha21/param["alpha12"],fill=rho_plus))+
  #   theme_classic()+scale_fill_gradientn(colours=density_col(100))+
  #   theme(legend.position = "bottom")+labs(x="Stress (S)",y=TeX(r'(Relative competitive ability \ $\frac{\alpha_{2,1}}{\alpha_{1,2}})'),fill="Global vegetation density")
  
  #pair 12
  density_pair=colorRampPalette(c("yellow","#CE7604"))
  d2$rho_12[which(d2$rho_1==0 | d2$rho_2==0)]=NA
  p4= ggplot(d2)+geom_tile(aes(x=S,y=alpha21/param["alpha12"],fill=rho_12))+
    theme_classic()+scale_fill_gradientn(colours=density_pair(100),na.value = "grey")+
    theme(legend.position = "bottom")+labs(x="Stress (S)",y=TeX(r'(Relative competitive ability \ $\frac{\alpha_{2,1}}{\alpha_{1,2}})'),fill="Pair sp1 & sp2")
  
  #some bifurcation diagrams
  
  a_values_bifu=Alpha_seq[c(1,round(length(Alpha_seq)/2),length(Alpha_seq))]
  d_bifu=filter(d2,round(alpha21,4) %in% round(a_values_bifu,4))
  for (a_bifu in 1:3){
    assign(paste0("p3_",a_bifu),
           ggplot(d_bifu %>% filter(.,round(alpha21,4)==round(a_values_bifu[a_bifu],4)) %>% melt(.,measure.vars=c("rho_1","rho_2"))%>%
                    mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive")))+
             geom_point(aes(x=S,y=value,color=variable),size=.75)+labs(x="Stress (S)",y="Density",color="")+
             the_theme+ scale_color_manual(values=color_rho[c(2,4)])+  theme(legend.text = element_text(size=11))
    )
  }
  p3=ggarrange(p3_3,p3_2,p3_1,nrow=3,common.legend = T,legend = "bottom")
  if (title){
    p1=p1+ggtitle(title_name)
  }
  
  p_tile=ggarrange(p1,p4,ncol=2)
  p_tot=ggarrange(p_tile,p3,ncol=2,widths = c(3,1))
  ggsave("../Figures/2_species/PA.pdf",width = 14,height = 6)
  return(p_tot)
}


# C) CA analysis ----


CA_2_species=function(landscape, param){
  
  # Variables : 1 = stress_tol, 2 = competitive, 0 = fertile, -1 = degraded   
  rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
  rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
  rho_0 = length(which((landscape == 0)))   / length(landscape) #fraction fertile
  rho_d = 1-rho_1-rho_2-rho_0 # fraction degraded
  
  # Neighbors :
  
  #using simcol package  
  neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
  neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
  neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
  neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
  
  
  r=param[1];d=param[2];f=param[3];beta=param[4];m=param[5];e=param[6]
  emax=param[7];cintra=param[8];cinter1=param[9];cinter2=param[10];S=param[11];delta=param[12]
  z=param[13]
  
  colonization_1 = beta*(delta * rho_1 + (1 - delta) * neigh_1 / z) * (emax* ( 1 - S * (1-e)) - (cintra*rho_1+cinter1*rho_2) ) 
  colonization_2 = beta*(delta * rho_2 + (1 - delta) * neigh_2 / z) * (emax* ( 1 - S ) - (cintra*rho_2+cinter2*rho_1) ) 
  
  
  # calculate regeneration, degradation & mortality rate
  death = m 
  regeneration = (r + f * (neigh_1) / z)
  degradation = d
  
  # Apply rules
  rnum = runif(length(landscape))# one random number between 0 and 1 for each cell
  landscape_update=landscape
  
  ## New stress_tol
  landscape_update[which(landscape ==0 & rnum <= colonization_1)] = 1
  
  ## New competitive
  landscape_update[which(landscape ==0 & rnum > colonization_1 & rnum <= colonization_2 + colonization_1)] = 2
  
  ## New fertile
  landscape_update[which(landscape == 1 & rnum <= death)] = 0
  landscape_update[which(landscape == 2 & rnum <= death)] = 0
  landscape_update[which(landscape == -1 & rnum <= regeneration)] = 0
  
  ## New degraded 
  landscape_update[which(landscape == 0 & rnum > colonization_2 + colonization_1 & rnum <= colonization_2 + colonization_1 + degradation)] = -1
  
  rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
  rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
  rho_0 = length(which((landscape == 0)))    / length(landscape) #fraction fertile
  rho_d = 1-rho_1-rho_2-rho_0 # fraction degraded
  
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
  
  
  return(list(State=c(rho_1, rho_2, rho_0, rho_d),Landscape=landscape_update))
  
}



Run_CA_2_species=function(time=seq(1,1000,1),param,landscape){
  
  d=tibble(rho_1=sum(landscape == 1) / length(landscape),rho_2=sum(landscape == 2) / length(landscape),
           rho_0=sum(landscape == 0) / length(landscape),rho_d=sum(landscape == -1) / length(landscape),time=1)
  
  for (k in 2:length(time)){
    
    param["dt"]=time[k]-time[k-1]
    output_CA=CA_2_species(landscape,param = param)
    landscape=output_CA$Landscape
    d=rbind(d,tibble(rho_1=output_CA$State[1],rho_2=output_CA$State[2],
                     rho_0=output_CA$State[3],rho_d=output_CA$State[4],
                     time=k))
    
  }
  
  return(list(state=d,landscape=landscape))
  
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







# D) Gillespie ----


Gillespie_tau_leaping_R=function(landscape,time,param){
  
  d2=tibble()
  #param
  r=param[1];d=param[2];f=param[3];beta=param[4];m=param[5];e=param[6]
  emax=param[7];cintra=param[8];cinter2=param[9];cinter1=param[10];S=param[11];delta=param[12]
  z=param[13];leap=param[14]
  
  Rate_landscape=array(0,c(nrow(landscape),ncol(landscape),6))
  
  rules_change=matrix(c(0,1,0,2,0,-1,1,0,2,0,-1,0),ncol=2,nrow=6,byrow = T)

  
  for (dt in 1:length(time)){ #for each time step
    
    #calculate all rates in the lattice
    
    #number of neighbors
    neigh_1= simecol::neighbors(x =landscape,state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_2= simecol::neighbors(x =landscape,state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_f= simecol::neighbors(x =landscape,state = 0, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    neigh_d= simecol::neighbors(x =landscape,state = -1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
    
    #global densities
    rho_1 = length(which((landscape == 1)))   / length(landscape) #fraction stress_tol
    rho_2 = length(which((landscape == 2)))   / length(landscape) #fraction competitive
    rho_0 = length(which((landscape == 0)))   / length(landscape) #fraction fertile
    rho_d = 1-rho_1-rho_2-rho_0 # fraction degraded
    
    # calculate the rates
    Rate_landscape[,,1] = beta*(delta * rho_1 + (1 - delta) * neigh_1 / z) * (emax* ( 1 - S * (1-e)) - (cintra*rho_1+cinter1*rho_2) ) 
    Rate_landscape[,,2] = beta*(delta * rho_2 + (1 - delta) * neigh_2 / z) * (emax* ( 1 - S ) - (cintra*rho_2+cinter2*rho_1) ) 
    Rate_landscape[,,3] = d
    Rate_landscape[,,4] = m
    Rate_landscape[,,5] = m
    Rate_landscape[,,6] = (r + f * neigh_1 / z)
    
    #calculate their propensity (i.e. sum of type of events)
    Propensity = c(sum(Rate_landscape[,,1][which(landscape==0)]), #colonization stress tol
                   sum(Rate_landscape[,,2][which(landscape==0)]), #colonization competitive
                   sum(Rate_landscape[,,3][which(landscape==0)]), #degradation
                   sum(Rate_landscape[,,4][which(landscape==1)]), #death stress_tol
                   sum(Rate_landscape[,,5][which(landscape==2)]), #death competitive
                   sum(Rate_landscape[,,6][which(landscape==-1)]) #regeneration
    )
    
    nb_events=sapply(1:length(Propensity),function(x){
      rpois(1,lambda = Propensity[x]*leap)
    }) #number of events per event type
    
    
    for (ev in 1:length(nb_events)){ #for each type of events
      
      patches=which(landscape==rules_change[ev,1])
      
      if (nb_events[ev] != 0 & length(patches) >= nb_events[ev]){
        
        if (length(Rate_landscape[,,ev][patches])==1){
          proba_sample=rep(Rate_landscape[,,ev][patches],length(patches))
        } else {
          proba_sample=Rate_landscape[,,ev][patches]
        }
        
        landscape[sample( patches,
                          #in prob we account that mortality is constant thus, we repeat the same proba to give to sample
                          prob = ,
                          nb_events[ev],replace = T)] = rules_change[ev,2]
      }
      
    }
    
    
    d2=rbind(d2,tibble(  rho_1 = length(which((landscape == 1)))   / length(landscape), #fraction stress_tol
                       rho_2 = length(which((landscape == 2)))   / length(landscape), #fraction competitive
                       rho_0 = length(which((landscape == 0)))   / length(landscape), #fraction fertile
                       rho_d = 1-rho_1-rho_2-rho_0 # fraction degraded
    ))
  }

  return(list(state=d2,landscape=landscape))

}
