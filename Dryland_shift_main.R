# Other) Recruitment rate & 4 grid matrices----
rm(list = ls())
source("./Dryland_shift_functions.R")
d = expand_grid(emax = 1, e = .1, S = seq(0, 1, length.out = 200), psi = c(0, .5, 1))

p = ggplot(d) +
    geom_line(aes(x = S, y = emax * (1 - S * (1 - e * psi)), color = as.factor(psi), group = psi),size=1) +
    scale_color_manual(values = color_Nsp(9)[c(2,5,8)]) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
        x = "Stress, S",
        y = TeX(r'($Recruitment rate,  (1-S (1-e \psi_i)))'), color = TeX("$\\psi_i$")
    )



ggsave("../Figures/Recruitment_rate.pdf", p, width = 6, height = 4)



# Step 1) Two species MF model ----
rm(list = ls())
source("./Dryland_shift_functions.R")
julia_setup()
de = diffeq_setup()
#





## 1) Multistability fixed traits ----



tspan = c(0, 5000)
t = seq(0, 5000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.4,length.out=N_sim)
branches=c("Degradation","Restoration")

for (branch in branches){
  
  if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
    state =Get_MF_initial_state(c(.4,.4,.1))
    S_seq=seq(0,1, length.out = N_sim)
  }else {
    state =Get_MF_initial_state(c(.005,.005,.49))
    S_seq=rev(seq(0,1, length.out = N_sim))
  }
  
  
  for (ccomp in c_seq){
    
    for (S in S_seq) { 
      
      julia_assign("state", state)
      param=Get_MF_parameters()
      param["cintra"]=.3
      param["S"] = S
      param["alpha_0"] = ccomp
      julia_assign("p", param)
      
      
      
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      
      
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = ccomp,
                                                 Branch=branch))
      
    }
  }
}
d2[d2 < 10^-4] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2
d2=d2[,-c(3:4)]
colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","Branches","Rho_plus")
write.table(d2,paste0("../Table/2_species/MF/Multistability_MF.csv"),sep=";")



## 2) Multistability trait difference ----



tspan = c(0, 3000)
t = seq(0,3000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=50
S_seq =seq(0,.82, length.out = 100)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=.3
psi1_seq=seq(0,1,length.out=N_sim)
f_seq=c(.9)
dispersal_scale=c(.1)
branches=c("Degradation","Restoration")

for (facil in f_seq){
  
  for (disp in dispersal_scale){
    
    for (S in S_seq) { 
      
      
      for (cinter in c_inter_seq){
        
        
        
        for (branch in branches){
          
          if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
            state = Get_MF_initial_state(c(.4,.4,.1))
          }else {
            state = Get_MF_initial_state(c(.005,.005,.49))
          }
          
          d2 = tibble()
          
          
          for (psi2 in psi_seq){
            
            for (psi1 in psi1_seq[which(psi1_seq>psi2)]){
              
              julia_assign("state", state)
              param=Get_MF_parameters()
              param["cintra"]=.3
              param["alpha_0"]=cinter
              param["S"] = S
              param["psi_1"]=psi1
              param["psi_2"]=psi2
              param["f"]=facil
              julia_assign("p", param)
              
              prob = julia_eval("ODEProblem(MF_two_species_varying_trait, state, tspan, p)")
              
              
              sol = de$solve(prob, de$Tsit5(), saveat = t)
              d = as.data.frame(t(sapply(sol$u, identity)))
              colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_0")
              
              d2 = rbind(d2, d[nrow(d), ] %>% add_column(Stress = S, Psi2 = psi2, Psi1 = psi1,alpha_0=cinter,
                                                         Branch=branch))
              
            }
          } # end trait value 2nd species
          
          d2[d2 < 10^-4] = 0
          d2$rho_plus = d2$rho_1 + d2$rho_2
          colnames(d2) = c("Sp1", "Sp2","Degraded","Fertile", "Stress", "Psi2","Psi1","alpha_0","Branches","Rho_plus")
          write.table(d2,paste0("../Table/2_species/MF/Test_interspe_comp_",
                                cinter,"_branch_",branch,
                                "_stress_",S,"_facilitation_",facil,".csv"),sep=";")
          
        } # end loop interspecific competition
        
      } # end loop branch
      
    } # end loop first species trait
    
  } #end loop dispersal
  
}#end facilitation loop






# Step 2) Pair approximation (PA) ----
rm(list = ls())
source("./Dryland_shift_functions.R")
julia_setup()
de = diffeq_setup()
#



## 1) Multistability fixed traits ----



tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.4,length.out=N_sim)
name_scena=c("local_C_local_F","global_C_global_F","local_C_global_F","global_C_local_F")
delta_seq=c(.1,.9)
branches=c("Degradation","Restoration")

for (branch in branches){
  
  if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
    state =Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
    S_seq=seq(0,1, length.out = N_sim)
  }else {
    state =Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
    S_seq=rev(seq(0,1, length.out = N_sim))
  }
  
  for (disp in delta_seq) {
    
    for (scena_ID in 1:4){ #for each scenario of species pairs
      
      for (ccomp in c_seq){
        
        for (S in S_seq) { #varying dispersal scale
          
          julia_assign("state", state)
          param=Get_PA_parameters()
          param["cintra"]=.3
          param["S"] = S
          param["alpha_0"] = ccomp
          param["delta"]=disp
          julia_assign("p", param)
          
          if (scena_ID==1){ 
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_local_F, state, tspan, p)")
            
          }else if (scena_ID==2){  #global C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==3){  #local C, global F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_local_C_global_F, state, tspan, p)")
            
          }else if (scena_ID==4){ #global C, local F
            
            julia_assign("p", param)
            prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
            
          }
          
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
          
          d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, alpha_0 = ccomp,
                                                     Scena=name_scena[scena_ID],
                                                     Delta=disp,Branch=branch))
          
        }
      }
      
      d2[d2 < 10^-4] = 0
      d2$rho_plus = d2$rho_1 + d2$rho_2
      d2=d2[,c(1:2,10:15)]
      colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "alpha_0","Scena","Delta","Branches","Rho_plus")
      write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Fixed_traits/Multistability_mapping_Scena_",
                            name_scena[scena_ID],"_branch_",branch,"_disp_",disp,".csv"),sep=";")
      d2=tibble()
    }
  }
}



d2 = tibble()
N_sim=100
S_seq = seq(0,1, length.out = N_sim)
c_seq=seq(0,.4,length.out=N_sim)
name_scena=c("global_C_local_F","global_C_global_F","local_C_local_F","local_C_global_F")
delta_seq=c(.1,.9)
branches=c("Degradation","Restoration")

for (branch in branches){
  
  for (disp in delta_seq) {
    
    for (scena_ID in 1:4){ #for each scenario of species pairs
      
      d2=rbind(d2,read.table(paste0("../Table/2_species/PA/Multistability_PA/Fixed_traits/Multistability_mapping_Scena_",name_scena[scena_ID],
                                    "_branch_",branch,"_disp_",disp,".csv"),sep=";"))
      
    }
  }
}
write.table(d2,paste0("../Table/2_species/PA/Multistability_fixed_traits_PA.csv"),sep=";")





## 2) Multistability and trait difference ----


tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=50
S_seq =seq(0,.82, length.out = 100)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=.3
psi1_seq=seq(0,1,length.out=N_sim)
f_seq=c(.9)
dispersal_scale=c(.1,.9)
branches=c("Degradation","Restoration")

for(scale_f in c("local","global")){

  for (facil in f_seq){
    
    for (disp in dispersal_scale){
      
      for (S in S_seq) { 
        
        for (cinter in c_inter_seq){
          
          for (branch in branches){
            
            if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
              state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
            }else {
              state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
            }
            
            d2 = tibble()
            
            
            for (psi2 in psi_seq){
              
              for (psi1 in psi1_seq[which(psi1_seq>psi2)]){
                
                julia_assign("state", state)
                param=Get_PA_parameters()
                param["delta"]=disp
                param["cintra"]=.3
                param["alpha_0"]=cinter
                param["S"] = S
                param["psi_1"]=psi1
                param["psi_2"]=psi2
                param["f"]=facil
                julia_assign("p", param)
                
                if (scale_f=="local"){
                  prob = julia_eval("ODEProblem(PA_two_species_varying_trait, state, tspan, p)")
                } else {
                  prob = julia_eval("ODEProblem(PA_two_species_varying_trait_global_facilitation, state, tspan, p)")
                }
                
                
                sol = de$solve(prob, de$Tsit5(), saveat = t)
                d = as.data.frame(t(sapply(sol$u, identity)))
                colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
                
                d2 = rbind(d2, d[nrow(d), ] %>% add_column(Stress = S, Psi2 = psi2, Psi1 = psi1,alpha_0=cinter,
                                                           Branch=branch))
                
              }
            } # end trait value 2nd species
            
            d2[d2 < 10^-4] = 0
            d2$rho_plus = d2$rho_1 + d2$rho_2
            d2=d2[,c(1,2,10:15)]
            colnames(d2) = c("Stress_tolerant", "Competitive", "Stress", "Psi2","Psi1","alpha_0","Branches","Rho_plus")
            write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                                  cinter,"_branch_",branch,
                                  "_stress_",round(S,4),"_delta_",disp,"_facilitation_",facil,"_scalefacilitation_",scale_f,".csv"),sep=";")
            
          } # end loop interspecific competition
          
        } # end loop branch
        
      } # end loop first species trait
      
    } #end loop dispersal
    
  }#end facilitation loop

}#end facilitation scale



## 3) Clustering between species fixed traits  ----


tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)

N_rep = 12
S_seq = c(0,.1,.73,.77)
alpha_seq = c(.2)
f_seq=.9
delta_seq=seq(0,1,length.out=N_rep)
cintra_seq=c(.3)


name_scena=c("global_C_global_F","global_C_local_F")

d_clustering=tibble() #initializing the tibble

for (scena_ID in 1:2){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (aii in cintra_seq){ #varying intraspecific competition strength
      
      for (f in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) { #varying competition
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          param["cintra"]=aii
          param["f"]=f
          param["delta"]=disp
          param["alpha_0"]=alpha0
          
          
          state=Get_PA_initial_state()
          julia_assign("state", state)
          
          #varying the global interspecific competition
          
          d2 = tibble()
          
          for (S in S_seq) { #varying the stress 
            
            param["S"] = S
            param["alpha_0"] = alpha0
            julia_assign("p", param)
            
            if (scena_ID==1){  #global C, global F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }else if (scena_ID==2){ #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
              
            }
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, alpha_0=alpha0))
            
          }
          d2[d2 < 10^-4] = 0
          colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
          d2$rho_plus = d2$rho_1 + d2$rho_2
          
          
          
          d_clustering=rbind(d_clustering,tibble(
            Rho_1=d2$rho_1,Rho_2=d2$rho_2,Rho_12=d2$rho_12,
            Rho_22=d2$rho_22,Rho_11=d2$rho_11,Rho_10=d2$rho_1 - d2$rho_11 - d2$rho_12 - d2$rho_1m,
            Rho_20=d2$rho_2 - d2$rho_22 - d2$rho_12 - d2$rho_2m,
            Rho_0=1-d2$rho_1-d2$rho_1-d2$rho_m,
            S   = d2$S,alpha_0 = d2$alpha_0,
            f=f,delta=disp,Scena=scena_ID,
            cintra=aii
          ))
          
          
        } #end competition loop
        
      } #end facilitation loop
      
    } #end h loop
    
  } #end dispersal loop
  
} #end scenario loop

write.table(d_clustering,"../Table/2_species/PA/Clustering_PA.csv",sep=";")





## 4) Threshold for invasion and extinction ----

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)
N_rep = 100
S_seq = c(seq(0,.75,length.out=100),seq(.75,.9,length.out=1400))
alpha_seq = c(.2)
f_seq=.9
delta_seq=seq(0,1,length.out=12)
aii=.3


name_scena=c("global_C_global_F","global_C_local_F")

d_invade=d_extinction=tibble() #initializing the tibble

for (scena_ID in 1:2){ #for each scenario of species pairs
  
  for (disp in delta_seq){ #varying dispersal scale
    
    for (traj in c("Restoration","Degradation")){ #varying intraspecific competition strength
      
      for (f in f_seq){ #varying facilitation strength
        
        for (alpha0 in alpha_seq) { #varying competition
          
          
          #Setting the parameters
          param=Get_PA_parameters()
          param["cintra"]=aii
          param["f"]=f
          param["delta"]=disp
          param["alpha_0"]=alpha0
          
          if (traj=="Restoration"){
            state=Get_PA_initial_state(ini =c(.05,.05,.49))
          } else{
            state=Get_PA_initial_state(ini =c(.4,.4,.1))
          }
          
          julia_assign("state", state)
          
          #varying the global interspecific competition
          
          d2 = tibble()
          
          for (S in S_seq) { #varying the stress 
            
            param["S"] = S
            param["alpha_0"] = alpha0
            julia_assign("p", param)
            
            if (scena_ID==1){  #global C, global F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_global_F, state, tspan, p)")
              
            }else if (scena_ID==2){ #global C, local F
              
              julia_assign("p", param)
              prob = julia_eval("ODEProblem(PA_two_species_global_C_local_F, state, tspan, p)")
              
            }
            
            sol = de$solve(prob, de$Tsit5(), saveat = t)
            d = as.data.frame(t(sapply(sol$u, identity)))
            
            colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
            
            d2 = rbind(d2, d[nrow(d),] %>% add_column(S = S, alpha_0=alpha0))
            
          }
          d2[d2 < 10^-4] = 0
          colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha_0")
          d2$rho_plus = d2$rho_1 + d2$rho_2
          
          
          
          if (traj=="Restoration"){
            
            d_invade=rbind(d_invade,tibble(
              Thresh_invasion=max(d2$S[which((d2$rho_1)>0)]),
              alpha_0 = alpha0,
              f=f,delta=disp,Scena=scena_ID,
              cintra=aii
            ))
            
          } else {
            
            d_extinction=rbind(d_extinction,tibble(
              Thresh_extinction=min(d2$S[which((d2$rho_1)==0)]),
              alpha_0 = alpha0,
              f=f,delta=disp,Scena=scena_ID,
              cintra=aii
            ))
            
            
          }         
          
        } #end competition loop
        
      } #end facilitation loop
      
    } #end trajectory loop
    
  } #end dispersal loop
  
} #end scenario loop
write.table(d_invade,"../Table/2_species/PA/Threshold_invasion.csv",sep=";")
write.table(d_extinction,"../Table/2_species/PA/Threshold_extinction.csv",sep=";")




## 5) Varying the trade-off shape ----

rm(list = ls())
source("./Dryland_shift_functions.R")
julia_setup()
de = diffeq_setup()

tspan = c(0, 2000) #to avoid long transient
t = seq(0, 2000, by = 1)
julia_library("DifferentialEquations")
julia_assign("tspan", tspan)


N_sim=50
S_seq =seq(0,.82, length.out = 100)
psi_seq=seq(0,1,length.out=N_sim)
c_inter_seq=.3
psi1_seq=seq(0,1,length.out=N_sim)
f_seq=c(.9)
trade_off_shape=c(.5,1.5)
branches=c("Degradation","Restoration")

for (facil in f_seq){
  
  for (tradeoff in trade_off_shape){
    
    for (S in S_seq) { 
      
      
      for (cinter in c_inter_seq){
        
        
        
        for (branch in branches){
          
          if (branch =="Degradation"){ #doing the two branches of the bifurcation diagram
            state = Get_PA_initial_state(Get_MF_initial_state(c(.4,.4,.1)))
          }else {
            state = Get_PA_initial_state(Get_MF_initial_state(c(.005,.005,.49)))
          }
          
          d2 = tibble()
          
          
          for (psi2 in psi_seq){
            
            for (psi1 in psi1_seq[which(psi1_seq>psi2)]){
              
              julia_assign("state", state)
              param=Get_PA_parameters()
              param["delta"]=.1
              param["cintra"]=.3
              param["alpha_0"]=cinter
              param["S"] = S
              param["psi_1"]=psi1
              param["psi_2"]=psi2
              param["f"]=facil
              param["shape"]=tradeoff
              julia_assign("p", param)
              
              prob = julia_eval("ODEProblem(PA_two_species_varying_trait_trade_off, state, tspan, p)")
              
              
              sol = de$solve(prob, de$Tsit5(), saveat = t)
              d = as.data.frame(t(sapply(sol$u, identity)))
              colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
              d2 = rbind(d2, d[nrow(d), ] %>% add_column(Stress = S, Psi2 = psi2, Psi1 = psi1,alpha_0=cinter,
                                                         Branch=branch,Tradeoff=tradeoff))
              
            }
          } # end trait value 2nd species
          
          d2[d2 < 10^-4] = 0
          d2$rho_plus = d2$rho_1 + d2$rho_2
          d2=d2[,c(1,2,10:16)]
          colnames(d2) = c("Sp1", "Sp1", "Stress", "Psi2","Psi1","alpha_0","Branches","Tradeoff","Rho_plus")
          write.table(d2,paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                                cinter,"_branch_",branch,
                                "_stress_",S,"_tradeoff_",tradeoff,"_facilitation_",facil,".csv"),sep=";")
          
        } # end loop interspecific competition
        
      } # end loop branch
      
    } # end loop first species trait
    
  } #end loop dispersal
  
}#end facilitation loop








# Step 3) N-species analysis ----
rm(list = ls())
source("./Dryland_shift_functions.R")


threshold=.1
d_tipping=d_richness=d_richness2=d_tot=d_position_shift=d_trait_shift=tibble()
# pdf(paste0("../Figures/N_species/MF/Dyn_Nspecies.pdf"),width = 6,height = 4)

for (i in c(5,15,25)){
  
  d_t=tibble()
  
  
  Nsp=i
  list_csv=list.files(paste0('../Table/N_species/MF/',i,'_sp/'))
  
  for ( k in 1:length(list_csv)){
    
    d2=read.table(paste0("../Table/N_species/MF/",i,"_sp/",list_csv[k]),sep=",")
    colnames(d2)=c(paste0("Sp_",1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
    
    d2[d2<10^(-4)]=0
    
    d_t=rbind(d_t,d2)
    
    
    d_richness=rbind(d_richness,
                     tibble(Nsp=Nsp,Random_ini=k,Competition=strsplit(list_csv[k],split = "_")[[1]][5],
                            Branch="Degradation",
                            Richness=length(which(colSums(d2[1:(nrow(d2)),])[1:Nsp] !=0))))
    
    
    # print(
    #   ggplot(d2%>%melt(., measure.vars=paste0("Sp_",1:Nsp)))+
    #     geom_line(aes(x=Stress,y=value,color=variable),size=.8)+
    #     ggtitle(paste0("Nsp = ",i,", aij = ",
    #                    strsplit(list_csv[k],split = "_")[[1]][5]))+
    #     the_theme+labs(y="",color=expression(paste(bar(psi),"    ")))
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
    
    d2$Rho_plus = sapply(1:nrow(d2),function(x){
      return(sum(d2[x,1:Nsp]))
    })
    
    
    
    d_tot=rbind(d_tot,d2[,-c(1:Nsp)]%>%add_column(., Nsp=Nsp,Competition=strsplit(list_csv[k],split = "_")[[1]][5],
                                                Div_species=length(which(colSums(d2)[1:Nsp]>0))    ))
    
    for (sp in 1:Nsp){
      
      
      
      if (any(abs(diff(d2[1:(nrow(d2)),sp]))>threshold)){ #if species shift
        tip=1
        nb_shift=length(which(round((abs(diff(d2[1:(nrow(d2)),sp]))),4)>threshold))
        
        if (nb_shift %% 2==1){ #case for the first species shifting as it is already present in the system
          nb_shift=nb_shift+1
        }
        
        
        d_tipping=rbind(d_tipping,tibble(Species=sp,Nsp=Nsp,Random_ini=k,Tipping=nb_shift/2,Branch="Degradation",
                                         Competition=strsplit(list_csv[k],split = "_")[[1]][5]))
        
        
        d_position_shift=rbind(d_position_shift,tibble(Species=sp,Nsp=Nsp,Random_ini=k,Branch="Degradation",
                                                       Stress_shift=d2$Stress[which(round((abs(diff(d2[1:(nrow(d2)),sp]))),4)>threshold)],
                                                       Competition=strsplit(list_csv[k],split = "_")[[1]][5]))
      }else {
        tip=0
      }
      
      d_trait_shift=rbind(d_trait_shift,tibble(Species=sp,Nsp=Nsp,Random_ini=k,Branch="Degradation",
                                               Tipping=tip,
                                               Competition=strsplit(list_csv[k],split = "_")[[1]][5]))
      
    }
  }
  
  d_t=d_t[order(d_t$Branch),]
  d_richness2=rbind(d_richness2,tibble(Nsp=Nsp,Competition=strsplit(list_csv[k],split = "_")[[1]][5],
                                       Richness_D=length(which(colSums(d_t[1:(nrow(d_t)),])[1:Nsp] !=0))))
}
# dev.off()

write.table(d_tipping,"../Table/N_species/MF/Multistability_tipping.csv",sep=";")
write.table(d_richness,"../Table/N_species/MF/Multistability_richness.csv",sep=";")
write.table(d_richness2,"../Table/N_species/MF/Multistability_richness2.csv",sep=";")
write.table(d_tot,"../Table/N_species/MF/Multistability_CSI.csv",sep=";")
write.table(d_position_shift,"../Table/N_species/MF/Position_shift_stress_gradient.csv",sep=";")
write.table(d_trait_shift,"../Table/N_species/MF/Trait_shift.csv",sep=";")





d_tipping=read.table("../Table/N_species/MF/Multistability_tipping.csv",sep=";")
d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")
d_richness2=read.table("../Table/N_species/MF/Multistability_richness2.csv",sep=";")
d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")
d_position_shift=read.table("../Table/N_species/MF/Position_shift_stress_gradient.csv",sep=";")

d_position_shift$Trait=sapply(1:nrow(d_position_shift),function(x){
  return(seq(1,0, length.out=d_position_shift$Nsp[x])[d_position_shift$Species[x]])
})



# Different metrics for analyzing tipping points in community.
# summarizing these metrics

d_tipping_summarized=d_tipping%>%
  group_by(., Nsp,Random_ini,Branch,Competition)%>%
  summarise(., Nb_different_species=length(unique(Species)),.groups = "keep",
            Nb_tipping=sum(Tipping))


d_tipping_mean=d_tipping_summarized%>%
  group_by(., Nsp,Branch,Competition)%>%
  summarise(., mean_nb_shift=mean(Nb_tipping),.groups = "keep",
            mean_number_shift_sp=mean(Nb_different_species))

d_richness_mean=d_richness%>%
  group_by(., Nsp,Branch,Competition)%>%
  summarise(., mean_richness=mean(Richness),.groups = "keep")


d_position_summarize=d_position_shift%>%
  group_by(.,Random_ini,Branch,Competition,Nsp)%>%
  summarise(.,.groups = "keep",Nb_shift=length(unique(Stress_shift)))

d_position_summarize_mean=d_position_summarize%>%
  group_by(.,Branch,Competition,Nsp)%>%
  summarise(.,.groups = "keep",Nb_shift_mean=mean(Nb_shift))


#Ploting

p1=ggplot(d_tipping_summarized)+
  geom_bar(aes(x=Competition,fill=as.factor(Nb_different_species),color=as.factor(Nb_different_species),
               y=Nb_different_species),
           color="transparent",position = "fill",stat="identity")+
  geom_point(data=d_tipping_mean,aes(x=Competition,y=mean_number_shift_sp/4),
             shape=1,size=3,color="white")+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),y="Fraction of replicates",fill="Unique number of species that shifts")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  scale_y_continuous(sec.axis=sec_axis(~.*4, name="Mean number of unique \n species that shift"))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))


p2=ggplot(d_tipping_summarized)+
  geom_bar(aes(x=Competition,fill=as.factor(Nb_tipping),color=as.factor(Nb_tipping),
               y=Nb_tipping),color="transparent",
           position = "fill",stat="identity")+
  geom_point(data=d_tipping_mean,aes(x=Competition,y=mean_nb_shift/5),
             shape=1,size=3,color="white")+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),y="Fraction of replicates",fill="Number of species shifts")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  scale_y_continuous(sec.axis=sec_axis(~.*5, name="Mean number of species shifts"))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))


p3=ggplot(d_richness)+
  geom_bar(aes(x=Competition,fill=as.factor(Richness),color=as.factor(Richness),y=Richness),
           color="transparent",position = "fill",stat="identity")+
  geom_point(data=d_richness_mean,aes(x=Competition,y=mean_richness/7),
             shape=1,size=3,color="white")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),
       y="Fraction of replicates",fill="Number of species \n with positive cover")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  scale_y_continuous(sec.axis=sec_axis(~.*7, name="Mean total number of species"))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))


p4=ggplot(d_position_summarize)+
  geom_bar(aes(x=Competition,fill=as.factor(Nb_shift),color=as.factor(Nb_shift),y=Nb_shift),
           color="transparent",position = "fill",stat="identity")+
  geom_point(data=d_position_summarize_mean,aes(x=Competition,y=Nb_shift_mean/4),
             shape=1,size=3,color="white")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),
       y="Fraction of replicates",fill="Number of community transitions")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  theme(strip.text.x = element_text(size=11))+
  scale_y_continuous(sec.axis=sec_axis(~.*4, name="Mean number of community transitions"))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))


ggsave("../Figures/N_species/MF/Nb_shift_diversity.pdf",
       ggarrange(
         ggarrange(p2,p4,ncol = 2,labels=LETTERS[1:2]),
         ggarrange(p1,p3,ncol=2,labels=LETTERS[3:4]),
         nrow=2,align = "hv"),width=19,height=10)



# Community index

d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

p=ggplot(d_tot)+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())

ggsave("../Figures/N_species/MF/CSI_aij_Nsp.pdf",p,width = 12,height = 7)




# Which species shifts and how frequently ?

d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")
N_replicate=150

d_richness$Trait=sapply(1:nrow(d_richness),function(x){
  return(seq(1,0, length.out=d_richness$Nsp[x])[d_richness$Species[x]])
})

p=ggplot(filter(d_richness), 
         aes(x=Trait,y=Richness/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .1)+
  facet_grid(Nsp~Competition)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Frequency of abrupt shift",fill="")+
  scale_fill_gradientn(colours = color_Nsp(100))

ggsave("../Figures/N_species/MF/Traits_which_shifts.pdf",p,width = 10,height = 6)





# Disitribution of abiotic stress conditions at which species shifts

d_position_shift=read.table("../Table/N_species/MF/Position_shift_stress_gradient.csv",sep=";")

d_position_shift$Trait=sapply(1:nrow(d_position_shift),function(x){
  return(seq(1,0, length.out=d_position_shift$Nsp[x])[d_position_shift$Species[x]])
})

p=ggplot(d_position_shift)+
  geom_violin(aes(x=Competition,y=Stress_shift,group=Competition))+
  geom_jitter(aes(x=Competition,y=Stress_shift,color=Trait,fill=Trait),size=.5,shape=21,alpha=.5,height = 0,width = .01)+
  the_theme+labs(y="Stress at which species shifts",color=TeX(r'(Trait of species, \ $\psi)'),
                 x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'))+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  facet_grid(Nsp~.,labeller = label_bquote(rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))

ggsave("../Figures/N_species/MF/Stress_at_which_species_shift.pdf",p,width = 8,height = 8)


#same but with community scale shifts

d_position_shift=read.table("../Table/N_species/MF/Position_shift_stress_gradient.csv",sep=";")

d_position_shift$Trait=sapply(1:nrow(d_position_shift),function(x){
  return(seq(1,0, length.out=d_position_shift$Nsp[x])[d_position_shift$Species[x]])
})


d_position_summarize=d_position_shift%>%
  group_by(.,Random_ini,Branch,Competition,Nsp)%>%
  summarise(.,.groups = "keep",Community_shift=unique(Stress_shift))

p=ggplot(d_position_summarize)+
  geom_violin(aes(x=Competition,y=Community_shift,group=Competition))+
  geom_jitter(aes(x=Competition,y=Community_shift),size=.5,shape=21,alpha=.5,height = 0,width = .01,color="black")+
  the_theme+labs(y="Stress at which community transition",
                 x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'))+
  guides(fill="none")+
  facet_grid(Nsp~.,labeller = label_bquote(rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.x = element_blank())+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025))



ggsave("../Figures/N_species/MF/Stress_at_which_community_shifts.pdf",p,width = 8,height = 8)


