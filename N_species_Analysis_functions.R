rm(list = ls())
x <- c(
  "tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2",
  "JuliaCall", "diffeqr", "simecol", "tseries", "phaseR","GGally",
  "ggquiver", "scales", "boot", "spatialwarnings","simecol","igraph","hillR"
)
lapply(x, require, character.only = TRUE)


the_theme <- theme_classic() + theme(
  legend.position = "bottom",
  strip.background = element_rect(fill = "#CCE8D8"),
  strip.text.y = element_text(size = 10, angle = -90),
  strip.text.x = element_text(size = 8), axis.text = element_text(size = 11), axis.title = element_text(size = 13),
  legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)

color_rho = c("coexistence" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")
color_Nsp=colorRampPalette(c("#077D10",as.character(color_rho[2]),as.character(color_rho[4]),"#2A39EF"))
set.seed(123)
Vec_densities_ass=runif(15)



dir.create("../Table/N_species", showWarnings = F)
dir.create("../Table/N_species/MF", showWarnings = F)
dir.create("../Figures/N_species", showWarnings = F)
dir.create("../Figures/N_species/MF", showWarnings = F)


# Species niche : area and range
Get_species_niche = function(d,Nsp){
  
  d_melt=melt(d,measure.vars = paste0("Sp_",1:Nsp))
  d_melt$value[d_melt$value<10^(-3)]=0
  delta_S=unique(round(abs(diff(d$Stress)),5))[1]
  
  
  d_niche=tibble()
  for (ini in unique(d$Random_ini)){
    for (f in unique(d$Facilitation)){
      for (disp in unique(d$Dispersal)){
        for (c in unique(d$Competition)){ #for each scenario
          for (sp in unique(d_melt$variable)){ #for each species
            
            d_sp=filter(d_melt,variable==sp,Competition==c,Facilitation==f,Random_ini==ini,Dispersal==disp)
            
            d_D=filter(d_sp,Branch=="Degradation")
            
            # Tipping_point_D = d_D$Stress[min(which(diff(as.vector((d_D$value)))==min(diff(as.vector((d_D$value))))))]
            
            d_R=filter(d_sp,Branch=="Restoration")
            # Tipping_point_R = d_R$Stress[min(which(diff(as.vector((d_R$value)))==min(diff(as.vector((d_R$value))))))]
            
            d_niche=rbind(d_niche,tibble(Niche_range_D = abs(diff(range(d_D$Stress[which(d_D$value!=0)]))), #range of niche degradation
                                         Niche_range_R = abs(diff(range(d_R$Stress[which(d_R$value!=0)]))), #range of niche restoration
                                         Delta_niche = abs(diff(range(d_D$Stress[which(d_D$value!=0)])))-abs(diff(range(d_R$Stress[which(d_R$value!=0)]))), #between degradation & restoration trajectories
                                         Species=sp,
                                         Dispersal=disp,
                                         Competition=c,
                                         Random_ini=ini,
                                         Facilitation=f))
            
            #length(which(d_R$value>0))*delta_S
            
            
            
            
            
            
            
          }
        }
      }
    }
  }
  return(d_niche)
}

# Number of abrupt shifts
Get_number_shifts = function(d,tresh,traits){
  
  d_melt=melt(d,measure.vars = paste0("Sp_",1:Nsp))%>%
    mutate(., Trait=rep(traits,each=nrow(d)))
  d_melt$value[d_melt$value<10^{-3}]=0
  
  d_shift=tibble()
  for (s in unique(d$Scena)){
    for (c in unique(d$Rela_c)){ #for each scenario
      for (sp in unique(d_melt$variable)){ #for each species
        d_sp=filter(d_melt,variable==sp,Rela_c==c,Scena==s)
        # plot(abs(diff(d_sp$value)),ylab=sp)
        # plot(abs((d_sp$value)),ylab=sp)
        if (length(which(abs(diff(d_sp$value))>tresh))>0){ #i.e. there is a shift
          shift=1
        } else {shift=0} #no shift
        d_shift=rbind(d_shift,tibble(Species=sp,
                                     Rela_c=c,
                                     Scena=s,
                                     Trait=unique(d_sp$Trait),
                                     Shift=shift))
      }
    }
  }
  return(d_shift)
}

# Hysteresis per species
Compute_hysteresis = function(d, Nsp = 15,tresh=0.01){
  
  d_hysteresis = tibble()
  
  for (scena in unique(d$Scena)) { 
    for (rela_c in unique(d$Rela_c)){
      for (com_compo in unique(d$Community_compo)){
        
        d2=filter(d,Rela_c==rela_c,Scena==scena,Community_compo==com_compo)

        for (sp in 1:Nsp) {
          d2_sp = d2[, c(paste0("Sp_", sp), "Stress", "Branch")]
          d2_sp[,1][d2_sp[,1]<10^{-4}]=0
          
          
          biomass_D = filter(d2_sp, Branch == "Degradation")
          biomass_R = filter(d2_sp, Branch == "Restoration")

          hysteresis_area = sum(biomass_D[, 1] - rev(biomass_R[, 1]))
          
          if (abs(min(diff(as.vector(t(biomass_D[,1])))))>tresh){
          
            hysteresis_range = biomass_D$Stress[min(which(diff(as.vector(t(biomass_D[,1])))==min(diff(as.vector(t(biomass_D[,1]))))))]  -
                               biomass_R$Stress[min(which(diff(as.vector(t(biomass_R[,1])))==min(diff(as.vector(t(biomass_R[,1]))))))]
          
          }else {hysteresis_range=0}
          
          d_hysteresis = rbind(d_hysteresis, tibble(Hysteresis_area = hysteresis_area, 
                                                    Hysteresis_range = hysteresis_range, 
                                                    Species = sp,
                                                    Community_compo=com_compo,Rela_c=rela_c,Scena=scena))
          
        }
      }
    }
  }
  return(d_hysteresis)
}

# Plot landscape
Plot_landscape = function(landscape,Nsp=15){
  color_CA = color_Nsp(Nsp)
  
  ggplot(melt(landscape)) +
    geom_tile(aes(x = Var1, y = Var2, fill = value)) +
    theme_transparent() +
    scale_fill_gradientn(colours = color_CA) +
    theme(panel.border = element_blank()) +
    theme(legend.position = "bottom") +
    labs(fill = "")
}

# Create co-occurrence matrix from z-scores and plot the co-occurrence network

Co_occurrence_matrix = function(df,traits){
  color_CA = color_Nsp(length(traits))
  color_edge=c("red","blue")
    
  d=df %>%
    pivot_wider(names_from = c("Var1"), values_from = value)%>%
    select(., -Var2)%>%
    replace(is.na(.), 0) #Get adjacency matrix
    
  d=d[-1,-1] #Removing fertile patches
  
  d[d > -1.96 & d < 1.96]=0;d[d < -1.96]=-1;d[d > 1.96] = 1  
  
  #  colnames(d)=round(as.numeric(traits[as.numeric(gsub("Sp_",replacement = "",x=colnames(d)))]),4)
  
  #getting colors for igraph plot = 
  inter=as.numeric(as.matrix(d))
  inter=inter[inter!=0] 
  inter[inter==1] = 2
  inter[inter==-1] = 1
  
  #vertex colors
  d_col=data.frame(Sp=1:length(traits),Traits=traits)
  d_col=d_col[order(d_col$Traits),]
  color_CA=color_CA[d_col$Sp]
  
  graph_d=graph_from_adjacency_matrix(as.matrix(d),weighted = T)
  
  V(graph_d)$color=color_CA
  V(graph_d)$label=""
  E(graph_d)$color=color_edge[inter]
  E(graph_d)$arrow.size=0
  E(graph_d)$width=2
  
  
  return(graph_d)    
}


Get_diversity_community=function(trait,densities){
  
  densities=as.numeric(densities) #to avoid errors
  trait=as.numeric(trait)
  
  
  trait=matrix(trait,length(trait),1)
  colnames(trait)="trait1";rownames(trait)=paste0("Sp",1:length(trait))
  
  densities = matrix(densities,1,length(densities))
  colnames(densities)=paste0("Sp",1:length(trait));rownames(densities)="comm1"
  
  
  diversity=d_name=c()
  for (type in c("FD","D")){ #we do that for functional diversity and species diversity
    for (q in 0:2){
      if (type=="D"){ div = hill_taxa(densities,q=q)}
      if (type=="FD"){ div = hill_func(densities,trait,q=q)["FD_q",1]}
      diversity=c(diversity,div)
      d_name= c(d_name,paste0(type,"_",q))

    }
  }
  
  #adding community-weighted mean of traits
  d_diversity=as_tibble(matrix(diversity,1,6))
  colnames(d_diversity)=d_name
  return(d_diversity)
}

