rm(list=ls())
x <- c(
  "tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2",
  "JuliaCall", "diffeqr", "simecol", "tseries", "phaseR",
  "ggquiver", "scales", "boot", "spatialwarnings"
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

dir.create("../Table/N_species", showWarnings = F)
dir.create("../Table/N_species/MF", showWarnings = F)
dir.create("../Figures/N_species", showWarnings = F)
dir.create("../Figures/N_species/MF", showWarnings = F)



# Species niche : area and range

Get_species_niche=function(d,Nsp,traits){
  
  d_melt=melt(d,measure.vars = paste0("Sp_",1:Nsp))%>%
    mutate(., Trait=rep(traits,each=nrow(d)))
  d_melt$value[d_melt$value<10^{-3}]=0
  
  d_niche=tibble()
  for (s in unique(d$Scena)){
    for (c in unique(d$Rela_c)){ #for each scenario
      for (sp in unique(d_melt$variable)){ #for each species
        
        d_sp=filter(d_melt,variable==sp,Rela_c==c,Scena==s)
        d_niche=rbind(d_niche,tibble(Niche_range = abs(diff(range(d_sp$Stress[which(d_sp$value!=0)]))), #range of niche
                                     Niche_area = sum(d_sp$value),
                                     Species=sp,
                                     Rela_c=c,
                                     Scena=s,
                                     Trait=unique(d_sp$Trait)))
      }
    }
  }
  return(d_niche)
}

# Number of abrupt shifts

Get_number_shifts=function(d,tresh,traits){
  
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

set.seed(123)
Vec_densities_ass=runif(15)


# Hysteresis per species

Compute_hysteresis = function(d, Nsp = 15,tresh=0.01) {
  
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




