rm(list = ls())
source("./2_Species_analysis_functions.R")
# julia_setup()
# de = diffeq_setup()
#  #

c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)
dispersal_scale=c(.1,.9)
f_seq=c(0,.45,.9)

for (f in f_seq){
  
  for (disp in dispersal_scale){
    
    for (Psi_sp1 in psi1_seq){
      d=tibble()  
      
      for (branch in c("Restoration","Degradation")){
        for (cinter in c_inter_seq){
          
          d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                               cinter,"_branch_",branch,
                               "_Psi1_",Psi_sp1,"_delta_",disp,"_facilitation_",f,".csv"),sep=";")
          d=rbind(d,d2)
          
        } # end loop interspecific competition
      } # end loop branch
      
      
      
      
      d[,1:2][d[,1:2] < 10^-4] = 0
      
      
      #COmputing CSI index
      set.seed(123)
      u=runif(2)
      d$CSI = sapply(1:nrow(d),function(x){
        return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
      })
      
      
      
      
      if (Psi_sp1 ==1) {
        
        d2=filter(d,Psi1==Psi_sp1)
        
        d2$state = sapply(1:nrow(d2), function(x) {
          if (d2[x, 1] > 0 & d2[x, 2] > 0) {
            return("Coexistence")
          }
          if (d2[x, 1] > 0 & d2[x, 2] == 0) {
            return("Stress_tolerant")
          }
          if (d2[x, 1] == 0 & d2[x, 2] > 0) {
            return("Species 2")
          }
          if (d2[x, 1] == 0 & d2[x, 2] == 0) {
            return("Desert")
          }
        })
        
        d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1),]
        
        all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
          if (d2$state[x] != d2$state[x+1]){
            return(paste0(d2$state[x],"/", d2$state[x+1]))
          }
          else {return(d2$state[x])}
        })
        
        d_state=d2%>%
          filter(., Branches=="Degradation")%>%
          select(.,-Branches)
        d_state$all_state=all_state
        
        
        color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
        
        appender <- function(string) {
          TeX(paste("$\\alpha_e = $", string))}
        
        
        
        
        p=ggplot(d_state%>%
                   mutate(all_state=recode_factor(all_state,
                                                  "Species 2/Coexistence"="Coexistence/Species 2",
                                                  "Desert/Species 2"="Species 2/Desert",
                                                  "Desert/Coexistence"="Coexistence/Desert",
                                                  "Desert/Stress_tolerant"="Stress_tolerant/Desert",
                                                  "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                                  "Stress_tolerant/Species 2"="Species 2/Stress_tolerant",
                                                  "Species 2/Coexistence"="Coexistence/Species 2"))) +
          geom_tile(aes(x=Stress,y=as.numeric(1-Psi2),fill=all_state))+
          theme_classic() +
          theme(legend.position = "bottom") +
          labs(x = "Stress (S)", y = TeX(r'(Trait difference \ |$\psi_1-\psi_2|)'), fill = "") +
          theme(legend.text = element_text(size = 11))+
          scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                     "Species 2" = "#ACD87B",
                                     "Coexistence/Species 2" = "#DDEFCA",
                                     "Stress_tolerant" = "#7BD8D3",
                                     "Stress_tolerant/Desert" ="#0F8E87",
                                     "Coexistence/Stress_tolerant"="#9BBBB9",
                                     "Coexistence/Desert"="#C19E5E",
                                     "Desert"=  "#696969",
                                     "Species 2/Stress_tolerant" = "#9B68A0"),
                            labels=c("Coexistence","Species 2","Coexistence/Species 2","Stress-tolerant","Stress-tolerant/Desert",
                                     "Coexistence/Stress-tolerant",
                                     "Coexistence/Desert","Desert","Species 2/Stress-tolerant"))+
          facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
          the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Multistability_trait_variation_sp_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 9,height = 4)
        
        
        
        d2t=transform(d2,
                      alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                              "alpha[e] == 0.5")),
                      Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                                  labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
        
        
        
        
        #Bifurcation diagram species
        p=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                   filter(., Psi2 %in% unique(d_state$Psi2)[c(1,75,100)])%>%
                   mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                   filter(., alpha_0!=.1))+
          geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
          labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
          the_theme +
          scale_color_manual(values = color_rho[c(2, 4)],labels=c("Species 2","Stress-tolerant")) +
          scale_linetype_manual(values=c(1,9))+
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
        
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Bifurcation_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        
        
        #Community index
        p=ggplot(d2%>%
                   filter(., Psi2 %in% unique(d2$Psi2)[c(1,75,100)])%>%
                   mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                   filter(., alpha_0!=.1))+
          geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
          labs(x = "Stress (S)", y = "Community index", color = "",linetype="") +
          the_theme +
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/CSI_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        
        
        
        #Computing hysteresis size : species scale
        d_hysteresis=tibble()
        for (a0 in unique(d2$alpha_0)){
          for (psi in unique(d2$Psi2)){
            
            d_a0=filter(d2,alpha_0==a0,Psi2==psi)
            d_a0=d_a0[order(d_a0$Branches,d_a0$Stress),]
            
            d_hysteresis=rbind(d_hysteresis,tibble(Competition=a0,Psi2=psi,
                                                   Hysteresis_com=abs(sum(d_a0$Rho_plus[1:100]-d_a0$Rho_plus[101:200])),
                                                   Hysteresis_STol=abs(sum(d_a0$Stress_tolerant[1:100]-d_a0$Stress_tolerant[101:200])),
                                                   Hysteresis_Comp=abs(sum(d_a0$Competitive[1:100]-d_a0$Competitive[101:200])),
                                                   Hysteresis_com_scaled=abs(sum(d_a0$Rho_plus[1:100]-d_a0$Rho_plus[101:200]))/(sum(d_a0$Rho_plus[1:100]+rev(d_a0$Rho_plus[101:200]))/2),
                                                   Hysteresis_STol_scaled=abs(sum(d_a0$Stress_tolerant[1:100]-d_a0$Stress_tolerant[101:200]))/(sum(d_a0$Stress_tolerant[1:100]+rev(d_a0$Stress_tolerant[101:200]))/2),
                                                   Hysteresis_Comp_scaled=abs(sum(d_a0$Competitive[1:100]-d_a0$Competitive[101:200]))/(sum(d_a0$Competitive[1:100]+rev(d_a0$Competitive[101:200]))/2)
            ))
            
            
          }
        }
        
        
        
        
        pdf(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Hysteresis_size_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_STol","Hysteresis_Comp")))+
                geom_point(aes(x = Psi2, y = value,color=variable),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition)))+
                scale_color_manual(values=as.character(color_rho[c(4,2)])))
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_STol_scaled","Hysteresis_Comp_scaled")))+
                geom_point(aes(x = Psi2, y = value,color=variable),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition)))+
                scale_color_manual(values=as.character(color_rho[c(4,2)])))
        
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_com")))+
                geom_point(aes(x = Psi2, y = value),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition))))
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_com_scaled")))+
                geom_point(aes(x = Psi2, y = value),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition))))
        
        
        dev.off()
        
        
        
        
        
      } 
      if (Psi_sp1==0){
        
        
        d2=filter(d,Psi1==Psi_sp1)
        
        d2$state = sapply(1:nrow(d2), function(x) {
          if (d2[x, 1] > 0 & d2[x, 2] > 0) {
            return("Coexistence")
          }
          if (d2[x, 1] > 0 & d2[x, 2] == 0) {
            return("Competitive")
          }
          if (d2[x, 1] == 0 & d2[x, 2] > 0) {
            return("Species 2")
          }
          if (d2[x, 1] == 0 & d2[x, 2] == 0) {
            return("Desert")
          }
        })
        
        d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1),]
        
        all_state =sapply(seq(1, nrow(d2) , by = 2),function(x){
          if (d2$state[x] != d2$state[x+1]){
            return(paste0(d2$state[x],"/", d2$state[x+1]))
          }
          else {return(d2$state[x])}
        })
        
        d_state=d2%>%
          filter(., Branches=="Degradation")%>%
          select(.,-Branches)
        d_state$all_state=all_state
        
        
        color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")
        
        appender <- function(string) {
          TeX(paste("$\\alpha_e = $", string))}
        
        
        
        p=ggplot(d_state%>%
                   mutate(all_state=recode_factor(all_state,
                                                  "Species 2/Coexistence"="Coexistence/Species 2",
                                                  "Desert/Species 2"="Species 2/Desert",
                                                  "Desert/Coexistence"="Coexistence/Desert",
                                                  "Desert/Competitive"="Competitive/Desert",
                                                  "Competitive/Coexistence" = "Coexistence/Competitive"
                   ),
                   Psi2=round(Psi2,5),
                   Stress=round(Stress,5))) +
          geom_tile(aes(x=Stress,y=abs(as.numeric(Psi2)),fill=all_state))+
          theme_classic() +
          theme(legend.position = "bottom") +
          labs(x = "Stress (S)", y = TeX(r'(Trait difference \ |$\psi_1-\psi_2|)'), fill = "") +
          theme(legend.text = element_text(size = 11))+
          scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                     "Species 2" = "#7BD8D3",
                                     "Coexistence/Species 2" = "#9BBBB9",
                                     "Species 2/Desert" ="#0F8E87",
                                     "Coexistence/Desert"="#C19E5E",
                                     "Desert"=  "#696969"))+
          facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
          the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Multistability_trait_variation_sp_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 9,height = 4)
        
        
        
        #Bifurcation diagrams
        
        d2t=transform(d2,
                      alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                              "alpha[e] == 0.5")),
                      Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                                  labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
        
        
        p=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                   filter(., Psi2 %in% unique(d_state$Psi2)[c(1,75,100)])%>%
                   mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                   filter(., alpha_0!=.1))+
          geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
          labs(x = "Stress (S)", y = "Density", color = "",linetype="") +
          the_theme +
          scale_color_manual(values = (as.character(color_rho[c(2, 4)])),labels=c("Competitive","Species 2")) +
          scale_linetype_manual(values=c(1,9))+
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
        
        
        
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Bifurcation_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        
        #COmmunity index
        p=ggplot(d2%>%
                   filter(., Psi2 %in% unique(d2$Psi2)[c(1,75,100)])%>%
                   mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                   filter(., alpha_0!=.1))+
          geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
          labs(x = "Stress (S)", y = "Community index", color = "",linetype="") +
          the_theme +
          theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
          facet_grid(Psi2~alpha_0, scales = "free",  
                     labeller = label_bquote(cols = alpha[e] == .(alpha_0), rows = abs(psi[1] - psi[2]) == .(Psi2)))
        
        ggsave(paste0("../Figures/2_species/PA/Multistability/Varying_traits/CSI_varying_traits_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        
        
        
        
        #Computing hysteresis size
        d_hysteresis=tibble()
        for (a0 in unique(d2$alpha_0)){
          for (psi in unique(d2$Psi2)){
            
            d_a0=filter(d2,alpha_0==a0,Psi2==psi)
            d_a0=d_a0[order(d_a0$Branches,d_a0$Stress),]
            
            d_hysteresis=rbind(d_hysteresis,tibble(Competition=a0,Psi2=psi,
                                                   Hysteresis_com=abs(sum(d_a0$Rho_plus[1:100]-d_a0$Rho_plus[101:200])),
                                                   Hysteresis_Comp=abs(sum(d_a0$Stress_tolerant[1:100]-d_a0$Stress_tolerant[101:200])),
                                                   Hysteresis_STol=abs(sum(d_a0$Competitive[1:100]-d_a0$Competitive[101:200])),
                                                   Hysteresis_com_scaled=abs(sum(d_a0$Rho_plus[1:100]-d_a0$Rho_plus[101:200]))/(sum(d_a0$Rho_plus[1:100]+rev(d_a0$Rho_plus[101:200]))/2),
                                                   Hysteresis_Comp_scaled=abs(sum(d_a0$Stress_tolerant[1:100]-d_a0$Stress_tolerant[101:200]))/(sum(d_a0$Stress_tolerant[1:100]+rev(d_a0$Stress_tolerant[101:200]))/2),
                                                   Hysteresis_STol_scaled=abs(sum(d_a0$Competitive[1:100]-d_a0$Competitive[101:200]))/(sum(d_a0$Competitive[1:100]+rev(d_a0$Competitive[101:200]))/2)
            ))
            
            
          }
        }
        
        
        pdf(paste0("../Figures/2_species/PA/Multistability/Varying_traits/Hysteresis_size_Psi1_",Psi_sp1,"_dispersal_",disp,"_facilitation_",f,".pdf"),p,width = 7,height = 4)
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_STol","Hysteresis_Comp")))+
                geom_point(aes(x = Psi2, y = value,color=variable),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition)))+
                scale_color_manual(values=as.character(color_rho[c(4,2)])))
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_STol_scaled","Hysteresis_Comp_scaled")))+
                geom_point(aes(x = Psi2, y = value,color=variable),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition)))+
                scale_color_manual(values=as.character(color_rho[c(4,2)])))
        
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_com")))+
                geom_point(aes(x = Psi2, y = value),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition))))
        
        print(ggplot(d_hysteresis%>%melt(., measure.vars=c("Hysteresis_com_scaled")))+
                geom_point(aes(x = Psi2, y = value),size=.5) +
                labs(x = TeX("$\\psi_2$"), y = "Hysteresis size", color = "",linetype="") +
                the_theme +
                theme(legend.text = element_text(size = 9),strip.text.y = element_text(size=9))+
                facet_grid(.~Competition, scales = "free",  
                           labeller = label_bquote(cols = alpha[e] == .(Competition))))
        
        
        dev.off()
        
        
        
      }
    }
  }
}
