source('./N_species_Analysis_functions.R')
source('./2_species_Analysis_functions.R')
library(grid)



# Multistability but with mean-field model ----

d2=read.table(paste0("../Table/2_species/MF/Multistability_MF.csv"),sep=";")

d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Stress_tolerant")
  }
  if (d2[x, 1] == 0 & d2[x, 2] > 0) {
    return("Competitive")
  }
  if (d2[x, 1] == 0 & d2[x, 2] == 0) {
    return("Desert")
  }
})

d2=d2[order(d2$alpha_0,d2$Stress),]

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

p=ggplot(d_state) +
  geom_tile(aes(x=Stress,y=alpha_0,fill=all_state))+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Stress_tolerant" = "#7BD8D3",
                             "Stress_tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress_tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ),labels=c("Stress_tolerant"="Stress-tolerant","Stress_tolerant/Desert"="Stress-tolerant/Desert","Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant"))+
  geom_hline(yintercept = round(unique(d2$alpha_0)[c(1,100)],4))+
  geom_text(data=tibble(x=c(1.025,1.025),y=c(0.01,.31),lab=c("b",'c')),aes(x=x,y=y,label=lab))+
  the_theme+theme(legend.text = element_text(size=10))


for (u in 1:2){
  d_plot=filter(d2,round(alpha_0,4) == round(unique(d2$alpha_0)[c(1,100)[u]],4))
  assign(paste0("p_1_",u),
         ggplot(d_plot%>%melt(., measure.vars=c("Stress_tolerant","Competitive")))+
           geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
           labs(x = "Stress (S)", y = "Cover", color = "",linetype="") +
           the_theme +
           scale_color_manual(values = color_rho[c(2, 4)]) +
           scale_linetype_manual(values=c(1,9))+
           theme(legend.text = element_text(size = 11))+
           scale_x_continuous(breaks = c(0,.5,1))+guides(color=F)
  )
  
  assign(paste0("p_2_",u),
         ggplot(d_plot%>%melt(., measure.vars=c("Rho_plus")))+
           geom_line(aes(x = Stress, y = value, linetype=Branches),lwd=.8) +
           labs(x = "Stress (S)", y = "Cover", color = "",linetype="") +
           the_theme +
           scale_color_manual(values = color_rho[c(2, 4)]) +
           scale_linetype_manual(values=c(1,9))+
           theme(legend.text = element_text(size = 11))+
           scale_x_continuous(breaks = c(0,.5,1))+guides(color=F)
  )
  
  
}

p2=ggarrange(p_1_2+xlab(""),p_2_2+xlab(""),p_1_1+xlab(""),p_2_1,nrow = 4,common.legend = T,legend = "bottom",labels = c(letters[2],"",letters[3],""))


ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_MF.pdf",ggarrange(p,p2,ncol = 2,widths = c(2,1),labels=c(letters[1],"")),
       width = 10,height = 6)




















# Multistability varying traits with mean-field model ----

c_inter_seq=c(0,.1, .2, .3)
psi1_seq=c(1,0)



for (Psi_sp1 in psi1_seq){
  d=tibble()  
  
  for (branch in c("Restoration","Degradation")){
    for (cinter in c_inter_seq){
      
      d2=read.table(paste0("../Table/2_species/MF/Varying_traits/Multistability_varying_trait_interspe_comp_",
                           cinter,"_branch_",branch,
                           "_Psi1_",Psi_sp1,".csv"),sep=";")
      
      d=rbind(d,d2)
      
    } # end loop interspecific competition
  } # end loop branch
  
  
  
  
  d2=d
  d2[,1:2][d2[,1:2] < 10^-4] = 0
  
  
  #COmputing CSI index
  set.seed(123)
  u=runif(2)
  d2$CSI = sapply(1:nrow(d2),function(x){
    return(u[1]*d2$Stress_tolerant[x]+u[2]*d2$Competitive[x])
  })
  
  
  
  if (Psi_sp1 ==1) {
    
    
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
    
    
    
    
    p1=ggplot(d_state%>%
               mutate(all_state=recode_factor(all_state,
                                              "Species 2/Coexistence"="Coexistence/Species 2",
                                              "Desert/Species 2"="Species 2/Desert",
                                              "Desert/Coexistence"="Coexistence/Desert",
                                              "Desert/Stress_tolerant"="Stress_tolerant/Desert",
                                              "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                              "Stress_tolerant/Species 2"="Species 2/Stress_tolerant",
                                              "Species 2/Coexistence"="Coexistence/Species 2"))) +
      geom_tile(aes(x=Stress,y=as.numeric(Psi2),fill=all_state))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = "Stress (S)", y = TeX(r'(Trait species 2, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Species 2" = "#ACD87B",
                                 "Coexistence/Species 2" = "#DDEFCA",
                                 "Stress_tolerant" = "#7BD8D3",
                                 "Stress_tolerant/Desert" ="#0F8E87",
                                 "Coexistence/Stress_tolerant"="#9BBBB9",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Desert"=  "#696969",
                                 "Species 2/Stress_tolerant" = "#C998CE"),
                        labels=c("Coexistence","Species 2","Coexistence/Species 2","Stress-tolerant","Stress-tolerant/Desert",
                                 "Coexistence/Stress-tolerant",
                                 "Coexistence/Desert","Desert","Species 2/Stress-tolerant"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
      ggtitle(TeX("$\\psi_1 = 1$"))
    
    
    
    
    
  } 
  if (Psi_sp1==0){
    
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
    
    
    
    p2=ggplot(d_state%>%
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
      labs(x = "Stress (S)", y = TeX(r'(Trait species 2, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Species 2" = "#7BD8D3",
                                 "Competitive" = "#ACD87B",
                                 "Coexistence/Species 2" = "#9BBBB9",
                                 "Species 2/Desert" ="#0F8E87",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Coexistence/Competitive" = "#DDEFCA",
                                 "Coexistence/Species 2"="#9BBBB9",
                                 "Desert"=  "#696969"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                             name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
      ggtitle(TeX("$\\psi_1 = 0$"))
    
    
  }
}

ggsave("../Figures/Final_figs/SI/Multistability_varying_traits_MF.pdf",
       ggarrange(p1+theme(axis.line.x = element_blank(),axis.text.x=element_blank(),axis.title.x = element_blank(),axis.ticks.x=element_blank()),
                 p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),
                 nrow = 2,heights=c(1,1),labels=letters[1:2],common.legend = T,legend = "bottom"),
       width = 10,height = 6)



#NintA for MF model and PA model ----

load(file="../Table/2_species/MF/Comparing_net_effects.RData")
d_net=d_all[[1]];d_RNE=d_all[[2]]


d_net$value[d_net$value>10]=NA
d_net$value[d_net$value < -10]=NA

p1=ggplot(d_RNE%>%
            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
            #mutate(., value=rescale(value,to=c(-1,1)))%>%
            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_wrap(.~variable)+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .1))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .32,x=1.05),label="b")+
  geom_hline(data = subset(d_RNE%>%
                             melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                             mutate(., value=rescale(value,to=c(-1,1)))%>%
                             mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                           variable == "Competitive"), aes(yintercept = .3))+
  geom_text(data = subset(d_RNE%>%
                            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
                            mutate(., value=rescale(value,to=c(-1,1)))%>%
                            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")),
                          variable == "Competitive"), aes(y = .12,x=1.05),label="c")+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")

alpha_for_plot=c(.1, .3)

for (alpha in 1:2){
  assign(paste0("p1_",alpha),ggplot(d_RNE%>%
                                     filter(., alpha_0==alpha_for_plot[alpha])%>%
                                     melt(., measure.vars=c("NintA_comp")))+#%>%
           #           mutate(., value=rescale(value,to=c(-1,1))))+
           geom_hline(yintercept = 0,linetype=9,lwd=.5)+
           geom_line(aes(x=S,y=value))+
           the_theme+theme(axis.text = element_text(size =9), axis.title = element_text(size = 10))+
           labs(x="Stress (S)",y=TeX("$\\NInt_{A}$")))
}

p_right=ggarrange(p1_2+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank())+
                    labs(x=""),p1_1,ggplot() + theme_void(),nrow=3,labels=letters[2:3],heights = c(1,1,.1))

p_tot_1=ggarrange(p1+theme(legend.position = "none",
                           axis.title.x = element_blank(),
                           axis.text.x = element_blank(),axis.line.x = element_blank(),
                           axis.ticks.x = element_blank()),p_right,widths = c(2.2,1),labels = c(letters[1],""))



load(file="../Table/2_species/PA/Comparing_net_effects.RData")
d_net=d_all[[1]];d_RNE=d_all[[2]]



d_RNE[,3:6][d_RNE[,3:6] > 2]=NA
d_RNE[,3:6][d_RNE[,3:6] < -1]=NA

p2=ggplot(d_RNE%>%
            melt(., measure.vars=c("NintA_st","NintA_comp"))%>%
            filter(.,Disp==.1,Scale=="LocalF_GlobalC")%>%
            #mutate(., value=rescale(value,to=c(-1,1)))%>%
            mutate(., variable=recode_factor(variable,"NintA_st"='Stress-tolerant',"NintA_comp"="Competitive")))+
  geom_tile(aes(x=S,y=alpha_0,fill=value))+
  facet_grid(.~factor(variable,levels=c("Stress-tolerant","Competitive")))+
  the_theme+labs(x="Stress (S)",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="")+
  scale_fill_gradient2(low="#F73030",mid="white",high="#185BB9")+
  geom_hline(data=tibble(variable="Competitive",yint=c(0.1,.3)),aes(yintercept=yint,variable="Competitive"))+
  geom_text(data=tibble(variable="Competitive",y=c(0.13,.33),x=c(1.04,1.04),label=c("f","e")),aes(y=y,x=x,label=label,variable="Competitive"))


alpha_for_plot=c(unique(d_RNE$alpha_0)[(50/3)+1], .3)

for (alpha in 1:2){
  assign(paste0("p2_",alpha),ggplot(d_RNE%>%
                                      filter(., alpha_0==alpha_for_plot[alpha],Scale=="LocalF_GlobalC",Disp==.1)%>%
                                      melt(., measure.vars=c("NintA_comp")))+#%>%
           #           mutate(., value=rescale(value,to=c(-1,1))))+
           geom_hline(yintercept = 0,linetype=9,lwd=.5)+
           geom_line(aes(x=S,y=value))+
           the_theme+theme(axis.text = element_text(size =9), axis.title = element_text(size = 10))+
           labs(x="Stress (S)",y=TeX("$\\NInt_{A}$")))
}

p_right=ggarrange(p2_2+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                             axis.line.x = element_blank(),
                             axis.ticks.x = element_blank())+
                    labs(x=""),p2_1,ggplot() + theme_void(),nrow=3,labels=letters[5:6])

p_tot_2=ggarrange(p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank()),p_right,widths = c(2.2,1),labels = c(letters[4],""))

p_tot=ggarrange(p_tot_1,p_tot_2,
                nrow=2,heights = c(1,1.4))
ggsave("../Figures/Final_figs/SI/NintA_MF_PA.pdf",p_tot,width = 7,height = 6)




# Niche expansion in MF model ----

d_niche = read.table("../Table/2_species/MF/Niche_expansion_varying_traits.csv",sep=";")
d_niche = do.call(data.frame,lapply(d_niche,function(x) replace(x, is.infinite(x), NA)))


d_niche=read.table("../Table/2_species/MF/Niche_expansion_MF.csv",sep=";")
alpha_seq=seq(0,.3,length.out=10)
d_niche$alpha_0=alpha_seq
p1=ggplot(d_niche%>%
            melt(.,measure.vars=c("Delta_niche_2")))+
  geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
  the_theme+
  scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
  labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="% of competitive species niche change")


#niche varying traits
d_niche = read.table("../Table/2_species/MF/Niche_expansion_varying_traits.csv",sep=";")
d_niche = do.call(data.frame,lapply(d_niche,function(x) replace(x, is.infinite(x), NA)))

d_niche$Delta_niche_2[1:(nrow(d_niche)/2)] = d_niche$Delta_niche_1[1:(nrow(d_niche)/2)] 



p3_1=ggplot(NULL)+
  geom_path(data=d_niche%>%filter(., Psi1==0)%>%
              melt(., measure.vars=c("Delta_niche_2"))%>%
              mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
              mutate(., Delta_psi=abs(Psi1-Psi2),
                     Facilitation=as.factor(Facilitation),
                     alpha_0=as.factor(alpha_0)),
            aes(x=Psi2,y=value,group=interaction(alpha_0,Facilitation),color=alpha_0),lwd=1)+
  geom_point(data=d_niche%>%filter(., Psi1==0)%>%
               melt(., measure.vars=c("Delta_niche_2"))%>%
               mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
               mutate(., Delta_psi=abs(Psi1-Psi2),
                      Facilitation=as.factor(Facilitation),
                      alpha_0=as.factor(alpha_0))%>%
               filter(., Psi2 %in% unique(.$Psi2)[seq(1,length(unique(.$Psi2)),by =5)]), #to only have 10 points
             aes(x=Psi2,y=value,shape=Facilitation,color=alpha_0),fill="white",size=2)+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  labs(y = "% of species 1 niche change", x = TeX(r'(Trait species 2, \ $\psi_2)'),
       color=TeX(r'(\ \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f)'))+ggtitle(TeX("$\\psi_1 = 0$"))+
  the_theme

p3_2=ggplot(NULL)+
  geom_path(data=d_niche%>%filter(., Psi1==1)%>%
              melt(., measure.vars=c("Delta_niche_2"))%>%
              mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
              mutate(., Delta_psi=abs(Psi1-Psi2),
                     Facilitation=as.factor(Facilitation),
                     alpha_0=as.factor(alpha_0)),
            aes(x=Psi2,y=value,group=interaction(alpha_0,Facilitation),color=alpha_0),lwd=1)+
  geom_point(data=d_niche%>%filter(., Psi1==1)%>%
               melt(., measure.vars=c("Delta_niche_2"))%>%
               mutate(., variable=recode_factor(variable,"Delta_niche_1"="Species 1","Delta_niche_2"="Species 2"))%>%
               mutate(., Delta_psi=abs(Psi1-Psi2),
                      Facilitation=as.factor(Facilitation),
                      alpha_0=as.factor(alpha_0))%>%
               filter(., Psi2 %in% unique(.$Psi2)[seq(1,length(unique(.$Psi2)),by =5)]), #to only have 10 points
             aes(x=Psi2,y=value,shape=Facilitation,color=alpha_0),fill="white",size=2)+
  geom_hline(yintercept = 0,linetype=9)+
  scale_shape_manual(values=c(21,22))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  labs(y = "% of species 2 niche change", x = TeX(r'(Trait species 2, \ $\psi_2)'),
       color=TeX(r'(\ \ $\alpha_e)'),
       shape=TeX(r'(Facilitation \ $\f)'))+ggtitle(TeX("$\\psi_1 = 1$"))+
  the_theme
p3=ggarrange(p3_1,p3_2,labels = letters[2:3],common.legend = T,legend = "bottom")

p_tot = ggarrange(ggarrange(ggplot() + theme_void(),p1,ggplot() + theme_void(),widths = c(.7,2,.7),ncol=3,labels=c("",letters[1],"")),
                     p3,nrow=2,labels=c("",""))

ggsave("../Figures/Final_figs/SI/Niche_expansion_MF.pdf",p_tot,width = 8,height = 8)

