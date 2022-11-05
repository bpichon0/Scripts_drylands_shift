source('./N_species_Analysis_functions.R')
source('./2_species_Analysis_functions.R')
library(grid)




# Figure 2: Scale facilitation, dispersal, multistability ----

d2=read.table(paste0("../Table/2_species/PA/Multistability_fixed_traits_PA.csv"),sep=";")
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

d2=d2[order(d2$alpha_0,d2$Scena,d2$Delta,d2$Stress),]

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


scale="global"
d2t=transform(d_state%>% 
                filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
                mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6)),
              Scena=factor(Scena,levels=c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")),
                           labels=c("Global facilitation","Local facilitation")))

color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")

p1=ggplot(d2t%>%
            mutate(all_state=recode_factor(all_state,
                                           "Desert/Coexistence"="Coexistence/Desert",
                                           "Stress_tolerant/Coexistence"="Coexistence/Stress_tolerant",
                                           "Desert/Stress_tolerant"="Stress_tolerant/Desert"))) +
  geom_tile(aes(x=Stress,y=alpha_0,fill=all_state))+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Competitive/Coexistence" = "#DDEFCA",
                             "Stress_tolerant" = "#7BD8D3",
                             "Stress_tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress_tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(Scena~Delta,labeller=label_bquote(cols = delta == .(Delta)))+
  the_theme+theme(strip.text.x = element_text(size=12))+
  geom_hline(data=tibble(Scena="Local facilitation",Delta=.1,y=c(0,.3)),aes(yintercept=y))+
  geom_text(data= tibble(Scena="Local facilitation",Delta=.1,y=c(0.015,.315),x=1.05,labels=c("b","c")), 
            aes(x=x, y=y,label=labels), alpha=1, colour="black")



for (i in 1:4){
  if (i%%2 ==0){
    
    assign(paste0("p2_",i),
           ggplot(d2%>%melt(., measure.vars=c("Rho_plus"))%>%
                    filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,.3)[abs(i/2)])%>%
                    mutate(.,variable=recode_factor(variable,"Rho_plus"="Global cover")))+
             
             geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
             facet_grid(.~alpha_0,labeller=label_bquote(cols=alpha[e]==.(alpha_0)))+
             the_theme+theme(strip.text.x = element_blank(),strip.background.x = element_blank())+
             scale_color_manual(values="black")+
             labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+guides(color=F)+
             scale_x_continuous(breaks = c(0,.5,1))
           
           
    )
    
  }else{
    
    
    assign(paste0("p2_",i),
           ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
                    filter(., Scena=="global_C_local_F",Delta==.1)%>%filter(., alpha_0 %in% c(0,0,.3)[i])%>%
                    mutate(.,variable=recode_factor(variable,"Stress_tolerant"="Stress-tolerant")))+
             
             geom_path(aes(x=Stress,y=value,linetype=Branches,color=variable,group=interaction(Branches,variable)))+
             facet_grid(.~alpha_0,labeller=label_bquote(cols=alpha[e]==.(alpha_0)))+
             the_theme+theme(strip.text.x = element_blank(),strip.background.x = element_blank())+
             scale_color_manual(values=c(as.character(color_rho[c(4)]),as.character(color_rho[c(2)])))+
             labs(x=TeX(r'(Stress, \ $S)'),y="",linetype="",color="")+guides(color=F)+
             scale_x_continuous(breaks = c(0,.5,1))
           
    )
    
  }
}


p_bottom=ggarrange(p2_1+ylab("Cover"),p2_2,p2_3,p2_4,common.legend = T,legend = "bottom",ncol = 4,align="hv",labels=c(letters[2],"",letters[3],""))

Fig_2=ggarrange(p1,p_bottom,nrow=2,heights = c(3,1),labels = c(letters[1],""))

ggsave("../Figures/Final_figs/Figure_2.pdf",Fig_2,width = 8,height = 9)


# Figure 3: Trait variation and multistability ----

c_inter_seq=c(0, .2, .3)
psi1_seq=c(1,0)
disp=.1
f=0.9

d=tibble()  
for (Psi_sp1 in c(0,1)){
  for (branch in c("Restoration","Degradation")){
    for (cinter in c_inter_seq){
      
      d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_traits/Multistability_varying_trait_interspe_comp_",
                           cinter,"_branch_",branch,
                           "_Psi1_",Psi_sp1,"_delta_",disp,"_facilitation_",f,".csv"),sep=";")
      
      d=rbind(d,d2)
      
    } # end loop interspecific competition
  } # end loop branch
}



d[,1:2][d[,1:2] < 10^-4] = 0


#COmputing CSI index
set.seed(123)
u=runif(2)
d$CSI = sapply(1:nrow(d),function(x){
  return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
})

for (Psi_sp1 in c(0,1)){

  if (Psi_sp1 ==1) {
    
    d2=filter(d,Psi1==1)
    
    
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
                                              "Stress_tolerant"="Species 1",
                                              "Desert/Stress_tolerant"="Species 1/Desert",
                                              "Stress_tolerant/Coexistence"="Coexistence/Species 1",
                                              "Stress_tolerant/Species 2"="Species 2/Species 1",
                                              "Species 2/Coexistence"="Coexistence/Species 2"))) +
      geom_tile(aes(x=Stress,y=as.numeric(Psi2),fill=all_state))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait species 2, \ $\psi_2)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Species 2" = "#ACD87B",
                                 "Coexistence/Species 2" = "#DDEFCA",
                                 "Species 1" = "#7BD8D3",
                                 "Species 1/Desert" ="#0F8E87",
                                 "Coexistence/Species 1"="#9BBBB9",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Desert"=  "#696969",
                                 "Species 2/Species 1" = "#C998CE"),
                        labels=c("Coexistence","Species 2","Coexistence/Species 2","Species 1","Species 1/Desert",
                                 "Coexistence/Species 1",
                                 "Coexistence/Desert","Desert","Species 2/Species 1"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      geom_hline(data=tibble(alpha_0=.2,y=c(.75)),aes(yintercept=y))+
      geom_text(data= tibble(alpha_0=.2,y=c(.8),x=1.05,labels=c("c")), 
                aes(x=x, y=y,label=labels), alpha=1, colour="black")
    
    
    #Bifurcation diagram
    p2_1=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
               filter(., Psi2 %in% unique(d_state$Psi2)[c(75)])%>%
               mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
               filter(., alpha_0==.2))+
      geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Cover", color = "",linetype="") +
      the_theme +
      scale_color_manual(values = color_rho[c(2, 4)],labels=c("Species 2","Stress-tolerant")) +
      scale_linetype_manual(values=c(1,9))+
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))
    
    
    
    #Community index
    p2_2=ggplot(d2%>%
                 filter(., Psi2 %in% unique(d_state$Psi2)[c(75)])%>%
                 mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
                 filter(., alpha_0==.2))+
      geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Community index", color = "",linetype="") +
      the_theme +
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))
    
  
    
    
    
    
    
  } 
  if (Psi_sp1==0){
    
    d2=filter(d,Psi1==0)
    
    
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
                                               "Species 2/Coexistence"="Coexistence/Species 1",
                                               "Desert/Species 2"="Species 1/Desert",
                                               "Desert/Coexistence"="Coexistence/Desert",
                                               "Desert/Competitive"="Species 2/Desert",
                                               "Species 2" = "Species 1",
                                               "Competitive/Coexistence" = "Coexistence/Species 2"
                ),
                Psi2=round(Psi2,5),
                Stress=round(Stress,5))) +
      geom_tile(aes(x=Stress,y=abs(as.numeric(Psi2)),fill=all_state))+
      theme_classic() +
      theme(legend.position = "bottom") +
      labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Trait species 1, \ $\psi_1)'), fill = "") +
      theme(legend.text = element_text(size = 11))+
      scale_fill_manual(values=c("Coexistence" = "#D8CC7B",
                                 "Species 1" = "#7BD8D3",
                                 "Coexistence/Species 1" = "#9BBBB9",
                                 "Species 1/Desert" ="#0F8E87",
                                 "Coexistence/Desert"="#C19E5E",
                                 "Desert"=  "#696969"))+
      facet_grid(.~alpha_0,labeller=label_bquote(cols = alpha[e] == .(alpha_0)))+
      the_theme+theme(legend.text = element_text(size=9),strip.text.x = element_text(size=13))+
      geom_hline(data=tibble(alpha_0=.3,y=c(.5)),aes(yintercept=y))+
      geom_text(data= tibble(alpha_0=.3,y=c(0.55),x=1.05,labels=c("d")), 
                aes(x=x, y=y,label=labels), alpha=1, colour="black")
    
    
    
    #Bifurcation diagrams
    
    d2t=transform(d2,
                  alpha_0 = factor(alpha_0, levels=c(0.1,.2,.3), labels=c("alpha[e] : 0", "alpha[e] : 0.1",
                                                                          "alpha[e] == 0.5")),
                  Psi2=factor(Psi2,levels=unique(d_state$Psi2)[c(1,75,100)],
                              labels=c("|psi[1] - psi[2]|","Stress-tolerant","Competitive")))
    
    
    p2_3=ggplot(d2%>%melt(., measure.vars=c("Stress_tolerant","Competitive"))%>%
               filter(., Psi2 %in% unique(d_state$Psi2)[c(50)])%>%
               mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
               filter(., alpha_0==.3))+
      geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Cover", color = "",linetype="") +
      the_theme +
      scale_color_manual(values = (as.character(color_rho[c(2, 4)])),labels=c("Competitive","Species 2")) +
      scale_linetype_manual(values=c(1,9))+
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))
    
    
    
    
    
    
    #COmmunity index
    p2_4=ggplot(d2%>%
               filter(., Psi2 %in% unique(d2$Psi2)[c(50)])%>%
               mutate(Psi2=round(abs(Psi1-Psi2),2))%>%
               filter(., alpha_0==.3))+
      geom_point(aes(x = Stress, y = CSI),size=.7,shape=21,alpha=.4) +
      labs(x = TeX(r'(Stress, \ $S)'), y = "Community index", color = "",linetype="") +
      the_theme +
      theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))
    
    
  }
}

Figure_3_top=ggarrange(p1+theme(axis.text.x = element_blank(),axis.title.x = element_blank(),
                            axis.line.x = element_blank(),
                            axis.ticks.x = element_blank(),axis.title.y.right = element_text(size = 11))+
            scale_y_continuous(sec.axis = sec_axis(trans = ~ (1-.x) ,
                                                   name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
              ggtitle(TeX(r'($\psi_1 = 1)')),
          p2+theme(strip.background.x = element_blank(),strip.text.x = element_blank(),axis.title.y.right = element_text(size = 11))+
            scale_y_continuous(sec.axis = sec_axis(trans = ~ .x ,
                                                   name = TeX(r'(Trait difference, \ |$\psi_1-\psi_2|)') ))+
            ggtitle(TeX(r'($\psi_2 = 0)')),common.legend = T,legend = "bottom",nrow=2,
          labels = letters[1:2])

Figure_3_bottom=ggarrange(ggarrange(p2_1+theme(legend.position = "none"),p2_2,ncol=2,labels = letters[3]),
                          ggarrange(p2_3+theme(legend.position = "none"),p2_4,ncol=2,labels = letters[4]))

Figure_3_tot=ggarrange(Figure_3_top,Figure_3_bottom,nrow=2,heights = c(3.5,1))

ggsave("../Figures/Final_figs/Figure_3.pdf",Figure_3_tot,width = 9,height = 9)


# Figure 4: Scaling-up to species-rich communities : species versus community shifts----

d_hysteresis=read.table("../Table/N_species/MF/Hysteresis_species.csv",sep=";")



p0=ggplot(NULL)+
  geom_point(data=d_hysteresis,aes(x=Trait,y=Hysteresis_scaled,color=as.factor(Competition)),size=1,width = .1,alpha=.5)+
  the_theme+labs(y="Hysteresis size",color=TeX("$\\alpha_e$"),x=TeX("$\\psi$"),fill=TeX("$\\alpha_e$"))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  ylim(0,.75)+ #to avoid outliers
  facet_wrap(.~Facilitation,labeller = label_bquote(cols=f==.(Facilitation)))+
  scale_x_continuous(breaks = c(0,.25,.5,.75,1))


d_ASS_com=read.table("../Table/N_species/MF/d_ASS_com.csv",sep=";")
d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")

d_sp=d_ASS_sp%>%
  group_by(.,Facilitation,Competition,Branch,Random_ini)%>%
  summarise(.,Tipping=sum(Tipping),.groups = "keep")

d_com=d_ASS_com%>%
  group_by(.,Facilitation,Competition,Branch,Random_ini)%>%
  summarise(.,Tipping=sum(Tipping),.groups = "keep")


d_all=d_com%>%
  rename(., Tipping_com=Tipping)%>%
  add_column(., Tipping_sp=d_sp$Tipping)

p1=ggplot(d_all)+
  geom_jitter(aes(x=Tipping_sp,y=Tipping_com,fill=as.factor(Competition),color=as.factor(Competition),shape=as.factor(Facilitation)),
              alpha=.7,size=.5,width=0.05, height=0.05)+
  geom_abline(slope = 1,intercept = 0,linetype=9)+
  facet_wrap(.~Branch,scales = "free")+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_shape_manual(values=c(21,22,24))+
  the_theme+
  labs(fill=TeX("$\\alpha_e$"),color=TeX("$\\alpha_e$"),shape="f",
       x="# Species scale tipping point",
       y="# Community scale tipping point")


d_ASS_sp=read.table("../Table/N_species/MF/d_ASS_sp.csv",sep=";")
N_replicate=40


#Species-specific tipping points size & frequency
p2=ggplot(d_ASS_sp, aes(x=Trait,y=Tipping/(N_replicate*3),fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f == .(Facilitation)))+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Fraction of tipping points",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))


d_ASS_sp=d_ASS_sp%>%
  group_by(.,Facilitation,Competition,Branch,Trait)%>%
  summarise(.,Size_shift=mean(Size_shift),.groups = "keep")

p3=ggplot(d_ASS_sp, aes(x=Trait,y=Size_shift,fill=as.factor(Competition))) + 
  geom_bar(position="stack", stat="identity")+
  facet_grid(Facilitation~Branch,labeller = label_bquote(rows = f == .(Facilitation)))+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Mean shift size (in fraction of cover)",fill=TeX("$\\alpha_e$ \ \ \ "))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))


Figure_4 = ggarrange(ggarrange(p2+theme(strip.background.y = element_blank(),strip.text.y = element_blank())
                               ,p3,ncol=2,labels=letters[1:2],common.legend = T,legend = "bottom"),
                     ggarrange(ggplot()+theme_void(),p0,ggplot()+theme_void(),ncol=3,labels = c("",letters[3],""),widths = c(1,8,1)),
                     nrow=2,heights = c(1.5,1),hjust = -3)
ggsave("../Figures/Final_figs/Figure_4.pdf",Figure_4,width = 9,height = 8)



# Figure 5 : Number of ASS ----
d=read.table('../Table/N_species/MF/Random_ini.csv',sep=";")

p1=ggplot(d%>%filter(., Branch=="Degradation"))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized),size=.8)+
  facet_grid(Competition~Facilitation,labeller = label_bquote(rows=alpha[e]==.(Competition), cols=f==.(Facilitation)))+
  the_theme+labs(y="Community index",color=expression(paste(bar(psi),"    ")))+ggtitle("Degradation")+
  scale_color_gradientn(colors = color_Nsp(Nsp),na.value = "black")




Nsp=15

d_hysteresis=tibble()
for (f in unique(d$Facilitation)){
  for (a0 in unique(d$Competition)){
    for (ini in unique(d$Random_ini)){
      
      d_h=filter(d,Facilitation==f,Competition==a0,Random_ini==ini)
      d_hysteresis=rbind(d_hysteresis,tibble(Facilitation=f,
                                             Competition=a0,
                                             Random_ini=ini,
                                             #We scale the hysteresis size with the mean cover to account for the effect of competition on global cover.
                                             Hysteresis=abs(sum(d_h$Rho_p[1:50]-rev(d_h$Rho_p[51:100])))/((sum(d_h$Rho_p[1:50]+rev(d_h$Rho_p[51:100])))/2))
      )
      
    }
  }
}


d_hysteresis_summary = d_hysteresis%>%
  group_by(., Competition,Facilitation)%>%
  summarise(.,.groups = "keep",Mean_hystersis=mean(Hysteresis),Sd_hysteresis=sd(Hysteresis))



p2=ggplot(NULL)+
  geom_jitter(data=d_hysteresis,aes(x=as.factor(Facilitation),y=Hysteresis,color=as.factor(Competition)),size=1,width = .1,alpha=.5)+
  geom_pointrange(data=d_hysteresis_summary,aes(x=as.factor(Facilitation),
                                                y=Mean_hystersis,
                                                ymin=Mean_hystersis-Sd_hysteresis,
                                                ymax=Mean_hystersis+Sd_hysteresis,
                                                fill=as.factor(Competition)),size=.5,color="black",shape=21)+
  the_theme+labs(y="Hysteresis size",color=TeX("$\\alpha_e$"),x="Facilitation, f",fill=TeX("$\\alpha_e$"))+
  scale_color_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))+
  scale_fill_manual(values=rev(c("#940000","#FF1F1F","#FFAFAF")))


ggsave("../Figures/Final_figs/Figure_5.pdf",ggarrange(p1,
                                                      ggarrange(ggplot()+theme_void(),p2,ggplot()+theme_void(),widths = c(.5,2,.5),ncol=3,labels = c("",letters[2],""))
                                                      ,nrow=2,heights = c(2,1),labels = c(letters[1],"")),width = 7,height = 8)

# Figure 6: Niche consequences ----

d_niche=read.table("../Table/2_species/PA/Niche_expansion_PA.csv",sep=";")
alpha_seq=seq(0,.3,length.out=10)
d_niche$alpha_0=alpha_seq
scale_facil="global_C_local_F"
p1=ggplot(d_niche%>%
           melt(.,measure.vars=c("Delta_niche_2"))%>%
           filter(., Scena==scale_facil,Local_disp==.1))+
  geom_tile(aes(x=Facilitation,y=alpha_0,fill=value))+
  the_theme+
  scale_fill_gradient2(low = "red",mid = "white",high = "blue")+
  labs(x="Facilitation ( f )",y=TeX(r'(Competition strength \ $\alpha_e)'),fill="% of competitive species niche change")


#niche varying traits
d_niche = read.table("../Table/2_species/PA/Niche_expansion_varying_traits.csv",sep=";")
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

Figure_6 = ggarrange(ggarrange(ggplot() + theme_void(),p1,ggplot() + theme_void(),widths = c(.7,2,.7),ncol=3,labels=c("",letters[1],"")),
                     p3,nrow=2,labels=c("",""))

ggsave("../Figures/Final_figs/Figure_6.pdf",Figure_6,width = 8,height = 8)

