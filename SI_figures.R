source('./Dryland_shift_functions.R')
library(grid)

# 2-species ----
## Advantage of CSI compared to Rho+----

list_f=list.files("../Table/2_species/PA/Multistability_PA/Frac_gradient/",pattern = "delta_0.1_facilitation_0.9_scalefacilitation_local")

d=tibble()
for (f in list_f){
  d2=read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/",f),sep=";")
  d=rbind(d,d2)
}

d[,1:2][d[,1:2] < 10^-4] = 0
set.seed(123)
u=runif(2)
d$CSI = sapply(1:nrow(d),function(x){
  return(u[1]*d$Stress_tolerant[x]+u[2]*d$Competitive[x])
})


d2_fil=filter(d,Psi1==1,Psi2==unique(d$Psi2)[41])
#Bifurcation diagram

p=ggplot(d2_fil%>%melt(., measure.vars=c("Rho_plus","CSI"))%>%
           mutate(Psi2=round(abs(Psi2),2)))+
  geom_line(aes(x = Stress, y = value, color = variable,linetype=Branches),lwd=.8) +
  labs(x = TeX(r'(Stress, \ $S)'), y = "", color = "",linetype="") +
  the_theme +
  scale_color_manual(values = c("black",alpha("red",.5)),labels=c("Global cover","Community index")) +
  scale_linetype_manual(values=c(1,9))+
  theme(axis.text.x = element_text(size = 9),axis.text.y = element_text(size = 9),strip.text.y = element_text(size=9))

ggsave("../Figures/Final_figs/SI/Using_CSI.pdf",p,width = 6,height = 4)

## Competition experienced ----
Nsp=15
trait=rev(seq(0,1,length.out=Nsp))
comp=comp2=c()
for (i in 1:Nsp){
  comp=c(comp,sum(sapply(1:Nsp,function(x){
    exp(-abs(trait[i]-trait[x]))
  })))
  comp2=c(comp2,sum(sapply(1:Nsp,function(x){
    trait[i]*exp(-abs(trait[i]-trait[x]))
  })))
}


p1=ggplot(tibble(Com=comp,Trait=trait,Comp2=comp2))+
  geom_line(aes(x=Trait,y=Com))+
  geom_point(aes(x=Trait,y=Com,color=Trait),size=3)+
  geom_text(data=tibble(x=c(0,trait[8],1),y=c(9.75,11.75,9.75),text=c("i","ii","iii")),aes(x=x,y=y,label=text),size=4)+
  the_theme+
  labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $\sum_j e^{-|\psi_i-\psi_j|})'),color="")+
  scale_color_gradientn(colours = color_Nsp(Nsp))+
  theme(legend.position = "none")

trait_subplot=c(0,trait[8],1)
for (i in 1:3){
  
  comp=sapply(1:Nsp,function(x){
    exp(-abs(trait_subplot[i]-trait[x]))
  })
  
  if (i != 3){
    assign(paste0("p_",i),
           ggplot(tibble(Com=comp,Trait=trait))+
             geom_line(aes(x=Trait,y=Com))+
             labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $e^{-|\psi_i-\psi_j|})'))+the_theme+
             theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.line.x = element_blank(),
                   axis.ticks.x = element_blank())+
             scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))
    )
  }else {
    assign(paste0("p_",i),
           ggplot(tibble(Com=comp,Trait=trait))+
             geom_line(aes(x=Trait,y=Com))+
             labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $e^{-|\psi_i-\psi_j|})'))+the_theme+
             scale_x_continuous(breaks = c(0,.5,1),labels = c("0","0.5","1"))
    )
  }
  
}



p2=ggplot(tibble(Com=comp,Trait=trait,Comp2=comp2))+
  geom_line(aes(x=Trait,y=Comp2))+
  geom_point(aes(x=Trait,y=Comp2,color=Trait),size=3)+
  the_theme+
  labs(x=TeX(r'(Trait species i, $\psi_i)'),y=TeX(r'(\ $\psi_i \sum_j e^{-|\psi_i-\psi_j|})'),color="")+
  scale_color_gradientn(colours = color_Nsp(15))

ggsave("../Figures/Final_figs/SI/Competition_experienced.pdf",
       ggarrange(
         ggarrange(p1,
                   ggarrange(p_1,p_2,p_3,nrow = 3,labels = c('i','ii','iii'),heights = c(1,1,1.4)),labels = c("a","","",""),widths = c(2,1)),
         p2,ncol=2,common.legend = T,legend = "none",labels = c("",letters[2]),widths = c(1.5,1)),
       width = 11,height = 4)



## Comparing CA & PA ----

#First layer = Simulations with PA (~Fig 2)

d2=read.table(paste0("../Table/2_species/PA/Multistability_fixed_traits_PA.csv"),sep=";")%>%
  filter(., Delta==.1)

d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Stress-tolerant")
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

#Second layer : CA simulations

list_f=list.files("../Table/2_species/CA/Justifying_use_PA/")
d=tibble()
for (i in list_f){
  d2=read.table(paste0("../Table/2_species/CA/Justifying_use_PA/",i),sep=",")
  mean_density=c(colMeans(d2[-c(1:1000),-1]),
                 str_split(i,pattern = "_")[[1]][3],
                 str_split(i,pattern = "_")[[1]][5],
                 gsub(pattern = ".csv",replacement = "",x=str_split(i,pattern = "_")[[1]][7]))
  d=rbind(d,mean_density)  
}
colnames(d)=c("Rho_1","Rho_2","Rho_0","Rho_d","Stress","Branches","alpha_0")
d=as_tibble(d)%>%
  mutate(., Rho_1=as.numeric(Rho_1),Rho_2=as.numeric(Rho_2),Rho_0=as.numeric(Rho_0),Rho_d=as.numeric(Rho_d),
         alpha_0=as.numeric(alpha_0),Stress=as.numeric(Stress))

d$state = sapply(1:nrow(d), function(x) {
  if (d[x, 1] > 0 & d[x, 2] > 0) {
    return("Coexistence")
  }
  if (d[x, 1] > 0 & d[x, 2] == 0) {
    return("Stress-tolerant")
  }
  if (d[x, 1] == 0 & d[x, 2] > 0) {
    return("Competitive")
  }
  if (d[x, 1] == 0 & d[x, 2] == 0) {
    return("Desert")
  }
})

d=d[order(d$alpha_0,d$Stress),]

all_state =sapply(seq(1, nrow(d) , by = 2),function(x){
  if (d$state[x] != d$state[x+1]){
    return(paste0(d$state[x],"/", d$state[x+1]))
  }
  else {return(d$state[x])}
})


d_state=d%>%
  filter(., Branches=="Degradation")%>%
  select(.,-Branches)
d_state$all_state=all_state

p=ggplot(NULL) +
  geom_tile(data=d2t,aes(x=Stress,y=alpha_0,fill=all_state))+
  geom_point(data=d_state,aes(x=Stress,y=alpha_0,fill=all_state),size=5,shape=21)+
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Coexistence/Competitive" = "#DDEFCA",
                             "Stress-tolerant" = "#BABEFF",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9))
ggsave("../Figures/Final_figs/SI/Comparizon_CA_PA.pdf",p,width = 7,height = 6)


## Multistability fixed traits MF model ----

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

d_state$Multistability=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c("Competitive","Coexistence","Desert","Stress_tolerant")){
    return("no_multi")
  } else{
    return("multi")
  }
})

p=ggplot(d_state) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',
                    color          = 'transparent',
                    pattern_spacing = 0.05, 
                    pattern_density = 0.01, 
                    pattern_fill    = alpha("black",.7),
                    pattern_colour  = alpha("black",.7))+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Stress (S)", y = TeX(r'(Strength of competition \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Stress_tolerant" = "#BABEFF",
                             "Stress_tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress_tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969",
                             "Competitive/Stress_tolerant"="#C998CE"
  ),labels=c("Stress_tolerant"="Stress-tolerant","Stress_tolerant/Desert"="Stress-tolerant/Desert",
             "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
             "Competitive/Stress_tolerant"="Competitive/Stress-tolerant"))+
  geom_hline(yintercept = round(unique(d2$alpha_0)[c(1,100)],4))+
  geom_text(data=tibble(x=c(1.025,1.025),y=c(0.01,.31),lab=c("c",'b')),aes(x=x,y=y,label=lab))+
  the_theme+theme(legend.text = element_text(size=10))+
  scale_pattern_manual(values=rev(c("transparent" ,"stripe")))


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


ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_MF.pdf",ggarrange(p+guides(pattern=F),p2,ncol = 2,widths = c(2,1),labels=c(letters[1],"")),
       width = 10,height = 6)




















## Multistability varying traits MF model ----

stress_seq=c(seq(0,.82, length.out = 100),seq(0.82,.9, length.out = 10)[-1])

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/MF/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",stress,"_facilitation_",.9,".csv"),sep=";"),
           read.table(paste0("../Table/2_species/MF/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",stress,"_facilitation_",.9,".csv"),sep=";"))
  d=rbind(d,d2)
  
  
  
}

d2=d
d2[,1:2][d2[,1:2] < 10^-3] = 0
d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Sp1")
  }
  if (d2[x, 1] == 0 & d2[x, 2] > 0) {
    return("Sp2")
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

d_state$multistab=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c( "Sp1/Sp2","Coexistence/Sp1","Sp2/Coexistence",
                                   "Sp2/Sp1")){
    return(1)
  } else {return(0)}
  
})

d_state$Sp1_reg=d2$Sp1[which(d2$Branches=="Restoration")]
d_state$Sp2_reg=d2$Sp2[which(d2$Branches=="Restoration")]

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      
      if (i>j){
        d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,
                                     Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                     Frac_multi_raw=sum(d_fil$multistab),
                                     Deg=d_fil$Stress[min(which(d_fil$Sp2==0))-1],
                                     Rest=d_fil$Stress[min(which(d_fil$Sp2_reg==0))]
        ))
      }
    }
  }
}

trait1_for_bifu=c(unique(d2$Psi1)[15],unique(d2$Psi1)[49],unique(d2$Psi1)[49])
trait2_for_bifu=c(unique(d2$Psi2)[13],unique(d2$Psi2)[1],unique(d2$Psi2)[49])

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),
               color="black",alpha=.3,breaks = c(.1,.15,.2))+
  geom_point(data=tibble(x=trait1_for_bifu,y=trait2_for_bifu),aes(x=x,y=y),shape=19)+
  geom_text(data=tibble(x=c(trait1_for_bifu[1],trait1_for_bifu[2]+.05,trait1_for_bifu[3]+.05),
                        y=c(trait2_for_bifu[1]+.13,trait2_for_bifu[2]+.05,trait2_for_bifu[3]+0.05),
                        txt=c("(ii)","(i)","(iii)")),aes(x=x,y=y,label=txt),size=4.5)+
  geom_segment(aes(x = 0.05, y = -0.02, xend = .95, yend = -0.02),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = .5, y = -0.06,label = "Sensitivity to competition"),stat = "unique")+
  geom_segment(aes(x = 1.02, y = 0, xend = 1.02, yend = .95),arrow = arrow(length = unit(0.3, "cm")))+
  geom_text(aes(x = 1.05, y = .45,label = "Trait similarity"),stat = "unique",angle=270)+
  the_theme+
  labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
       y=TeX(r'(Species 2 trait, $\psi_2$)'),
       fill="Fraction of community bistability \n along the stress gradient")+
  scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100))+
  guides(shape="none")


for (i in 1:3){
  
  if (i %in% c(1,3)){
    assign(paste0("p2_",i),
           ggplot(d2%>%
                    filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i])%>%
                    melt(.,measure.vars=c("Sp1","Sp2"))%>%
                    mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
             geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
             scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank(),
                   legend.position = "none")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             new_scale_color() +
             
             geom_line(data=d_state%>%
                         filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i]),
                       aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
             scale_color_manual(values=c("white","black"))+
             guides(color=F)
    )
    
  }else {
    assign(paste0("p2_",i),
           ggplot(d2%>%
                    filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i])%>%
                    melt(.,measure.vars=c("Sp1","Sp2"))%>%
                    mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
             geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
             scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
             the_theme+
             labs(x="Stress (S)",y="Species cover",color="",linetype="")+
             scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
             theme(legend.text = element_text(size=11))+
             new_scale_color() +
             
             geom_line(data=d_state%>%
                         filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i]),
                       aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
             scale_color_manual(values=c("white","black"))+
             guides(color=F)
           
           
    )
    
  }
}





pbifu=ggarrange(p2_2,p2_1,p2_3,
                ncol=3,
                labels = c("(i)","(ii)","(iii)"),common.legend = T,legend = "bottom",widths = c(1.2,1,1),hjust = c(-5,-1.3,-1.7),vjust=4)
p_tot=ggarrange(
  ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.15,1,.15)),
  pbifu,nrow=2,heights =  c(2,1),labels = c(letters[1:2]),hjust = c(-9,0))

ggsave(filename = "../Figures/Final_figs/SI/Multistability_varying_traits_MF.pdf",p_tot,width = 7,height = 7)






## Scaling-up 15sp, equal initial conditions ----


d_mf_equal=read.table("../Table/N_species/MF/MF_equal_ini.csv",sep=",")
d_mf_equal=d_mf_equal[,c(1:15,(ncol(d_mf_equal)-2):ncol(d_mf_equal))]
colnames(d_mf_equal)=c(paste0("Sp_",1:15),"alpha_0", "Branches", "Stress")
d_melt=d_mf_equal%>%
  melt(., measure.vars=paste0("Sp_",1:15))

d_melt$variable=as.character(d_melt$variable)

d_melt$Trait=sapply(1:nrow(d_melt),function(x){
  seq(0,1,length.out=15)[as.numeric(strsplit(d_melt$variable[x],split = "_")[[1]][2])]
})

p1=ggplot(d_melt%>%
            mutate(., Branches=recode_factor(Branches,"1"="Degradation","2"="Restoration"))) +
  geom_line(aes(x=Stress,y=value,color=Trait,group=interaction(Trait,alpha_0,Branches),linetype=as.factor(Branches)))+
  facet_wrap(.~alpha_0,labeller = label_bquote(cols=alpha[e]==.(alpha_0)),scales = "free")+
  the_theme+
  labs(x="Stress, (S)", y="Species cover",linetype="",color="")+
  scale_color_gradientn(colours = rev(color_Nsp(15)))+
  theme(strip.text.x = element_text(size=12))

ggsave("../Figures/Final_figs/SI/Equal_ini_MF.pdf",p1,width=8,height=3)





## Multistability along competition gradient with local competition ----
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


scale="local"
d2t=transform(d_state%>% 
                filter(., Scena %in% c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")))%>%
                mutate(., Stress=round(Stress,6),alpha_0=round(alpha_0,6)),
              Scena=factor(Scena,levels=c(paste0(scale,"_C_global_F"),paste0(scale,"_C_local_F")),
                           labels=c("Global facilitation","Local facilitation")))

color_rho = c("Coexistence" = "#D8CC7B", "Competitive" = "#ACD87B", "Desert" = "#696969", "Stress_tolerant" = "#7BD8D3")


d2t$Multistability=sapply(1:nrow(d2t),function(x){
  if (d2t$all_state[x] %in% c("Competitive","Coexistence","Desert","Stress_tolerant")){
    return("no_multi")
  } else{
    return("multi")
  }
})


Fig_2_SI=ggplot(d2t%>%
                    mutate(.,all_state=recode_factor(all_state,
                                                     "Competitive/Coexistence"="Coexistence/Competitive",
                                                     "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
                                                     "Stress_tolerant"="Stress-tolerant",
                                                     "Stress_tolerant/Desert"="Stress-tolerant/Desert"))%>%
                    mutate(., Delta=recode_factor(Delta,"0.1"="0.1, local dispersal","0.9"="0.9, global dispersal"))%>%
                    mutate(., Delta=as.character(Delta))) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',   
                    color          = 'transparent',
                    pattern_spacing = 0.05, 
                    pattern_density = 0.01, 
                    pattern_fill    = alpha("black",.7),
                    pattern_colour  = alpha("black",.7))+
  
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of competition, \ $\alpha_e)'), fill = "") +
  theme(legend.text = element_text(size = 11))+
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive" = "#ACD87B", 
                             "Coexistence/Competitive" = "#DDEFCA",
                             "Stress-tolerant" = "#BABEFF",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"))+
  facet_grid(Scena~Delta,labeller = label_bquote(cols = delta==.(Delta)))+
  the_theme+theme(strip.text.x = element_text(size=12),legend.text = element_text(size=9),strip.text.y = element_text(size=12))+
  scale_pattern_manual(values=rev(c("transparent" ,"stripe")))+
  guides(pattern="none")

ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_local_competition.pdf",Fig_2_SI,width = 8,height = 6)





## Multistability along competition gradient with global facilitation ----

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


d2t$Multistability=sapply(1:nrow(d2t),function(x){
  if (d2t$all_state[x] %in% c("Competitive","Coexistence","Desert","Stress_tolerant")){
    return("no_multi")
  } else{
    return("multi")
  }
})

set.seed(432)
rand_vec=runif(2)

d2$CSI=d2$Stress_tolerant*rand_vec[1]+d2$Competitive*rand_vec[2]

p1=ggplot(d2t%>%
            mutate(.,all_state=recode_factor(all_state,
                                             "Desert/Coexistence"="Coexistence/Desert",
                                             "Coexistence/Stress_tolerant"="Coexistence/Stress-tolerant",
                                             "Stress_tolerant"="Stress-tolerant",
                                             "Stress_tolerant/Desert"="Stress-tolerant/Desert",
                                             "Competitive/Stress_tolerant"="Competitive/Stress-tolerant"))%>%
            mutate(., Delta=recode_factor(Delta,"0.1"="0.1, local dispersal","0.9"="0.9, global dispersal"))%>%
            filter(., Scena=="Global facilitation")) +
  geom_tile_pattern(aes(x=Stress,y=alpha_0,fill=all_state,pattern=Multistability),              
                    colour          = 'transparent',
                    pattern_spacing = 0.05, 
                    pattern_density = 0.01, 
                    pattern_fill    = "black",
                    pattern_colour  = "black")+
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = TeX(r'(Stress, \ $S)'), y = TeX(r'(Strength of interspecific competition, \ $\alpha_e)'), fill = "") +
  scale_fill_manual(values=c("Coexistence" = "#D8CC7B", 
                             "Competitive/Stress-tolerant"="#C998CE",
                             "Competitive" = "#ACD87B", 
                             "Competitive/Coexistence" = "#DDEFCA",
                             "Stress-tolerant" = "#BABEFF",
                             "Stress-tolerant/Desert" ="#0F8E87",
                             "Coexistence/Stress-tolerant"="#9BBBB9",
                             "Coexistence/Desert"="#C19E5E",
                             "Desert"=  "#696969"
  ))+
  facet_grid(.~Delta,labeller = label_bquote(cols= delta==.(Delta) ))+
  the_theme+theme(strip.text.x = element_text(size=12),strip.text.y = element_text(size=11))+
  theme(legend.text = element_text(size = 9))+
  geom_text(data= tibble(Scena="Local facilitation",Delta="0.1, local dispersal",y=c(0.015,.415),x=1.05,labels=c("b","c")), 
            aes(x=x, y=y,label=labels), alpha=1, colour="white")+
  scale_pattern_manual(values=rev(c("none" ,"stripe")))+
  guides(pattern="none")




ggsave("../Figures/Final_figs/SI/Multistability_fixed_traits_global_facilitation.pdf",p1,width = 9,height = 5)


## Multistability varying traits with global facilitation ----

param_space=expand.grid(Scale_facil=c("local","global"),Disp=c(.1,.9))

for (num in 1:nrow(param_space)){

  list_f=list.files("../Table/2_species/PA/Multistability_PA/Frac_gradient",
                    pattern =paste0("delta_",param_space$Disp[num],"_facilitation_0.9_scalefacilitation_",param_space$Scale_facil[num]) )
  d=tibble()
  for (f in list_f){
    d=rbind(d,read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/",f),sep=";")%>%
              add_column(., Scale_facil=param_space$Scale_facil[num],Disp=param_space$Disp[num]))
  }


  d2=d
  d2$Stress=round(d2$Stress,4)
  d2[,1:2][d2[,1:2] < 10^-3] = 0
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
  d2=d2[order(d2$Psi2,d2$Stress,d2$alpha_0,d2$Psi1,d2$Disp,d2$Scale_facil),]

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

  d_state$multistab=sapply(1:nrow(d_state),function(x){
    if (d_state$all_state[x] %in% c( "Species 2/Stress_tolerant","Coexistence/Stress_tolerant","Species 2/Coexistence",
                                     "Stress_tolerant/Species 2")){
      return(1)
    } else {return(0)}

  })

  d_state$Stress_tolerant_reg=d2$Stress_tolerant[which(d2$Branches=="Restoration")]
  d_state$Competitive_reg=d2$Competitive[which(d2$Branches=="Restoration")]

  d_final=tibble()
  for (i in unique(d_state$Psi1)){
    for (j in unique(d_state$Psi2)){
      for (a0 in unique(d_state$alpha_0)){
        d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)

        if (i>j){
          d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,Disp=param_space$Disp[num],
                                       Scale_facil=param_space$Scale_facil[num],
                                       Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                       Frac_multi_raw=sum(d_fil$multistab),
                                       Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
                                       Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))]
          ))
        }


      }
    }
  }

  assign(paste0("p_",num),
         ggplot(d_final)+
           geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
           the_theme+
           labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
                y=TeX(r'(Species 2 trait, $\psi_2$)'),
                fill="Fraction of community bistability \n along the stress gradient")+
           scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100))+
           guides(shape="none")+
           ggtitle(paste0(ifelse(param_space$Scale_facil[num]=="global","Global","Local")," facilitation, ",
                          ifelse(param_space$Disp[num]==.1,"local","global")," dispersal"))
         
         
         )
  
  
  
}

p=ggarrange(p_1+labs(fill="")+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()),
            p_2+labs(fill="")+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank(),
                    axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            p_3+labs(fill="")+scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100),breaks=c(0,.1,.2)),
            p_4+labs(fill="")+scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100),breaks=c(0,.1,.2))+
              theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            nrow=2,ncol=2,heights = c(1,1.1))

ggsave("../Figures/Final_figs/SI/Multistability_traits_scale_facil_disp.pdf",p,width = 7,height = 7)

## Mechanisms for global facilitation: clustering, pairs, thresholds ----



#invasion and extinction
d_invade=read.table("../Table/2_species/PA/Threshold_invasion.csv",sep=";")
d_extinction=read.table("../Table/2_species/PA/Threshold_extinction.csv",sep=";")
colnames(d_invade)=colnames(d_extinction)
d_threshold=rbind(d_invade%>%add_column(., Branch="Restoration"),
                  d_extinction%>%add_column(., Branch="Desertification"))

d_ribbon=tibble(x=d_invade$delta,ymax=d_invade$Thresh_extinction,
                ymin=d_extinction$Thresh_extinction,Scena=d_invade$Scena)

p1=ggplot(d_threshold%>%filter(., Scena==1))+
  geom_point(aes(x=delta,y=Thresh_extinction,shape=Branch),size=3)+
  geom_line(aes(x=delta,y=Thresh_extinction,group=Branch),lwd=.5)+
  #geom_ribbon(data=d_ribbon%>%filter(., Scena==1),aes(x=x,ymin=ymin,ymax=ymax),fill="#0F8E87", alpha=0.15)+
  labs(x=TeX('$\\delta$'),y="Threshold of stress",color="",shape="")+
  the_theme+
  scale_shape_manual(values=c(1,8))+
  guides(shape=guide_legend(ncol=1))+
  ggtitle(TeX(r'($\alpha_e = 0.2)'))




# clustering of stress-tolerant species
d_clustering = read.table("../Table/2_species/PA/Clustering_PA.csv",sep=";")


d_clustering$c12=d_clustering$c11=d_clustering$c22=0
name_scena=c("global_C_global_F","global_C_local_F")

for (nr in 1:nrow(d_clustering)){
  
  rho_2=d_clustering$Rho_2[nr]
  rho_1=d_clustering$Rho_1[nr]
  rho_12=d_clustering$Rho_12[nr]
  rho_11=d_clustering$Rho_11[nr]
  rho_22=d_clustering$Rho_22[nr]
  
  if (rho_2 < 10^-2) rho_2 = 0#if (rho_2 < 10^-4) rho_2 = 0
  if (rho_1 < 10^-2) rho_1 = 0#if (rho_1 < 10^-4) rho_1 = 0
  
  if (rho_1 > 0 && rho_2 > 0){
    q2_1 = rho_12 / rho_1 #average number of species 2 around species 1 
    c21 = q2_1 / rho_2
    q1_2 = rho_12 / rho_2 #average number of species 1 around species 2
    c12 = q1_2 / rho_1
    
  } else{
    q1_2 = 0
    c12 = NaN
    q2_1 = 0
    c21 = NaN
  }
  
  if (rho_2 > 0){
    q2_2 = rho_22 / rho_2 #average number of species 2 around species 2 
    c22 = q2_2 / rho_2
  } else{
    q2_2 = 0
    c22 = NaN
  }
  
  if (rho_1 > 0){
    q1_1 = rho_11 / rho_1 #average number of species 1 around species 1 
    c11 = q1_1 / rho_1
  } else{
    q1_1 = 0
    c11 = NaN
  }
  
  d_clustering$c12[nr]=c12
  d_clustering$c22[nr]=c22
  d_clustering$c11[nr]=c11
  d_clustering$Scena[nr]=name_scena[as.numeric(d_clustering$Scena[nr])]
  
  
}


d_clustering=mutate(d_clustering, Scena=recode_factor(Scena,"global_C_global_F"="Global facilitation","global_C_local_F"="Local facilitation"))


d_clust_mean=d_clustering%>%
  filter(., alpha_0==.2, S %in% c(.73),Scena=="Global facilitation",
         delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
  group_by(., cintra,alpha_0,S,delta,Scena)%>%
  summarise(.,.groups ="keep",c11=mean(c11) )


p2=ggplot(d_clust_mean)+
  geom_point(aes(x=delta,y=c11),size=3,shape=1)+
  geom_line(aes(x=delta,y=c11,group=interaction(S)),lwd=.5)+
  the_theme+labs(x=TeX('$\\delta$'),color="",y="Stress-tolerant clustering")+
  geom_hline(yintercept = 1)+
  ggtitle(TeX(r'($\S = 0.73, \alpha_e = 0.2)'))


# density of pairs
d_pairs=d_clustering%>%
  filter(., alpha_0==.2, S %in% c(0),Scena=="Global facilitation",
         delta %in% unique(d_clustering$delta)[seq(1,length(unique(d_clustering$delta)),length.out=15)])%>%
  melt(., measure.vars=c("Rho_20","Rho_10"))%>%
  mutate(.,variable=recode_factor(variable,"Rho_10"="Pair Stress-tol. sp./Fertile",
                                  "Rho_20"="Pair Competitive sp./Fertile"))
p3=ggplot(d_pairs)+
  geom_point(aes(x=delta,y=value,color=variable),size=3,shape=1)+
  geom_line(aes(x=delta,y=value,color=variable,group=interaction(variable)),lwd=.5)+
  ggtitle(TeX(r'($\S = 0, \alpha_e = 0.2)'))+
  the_theme+
  labs(x=TeX('$\\delta$'),y=expression(paste("Pair cover")),color="")+
  scale_color_manual(values=c("#858BE0",as.character(color_Nsp(5)[c(2)])))+
  theme(legend.box = "vertical")+
  guides(color=guide_legend(ncol=1))

p_bottom=ggarrange(p3,
                   ggarrange(p2,ggplot()+theme_void(),nrow=2,heights = c(1,.25)),
                   p1,ncol=3,labels = letters[1:3],align = "h")

ggsave(filename = "../Figures/Final_figs/SI/Mechanisms_global_facilitation.pdf",width = 7,height = 4)


## Multistability, varying the trade-off shape ----

#first: explaining trade-off


d=expand.grid(Psi=seq(0,1,length.out=20),Gamma=c(.5,1,1.5))
d$Psi2=sapply(1:nrow(d),function(x){
  return(d$Psi[x]^d$Gamma[x])
})

p=ggplot(d)+
  geom_line(aes(x=Psi,y=Psi2,color=as.factor(Gamma)))+
  geom_point(aes(x=Psi,y=Psi2,color=as.factor(Gamma)),shape=1,size=3,fill="white")+
  geom_segment(aes(x = .5, y = .15,
                   xend = 1, yend = .15),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(aes(x = .75, y = .22,label="Stress-tolerance"),stat = "unique")+
  geom_text(aes(x = .75, y = .08,label="Facilitation"),stat = "unique")+
  geom_segment(aes(x = 0.5, y = .8,
                   xend = 0, yend = .8),
               arrow = arrow(length = unit(0.2, "cm")))+
  geom_text(aes(x = .25, y = .9,label="Tolerance to competition"),stat = "unique")+
  the_theme+
  labs(x=TeX("$\\psi$"),y=TeX("$\\psi^{\\gamma}$"),color=TeX("$\\gamma$"))+
  scale_color_viridis_d()
ggsave( "../Figures/Final_figs/SI/Trade_off_shape.pdf",p,width = 6,height = 4)



#second: multistability
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()

for (tradeoff in c(.5,1.5)){
  for (stress in stress_seq){
    
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Degradation_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Restoration_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"))
    d=rbind(d,d2)
    
    
    
  }
}

d2=d
colnames(d2)[1:2]=c("Sp1","Sp2")
d2[,1:2][d2[,1:2] < 10^-2] = 0


d_state=d2%>%
  filter(., Branches=="Degradation")%>%
  select(.,-Branches)

d_state$multistab=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c( "Sp1/Sp2","Coexistence/Sp1","Coexistence/Sp2",
                                   "Sp2/Coexistence","Sp1/Coexistence",
                                   "Sp2/Sp1")){
    return(1)
  } else {return(0)}
  
})

d_state$Sp1_reg=d2$Sp1[which(d2$Branches=="Restoration")]
d_state$Sp2_reg=d2$Sp2[which(d2$Branches=="Restoration")]

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      for (tradeoff in unique(d_state$Tradeoff)){
        d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0,Tradeoff==tradeoff)
        
        if (i>j){
          d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,Tradeoff=tradeoff,
                                       Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                       Frac_multi_raw=sum(d_fil$multistab),
                                       Deg=d_fil$Stress[min(which(d_fil$Sp2==0))-1],
                                       Rest=d_fil$Stress[min(which(d_fil$Sp2_reg==0))]
          ))
        }
      }
    }
  }
}

trait1_for_bifu=c(unique(d_state$Psi1)[49],unique(d_state$Psi1)[13])
trait2_for_bifu=c(unique(d_state$Psi2)[1],unique(d_state$Psi2)[12])

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
  geom_contour(aes(x=Psi1,y=Psi2,z=Frac_multi_raw/length(unique(d_state$Stress))),
               color="black",alpha=.3,breaks = c(0,.2,.4))+
  geom_point(data=tibble(x=trait1_for_bifu,y=trait2_for_bifu),aes(x=x,y=y),shape=19)+
  geom_text(data=tibble(x=c(trait1_for_bifu[1],trait1_for_bifu[2]+.02),
                        y=c(trait2_for_bifu[1]+.06,trait2_for_bifu[2]+.12),
                        txt=c("(i)","(ii)")),aes(x=x,y=y,label=txt),size=4.5)+
  the_theme+
  facet_wrap(.~Tradeoff,labeller = label_bquote(cols = gamma == .(Tradeoff)))+
  labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
       y=TeX(r'(Species 2 trait, $\psi_2$)'),
       fill="Fraction of community bistability \n along the stress gradient")+
  scale_fill_gradientn(colors=colorRampPalette(c("#FFFFFF","#E2BFE2","#D48ACF","#650328"))(100))+
  guides(shape="none")+
  theme(strip.text.x = element_text(size=13))


for (i in 1:2){
  
  assign(paste0("p2_",i),
         ggplot(d2%>%
                  filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==0.5)%>%
                  melt(.,measure.vars=c("Sp1","Sp2"))%>%
                  mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
           geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
           scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
           the_theme+
           labs(x="Stress (S)",y="Species cover",color="",linetype="")+
           scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
           theme(legend.text = element_text(size=11))+
           guides(color = guide_legend(override.aes = list(size = 1.5)))+
           new_scale_color() +
           
           geom_line(data=d_state%>%
                       filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==0.5),
                     aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
           scale_color_manual(values=c("white","black"))+
           guides(color=F)
  )
  
  assign(paste0("p3_",i),
         ggplot(d2%>%
                  filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==1.5)%>%
                  melt(.,measure.vars=c("Sp1","Sp2"))%>%
                  mutate(., variable=recode_factor(variable,"Sp1"="Sp. 1","Sp2"="Sp. 2")))+
           geom_line(aes(x=Stress,y=value,color=variable,linetype=Branches))+
           scale_color_manual(values=as.character(color_Nsp(5)[c(2,5)]))+
           the_theme+
           labs(x="Stress (S)",y="Species cover",color="",linetype="")+
           scale_y_continuous(breaks = c(0,.2,.4,.6,.8))+
           theme(legend.text = element_text(size=11))+
           guides(color = guide_legend(override.aes = list(size = 1.5)))+
           new_scale_color() +
           
           geom_line(data=d_state%>%
                       filter(., Psi2==trait2_for_bifu[i],Psi1==trait1_for_bifu[i],Tradeoff==1.5),
                     aes(x=Stress,y=.85,color=as.factor(multistab),group=as.factor(multistab)),alpha=.3,lwd=2)+
           scale_color_manual(values=c("white","black"))+
           guides(color=F)
  )  
    
  
}





pbifu=ggarrange(ggarrange(ggarrange(p2_1,p2_2+ylab(""), labels = c("(i)","(ii)"),
                          ncol=2,widths = c(1.2,1),hjust = c(-4.5,-3.5),vjust=4,
                          legend="none"),
                ggarrange(p3_1+ylab(""),p3_2+ylab(""), labels = c("(i)","(ii)"),
                          ncol=2,widths = c(1,1),hjust = c(-4.5,-6),vjust=4,
                          legend="none"),
                ncol=2,
                legend = "none",
                widths = c(1.1,1),labels = c("b","c")),
                ggarrange(ggplot()+theme_void(),
                          get_legend(p2_1),
                          ggplot()+theme_void(),ncol=3,widths = c(.5,2,.5)),nrow=2,heights = c(4,1))


p_tot=ggarrange(
  ggarrange(ggplot()+theme_void(),p1,ggplot()+theme_void(),ncol=3,widths = c(.15,1,.15)),
  pbifu,nrow=2,heights =  c(2,1),labels = c(letters[1],""),hjust = c(-9,0))

ggsave(filename = "../Figures/Final_figs/SI/Trade_off_multistability_PA.pdf",plot = p_tot,width = 9,height = 7)





## Comparing degradation and restoration points trade-off ----

#convex and concave trade-off
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()

for (tradeoff in c(.5,1.5)){
  for (stress in stress_seq){
    
    
    
    d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Degradation_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"),
             read.table(paste0("../Table/2_species/PA/Multistability_PA/Varying_tradeoff/Test_interspe_comp_",
                               .3,"_branch_Restoration_stress_",round(stress,4),"_tradeoff_",tradeoff,"_facilitation_",.9,".csv"),sep=";"))
    d=rbind(d,d2)
    
    
    
  }
}

d2=d
colnames(d2)[1:2]=c("Sp1","Sp2")
d2[,1:2][d2[,1:2] < 10^-2] = 0
d2$state = sapply(1:nrow(d2), function(x) {
  if (d2[x, 1] > 0 & d2[x, 2] > 0) {
    return("Coexistence")
  }
  if (d2[x, 1] > 0 & d2[x, 2] == 0) {
    return("Sp1")
  }
  if (d2[x, 1] == 0 & d2[x, 2] > 0) {
    return("Sp2")
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

d_state$multistab=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c( "Sp1/Sp2","Coexistence/Sp1","Coexistence/Sp2",
                                   "Sp2/Coexistence","Sp1/Coexistence",
                                   "Sp2/Sp1")){
    return(1)
  } else {return(0)}
  
})

d_state$Sp1_reg=d2$Sp1[which(d2$Branches=="Restoration")]
d_state$Sp2_reg=d2$Sp2[which(d2$Branches=="Restoration")]

d_final1=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      for (tradeoff in unique(d_state$Tradeoff)){
        d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0,Tradeoff==tradeoff)
        if (i>j){
          d_final1=rbind(d_final1,tibble(Psi1=i,Psi2=j,alpha_0=a0,Tradeoff=tradeoff,
                                       Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                       Frac_multi_raw=sum(d_fil$multistab),
                                       Deg=d_fil$Stress[min(which(d_fil$Sp2==0))-1],
                                       Rest=d_fil$Stress[min(which(d_fil$Sp2_reg==0))]
          ))
        }
      }
    }
  }
}



#linear trade-off
stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",round(stress,4),"_delta_",.1,"_facilitation_0.9_scalefacilitation_local.csv"),sep=";"),
           read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",round(stress,4),"_delta_",.1,"_facilitation_0.9_scalefacilitation_local.csv"),sep=";"))
  d=rbind(d,d2)
  
  
  
}



d2=d
d2$Stress=round(d2$Stress,4)
d2[,1:2][d2[,1:2] < 10^-3] = 0
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

d_state$multistab=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c( "Species 2/Stress_tolerant","Coexistence/Stress_tolerant","Species 2/Coexistence",
                                   "Stress_tolerant/Species 2")){
    return(1)
  } else {return(0)}
  
})

d_state$Stress_tolerant_reg=d2$Stress_tolerant[which(d2$Branches=="Restoration")]
d_state$Competitive_reg=d2$Competitive[which(d2$Branches=="Restoration")]

d_final2=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      if (i>j){
        d_final2=rbind(d_final2,tibble(Psi1=i,Psi2=j,alpha_0=a0,
                                     Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                     Frac_multi_raw=sum(d_fil$multistab),
                                     Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
                                     Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))]
        ))
      }
    }
  }
}

#now comparing degradation and restoration points
#making sure that they are all sorted in the same way

d_final2=d_final2[order(d_final2$Psi1,d_final2$Psi2,d_final2$alpha_0),]
d_convex=filter(d_final1,Tradeoff==.5)
d_concave=filter(d_final1,Tradeoff==1.5)

d_final2=d_final2[order(d_final2$Psi1,d_final2$Psi2,d_final2$alpha_0),]


#now making the difference between non-linear trade-off and linear trade-off

d_convex$Deg=d_convex$Deg-d_final2$Deg
d_concave$Deg=d_concave$Deg-d_final2$Deg
d_convex$Rest=d_convex$Rest-d_final2$Rest
d_concave$Rest=d_concave$Rest-d_final2$Rest
d_convex$Frac_multi_raw=d_convex$Frac_multi_raw-d_final2$Frac_multi_raw
d_concave$Frac_multi_raw=d_concave$Frac_multi_raw-d_final2$Frac_multi_raw

d_tot=rbind(d_convex,d_concave)

#now plotting

for (i in 1:2){
  assign(paste0("p1_",i),
         
         ggplot(d_tot%>%filter(., Tradeoff==c(.5,1.5)[i]))+
           geom_tile(aes(x=Psi1,y=Psi2,fill=Frac_multi_raw/length(unique(d_state$Stress))))+
           the_theme+
           facet_grid(.~Tradeoff,labeller = label_bquote(cols = gamma == .(Tradeoff)))+
           labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
                y=TeX(r'(Species 2 trait, $\psi_2$)'),
                fill="")+
           scale_fill_gradient2(low = "red",mid="white",high="blue",midpoint = 0)+
           guides(shape="none")+
           theme(strip.text.x = element_text(size=13),strip.background.x = element_blank())
  )
}





param_space=expand.grid(point=c("Degradation point","Restoration point"),
                        shape=c(.5,1.5))

for (i in 1:4){
  
  assign(paste0("p2_",i),
         ggplot(d_tot%>%
                  melt(., measure.vars=c("Rest","Deg"))%>%
                  mutate(., variable=recode_factor(variable,
                                                   "Rest" = "Restoration point",
                                                   "Deg" = "Degradation point"))%>%
                  filter(., variable==param_space$point[i],Tradeoff==param_space$shape[i]))+
           geom_tile(aes(x=Psi1,y=Psi2,fill=value))+
           the_theme+
           facet_grid(.~Tradeoff,labeller = label_bquote(cols = gamma == .(Tradeoff)))+
           labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),
                y=TeX(r'(Species 2 trait, $\psi_2$)'),
                fill="")+
           scale_fill_gradient2(low = "red",mid=alpha("white",.5),high="blue")+
           guides(shape="none")+
           theme(strip.text.x = element_text(size=13),strip.background.x = element_blank())
  )
  
}


p_tot=ggarrange(
  ggarrange(p1_1+ggtitle("Community bistability"),p1_2,labels = c(letters[1],""),align = "hv"),
  ggarrange(p2_1+ggtitle("Extinction point"),p2_3,labels = c(letters[2],""),align = "hv"),
  ggarrange(p2_2+ggtitle("Restoration point"),p2_4,labels = c(letters[3],""),align = "hv"),
  nrow=3
)


ggsave("../Figures/Final_figs/SI/Comparizon_linear_non_linear_tradeoff.pdf",width = 7,height = 10)






## Restoration & extinction points trait space----


stress_seq=seq(0,.82, length.out = 100)

d=tibble()  
d_bistab=tibble()


for (stress in stress_seq){
  
  
  
  d2=rbind(read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Degradation_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"),
           read.table(paste0("../Table/2_species/PA/Multistability_PA/Frac_gradient/Test_interspe_comp_",
                             .3,"_branch_Restoration_stress_",stress,"_delta_",.1,"_facilitation_",.9,".csv"),sep=";"))
  d=rbind(d,d2)
  
  
  
}

d2=d
d2[,1:2][d2[,1:2] < 10^-3] = 0
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

d_state$multistab=sapply(1:nrow(d_state),function(x){
  if (d_state$all_state[x] %in% c( "Species 2/Stress_tolerant","Coexistence/Stress_tolerant","Species 2/Coexistence",
                                   "Stress_tolerant/Species 2")){
    return(1)
  } else {return(0)}
  
})

d_state$Stress_tolerant_reg=d2$Stress_tolerant[which(d2$Branches=="Restoration")]
d_state$Competitive_reg=d2$Competitive[which(d2$Branches=="Restoration")]

d_final=tibble()
for (i in unique(d_state$Psi1)){
  for (j in unique(d_state$Psi2)){
    for (a0 in unique(d_state$alpha_0)){
      d_fil=filter(d_state,Psi1==i,Psi2==j,alpha_0==a0)
      
      if (i>j){
        d_final=rbind(d_final,tibble(Psi1=i,Psi2=j,alpha_0=a0,
                                     Frac_multi=sum(d_fil$multistab)/length(which(d_fil$state!="Desert")),
                                     Frac_multi_raw=sum(d_fil$multistab),
                                     Deg=d_fil$Stress[min(which(d_fil$Competitive==0))-1],
                                     Deg2=d_fil$Stress[max(which(d_fil$Stress_tolerant>0))+1],
                                     Rest=d_fil$Stress[min(which(d_fil$Competitive_reg==0))],
                                     Niche_1 = length(which(d_fil$Stress_tolerant>0))/length(which(d_fil$state!="Desert"))  ,
                                     Niche_2 = length(which(d_fil$Competitive>0))/length(which(d_fil$state!="Desert"))
        ))
      }
    }
  }
}



# Niche least competitive species
p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Niche_1))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_2$)'),
                 fill="Fraction of vegetation states \n where species 1 has a positive cover")+
  scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100))+
  guides(shape=F)


ggsave("../Figures/Final_figs/SI/Niche_least_competitive_species.pdf",p1,width = 6,height = 5)


#Degradation and restoration points

p1=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Deg))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_1$)'),fill="Species 2 extinction point")+
  scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100))+
  guides(shape=F)



#Restoration point traits
p2=ggplot(d_final)+
  geom_tile(aes(x=Psi1,y=Psi2,fill=Rest))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y=TeX(r'(Species 2 trait, $\psi_2$)'),fill="Species 2 restoration point")+
  scale_fill_gradientn(colors=colorRampPalette(c("#F9F4E7","#E6CC8B","#E2A472","#B50F02"))(100),
                       breaks=c(0,.1,.2))+
  geom_line(data=tibble(x=c(0,1,1,1),y=c(0,0,0,1),group=c(1,1,2,2)),aes(x=x,y=y,group=group),lwd=1)+
  geom_text(data=tibble(x=c(.5,1.05),
                        y=c(0.05,.55),
                        txt=c("(c)","(d)")),aes(x=x,y=y,label=txt),size=4.5)+
  guides(shape=F)



p3=ggplot(d_final%>%filter(Psi2==0))+
  geom_smooth(aes(x=Psi1,y=Rest),se = F,color="black")+ggtitle(TeX("$\\psi_2 = 0$"))+
  the_theme+labs(x=TeX(r'(Species 1 trait, $\psi_1$)'),y="Species 2 restoration point")

p4=ggplot(d_final%>%filter(Psi1==unique(d_final$Psi1)[49]))+
  geom_smooth(aes(x=Psi2,y=Rest),se = F,color="black")+ggtitle(TeX("$\\psi_1 = 1$"))+
  the_theme+labs(x=TeX(r'(Species 2 trait, $\psi_2$)'),y="Species 2 restoration point")

p_bottom=ggarrange(p3,p4,labels = letters[3:4])
p_tot=ggarrange(ggarrange(p1,p2,ncol=2,labels = letters[1:2]),
                p_bottom,nrow=2,heights = c(1.5,1))

ggsave("../Figures/Final_figs/SI/Restoration_degration.pdf",p_tot,width = 9,height = 7)










## Landscapes along dispersal gradient ----

list_landscape=list.files("../Table/2_species/CA/",pattern = "Dispersal")
for (i in 1:length(list_landscape)){
  landscape=read.table(paste0("../Table/2_species/CA/",list_landscape[i]),sep=",")
  assign(paste0("p_land_",i),Plot_landscape(landscape,Nsp = 2)+theme(legend.position = "none",
                                                                     plot.title = element_text(size=12)))
}
p=ggarrange(p_land_1+ggtitle(TeX("$\\delta = 0$")),p_land_2+ggtitle(TeX("$\\delta = 0.1$")),
          p_land_3+ggtitle(TeX("$\\delta = 0.2$")),p_land_4+ggtitle(TeX("$\\delta = 0.3$")),
          p_land_5+ggtitle(TeX("$\\delta = 0.7$")),p_land_6+ggtitle(TeX("$\\delta = 1$")),nrow=2,ncol=3)


ggsave("../Figures/Final_figs/SI/Dynamics_landscape_dispersal.pdf",p,width = 7,height = 4)

## Propensity of colonization along dispersal gradient ----

list_file_prop=list.files('../Table/2_species/CA/Propensity_colonization')

d_propensity=tibble()
N_rep=1

for (f in list_file_prop){
  
  
  
  landscape=read.table(paste0("../Table/2_species/CA/Propensity_colonization/",f),sep=",")
  save_landscape=landscape
  
  delta=as.numeric(strsplit(x = f,split = "_")[[1]][4])
  
  #propensity colonization fertile species 1 = \beta (\delta * \rho_1 + (1-\delta) * q_{1|0}) in each fertile sites
  #yet this metric depends on the vegetation cover and the number sites colonized. 
  #Therefore, we take the mean (dividing by the total number of site 1)
  
  landscape[landscape != 1] = 0
  neigh_1= simecol::neighbors(x =as.matrix(landscape),state = 1, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
  propensity_species_1 = sum((delta * sum(landscape)/length(landscape) + (1-delta) * neigh_1[which(save_landscape==0)]))
  propensity_species_1=propensity_species_1/(length(which(save_landscape==1))*length(which(save_landscape==0)))
  
  #propensity colonization fertile species 2
  
  landscape=save_landscape
  landscape[landscape != 2] = 0
  neigh_2= simecol::neighbors(x =as.matrix(landscape),state = 2, wdist =  matrix( c(0, 1, 0,1, 0, 1, 0, 1, 0), nrow = 3),bounds = 1)
  propensity_species_2 = sum(delta * sum(landscape)/length(landscape) + (1-delta) * neigh_2[which(save_landscape==0)])
  propensity_species_2=propensity_species_2/(length(which(save_landscape==2))*length(which(save_landscape==0)))
  
  d_propensity=rbind(d_propensity,tibble(Delta=delta,Replicate=N_rep,
                                         Sp1_cover=length(which(save_landscape==1))/(length(landscape)),
                                         Sp2_cover=length(which(save_landscape==2))/(length(landscape)),
                                         Prop_sp1=propensity_species_1,
                                         Prop_sp2=propensity_species_2))
  
  if (N_rep==30){
    N_rep=0
  }
    
  N_rep=N_rep+1
}

d_averaged=d_propensity%>%
  group_by(., Delta)%>%
  summarise(., .groups = "keep",Prop1_mean=mean(Prop_sp1),Prop2_mean=mean(Prop_sp2),
            mean_Sp1=mean(Sp1_cover),mean_Sp2=mean(Sp2_cover))

ggplot(NULL)+
  geom_point(data=d_averaged%>%
               melt(., measure.vars=c("Prop1_mean","Prop2_mean")),
             aes(x=Delta,y=value,color=variable))+
  geom_line(data=d_averaged%>%
               melt(., measure.vars=c("Prop1_mean","Prop2_mean")),
             aes(x=Delta,y=value,color=variable))+
  the_theme+
  scale_y_log10()
  


  
  






# N-species ----
## Main fig Nspecies with mean-field model ----
d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

p1=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.35),Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.1,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())

p2=ggplot(d_tot%>%filter(., Nsp %in% c(5,25),Competition %in% c(.225,.35),Branch==1))+
  geom_point(aes(x=Stress,y=Rho_plus,color=Psi_normalized,fill=Psi_normalized),size=.3,shape=21)+
  the_theme+labs(y="Vegetation cover",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S",color="Mean community trait  ")+
  facet_grid(Competition~Nsp,labeller = label_bquote(rows= alpha[e] ==.(Competition),cols="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=10),strip.text.y = element_text(size=11),
        panel.background = element_blank(),strip.background.y = element_blank())


d_richness=read.table("../Table/N_species/MF/Multistability_richness.csv",sep=";")


d_richness_mean=d_richness%>%
  group_by(., Nsp,Branch,Competition)%>%
  summarise(., mean_richness=mean(Richness),.groups = "keep")

p3=ggplot(d_richness%>%filter(., Branch=="Degradation"))+
  geom_bar(aes(x=Competition,fill=as.factor(Richness),color=as.factor(Richness),y=Richness),
           color="transparent",position = "fill",stat="identity")+
  geom_point(data=d_richness_mean,aes(x=Competition,y=mean_richness/7),
             shape=1,size=3,color="white")+
  facet_wrap(.~Nsp,labeller = label_bquote(cols="Number of species"==.(Nsp)))+
  labs(x=TeX(r'(Strength of interspecific competition, \ $\alpha_e)'),
       y="Fraction of initial conditions",fill="Number of species \n with positive cover")+
  scale_fill_viridis_d()+
  scale_color_viridis_d()+
  the_theme+
  guides(color="none")+
  scale_y_continuous(sec.axis=sec_axis(~.*7, name=TeX(r'(Mean community \n trait, \ $\bar \psi)')))+
  theme(strip.text.x = element_text(size=11))+
  scale_x_continuous(breaks = seq(0.225,.35,by=.025),labels = c("0.225","0.25",'0.275','0.3','0.325','0.35'))


Fig_6_SI=ggarrange(ggarrange(p2,p1,common.legend = T,legend = "bottom",ncol=2,labels = letters[1:2]),
                ggarrange(ggplot()+theme_void(),p3,ggplot()+theme_void(),ncol=3,widths = c(.05,1.5,.05),labels = c("",letters[3],"")),nrow=2,heights = c(1,1.2))

ggsave("../Figures/Final_figs/SI/N_species_CSI_cover_MF.pdf",Fig_6_SI,width = 9,height = 8)





## Example succession species


## Example of species succession gradient stress Nspecies ----


Fig_to_plot=tibble(Nsp=c(5,15,15,15,25),Number_plot=c(3,169,214,737,397))


for (i in 1:nrow(Fig_to_plot)){

  Nsp=Fig_to_plot$Nsp[i]
  list_csv=list.files(paste0('../Table/N_species/MF/',Nsp,'_sp/'))
  d2=read.table(paste0("../Table/N_species/MF/",Nsp,"_sp/",list_csv[Fig_to_plot$Number_plot[i]+1]),sep=",")
  colnames(d2)=c(paste0(1:Nsp),"Fertile","Degraded","Random_ini","Competition","Dispersal","Facilitation","Branch","Stress")
  
  d2[d2<10^(-4)]=0
  
    
  
  assign(paste0("p_",i),
    ggplot(d2%>%
             melt(., measure.vars=paste0(1:Nsp))%>%
             mutate(., Trait=sapply(1:nrow(.),function(x){
               return(seq(1,0,length.out=Nsp)[as.numeric(.$variable[x])])
             })))+
      geom_line(aes(x=Stress,y=value,color=Trait,group=variable),size=.8)+
      the_theme+labs(x="Stress, S",y="Species cover",color=TeX("$\\psi$   "))+
      scale_color_gradientn(colors=color_Nsp(Nsp))
  )
}

p=ggarrange(p_1+ggtitle(TeX(r'(N species = 5, \ $\alpha_e = 0.3)'))+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank()),
            p_2+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.25)'))+
              theme(axis.title.x = element_blank(),axis.text.x = element_blank(),axis.ticks.x=element_blank(),axis.line.x = element_blank(),
                    axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            p_3+ggtitle(TeX(r'(N species = 15, \ $\alpha_e = 0.325)')),
            p_5+ggtitle(TeX(r'(N species = 25, \ $\alpha_e = 0.25)'))+
              theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y=element_blank(),axis.line.y = element_blank()),
            nrow=2,ncol=2,heights = c(1,1.1),common.legend = T,legend="bottom")

ggsave("../Figures/Final_figs/SI/Example_bifu_Nspecies.pdf",width = 7,height = 5)



## Which species shifts and how frequently ? ----

Trait_shift=read.table("../Table/N_species/MF/Trait_shift.csv",sep=";")
N_replicate=150


Trait_shift_sum=Trait_shift%>%
  group_by(., Nsp,Branch,Competition,Species)%>%
  summarise(., .groups = "keep",Sum_tipping=sum(Tipping))

Trait_shift_sum$Trait=sapply(1:nrow(Trait_shift_sum),function(x){
  return(seq(1,0, length.out=Trait_shift_sum$Nsp[x])[Trait_shift_sum$Species[x]])
})

p=ggplot(filter(Trait_shift_sum), 
         aes(x=Trait,y=Sum_tipping/N_replicate,fill=Trait)) + 
  geom_bar( stat="identity",width = .05)+
  facet_grid(Nsp~Competition)+
  the_theme+
  labs(x=TeX("$\\psi$"),y="Frequency of abrupt shift",fill="")+
  scale_fill_gradientn(colours = color_Nsp(100))

ggsave("../Figures/Final_figs/SI/Traits_which_shifts.pdf",p,width = 10,height = 6)



## All CSI/Rho_+ bifurcation competition/Nsp ----

# Community index

d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

p=ggplot(d_tot%>%filter(., Branch==1))+
  geom_point(aes(x=Stress,y=CSI,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=13),strip.text.y = element_text(size=12),
        panel.background = element_blank(),strip.background.x = element_blank())

ggsave("../Figures/Final_figs/SI/CSI_aij_Nsp.pdf",p,width = 12,height = 7)


#community scale vegetation cover

d_tot=read.table("../Table/N_species/MF/Multistability_CSI.csv",sep=";")

p=ggplot(d_tot%>%filter(., Branch==1))+
  geom_point(aes(x=Stress,y=Rho_plus,color=Psi_normalized,fill=Psi_normalized),size=.5,shape=21)+
  the_theme+labs(y="Community index",color="")+
  scale_color_gradientn(colors = color_Nsp(100),na.value = "black")+
  scale_fill_gradientn(colors = color_Nsp(100),na.value = "black")+
  guides(fill="none")+
  labs(x="Stress, S")+
  facet_grid(Nsp~Competition,labeller = label_bquote(cols= alpha[e] ==.(Competition),rows="N species"==.(Nsp)))+
  theme(strip.text.x = element_text(size=13),strip.text.y = element_text(size=12),
        panel.background = element_blank(),strip.background.x = element_blank())

ggsave("../Figures/Final_figs/SI/Vegetation_cover_aij_Nsp.pdf",p,width = 12,height = 7)







