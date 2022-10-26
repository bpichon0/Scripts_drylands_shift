source('./N_species_Analysis_functions.R')

# Step 1 : N species MF -----
## 0) First test with 15 species ----

d_15sp=as_tibble(read.table("../Table/N_species/MF/15_species.csv",sep=","))
colnames(d_15sp)=c(paste0("Sp_",1:15),"Fertile","Degraded","Stress","alpha_0","Community_compo","Ini","Branch")  
d_15sp$Rho_p=rowSums(d_15sp[,1:15])
d_15sp$CSI = sapply(1:nrow(d_15sp),function(x){
  set.seed(123)
  u=runif(15)
  return(sum(d_15sp[x,1:15]*u))
  }
)



p=ggplot(melt(d_15sp,measure.vars = c("Rho_p")))+
  geom_line(aes(x=Stress,value,group=interaction(variable,Branch,Ini)),size=.5)+
  the_theme+facet_wrap(.~Branch)+
  labs(x="Stress (S)",y="Densities",color="")

ggsave("../Figures/N_species/MF/15_species_scenario_community.pdf",p,width = 8,height = 5)


ggplot(melt(d_15sp,measure.vars = paste0("Sp_",1:15)))+
  geom_line(aes(x=Stress,value,group=interaction(variable,Branch,Ini)),size=.5)+
  the_theme+facet_wrap(.~Branch)+
  labs(x="Stress (S)",y="Densities",color="")

p=ggplot(melt(d_15sp,measure.vars = c(paste0("Sp_",1:15))))+
  geom_line(aes(x=Stress,value,color=variable,group=interaction(variable,Branch,Ini)),size=.5)+
  the_theme+facet_grid(Branch~alpha_0,scales = "free")+
  labs(x="Stress (S)",y="Densities",color="")+
  scale_color_manual(values=c(color_Nsp(15),"gray50"),labels=c(paste0("Sp",1:15),"All"))

ggsave("../Figures/N_species/MF/15_species_scenario_species.pdf",p,width = 8,height = 5)


p=ggplot(melt(d_15sp,measure.vars = "CSI"))+
  geom_point(aes(x=Stress,value),size=1,alpha=.4,shape=21)+
  the_theme+facet_grid(Branch~.,scales = "free")+
  labs(x="Stress (S)",y="Densities",color="")
ggsave("../Figures/N_species/MF/CSI_15_species.pdf",p,width = 8,height = 5)

## 1) Varying facilitation and competition ----


Nsp=15
d = read.table("../Table/N_species/MF/Varying_competition_facilitation.csv",sep=",")
colnames(d) = c(paste0("Rho_",1:Nsp),"Rho_0","Rho_d","Stress","alpha_0","Facilitation")
d$Rho_p=rowSums(d[,1:Nsp])
d$CSI = sapply(1:nrow(d),function(x){
  set.seed(123)
  u=runif(Nsp)
  return(sum(d[x,1:Nsp]*u))
}
)

ggplot(d%>%melt(., measure.vars=c("CSI")))+
  geom_point(aes(x=Stress,y=value))+
  the_theme+
  facet_grid(alpha_0~Facilitation)

ggplot(d%>%melt(., measure.vars=paste0("Rho_",1:Nsp)))+
  geom_line(aes(x=Stress,y=value,group=interaction(Facilitation,alpha_0,variable)))+
  the_theme+
  facet_grid(alpha_0~Facilitation)
