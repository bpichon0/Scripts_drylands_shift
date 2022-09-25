# Other) Recruitment rate ----
rm(list = ls())
library(tidyverse)
d = expand_grid(emax = 1, e = .5, S = seq(0, 1, length.out = 200), psi = c(0, .5, 1))

p = ggplot(d) +
    geom_line(aes(x = S, y = emax * (1 - S * (1 - e * psi)), color = as.factor(psi), group = psi)) +
    scale_color_manual(values = c("blue", "green", "red")) +
    theme_classic() +
    theme(legend.position = "bottom") +
    labs(
        x = "S",
        y = TeX("$\\epsilon_{max} (1-S (1-e \\psi))$"), color = TeX("\\psi")
    )



ggsave("../Figures/Recruitment_rate.pdf", p, width = 6, height = 4)



# Step 1) Two species MF model ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#



## 1) No difference of competition, influence of tolerance ----


param = c(r = 0.05, d = 0.1, f = 0.9, beta = 0.8, m = 0.1, e = 0.1, emax = 1.2, c = 0.3, S = 0)
state = c(rho_1 = 0.4, rho_2 = 0.4, rho_d = 0.1, rho_0 = 0.1)
tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)

MF_two_species_julia = julia_eval("
function MF_two_species(du,u,p,t)
  r,d,f,beta,m,e,emax,c,S=p
  rho_1,rho_2,rho_d,rho_0=u

  du[1] = rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c*(rho_2 + rho_1))) - rho_1 * m
  du[2] = rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c*(rho_2 + rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1))-rho_0 * ( beta * rho_1 *  (emax *( 1 -S*(1-e)) - c*(rho_2 + rho_1)))  -
      rho_0 * ( beta * rho_2 *  (emax *( 1 -S) - c*(rho_2 + rho_1))) + rho_2 * m + rho_1 * m


end")




d2 = tibble()
S_seq = seq(0, 1, length.out = 100)
for (S in S_seq) {
  for (e in c(.1, .4)) {
    param[9] = S
    param[6] = e
    julia_assign("p", param)
    prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
    sol = de$solve(prob, de$Tsit5(), saveat = t)
    d = as.data.frame(t(sapply(sol$u, identity)))
    d2 = rbind(d2, d[nrow(d), ] %>% add_column(e = e, S = S))
  }
}
colnames(d2) = c("rho_1", "rho_2", "rho_d", "rho_0", "e", "S")
d2$rho_plus = d2$rho_1 + d2$rho_2

p1 = ggplot(d2 %>% melt(., id.vars = c("S", "e"))) +
  geom_point(aes(x = S, y = value, color = as.factor(variable), group = variable), size = 1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Stress", y = "Density") +
  scale_color_manual(
    name = "", values = color_rho,
    labels = c(TeX("$\\rho_1$"), TeX("$\\rho_2$"), TeX("$\\rho_d$"), TeX("$\\rho_0$"), TeX("$\\rho_+$"))
  ) +
  facet_wrap(. ~ e, labeller = label_both) +
  theme(legend.text = element_text(size = 14))

ggsave("../Figures/2_species/MF/Bifu_kefi_2_species_no_compet.pdf", p1, width = 7, height = 4)
















## 2) State diagram with competitive ability and stress for different type of competition ----

# We now try to replicate Danet figure with abiotic stress in x-axis and differential competitive ability in y-axis


# Intra + interspecific competition


# preparing the 3 types of simulations


state =Get_MF_initial_state()
tspan = c(0, 4000)
t = seq(0, 4000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)


for_sim = tibble(
  scena = c("inter=0", "low_inter", "inter=intra", "intra=0"),
  cinter2 = c(0, .05, .1, .15),
  cintra = c(.25, .2, .1, 0)
)[-1,]


for (i in 1:nrow(for_sim)) {
  
  param = c(
    r = 0.05, d = 0.1, f = 0.9, beta = 0.8, m = 0.1, e = .1, emax = 1.2, cintra = for_sim$cintra[i],
    cinter1 = for_sim$cinter2[i], cinter2 = for_sim$cinter2[i], S = 0
  )
  
  # d2=read.table(paste0("../Table/2_species/2_species_",for_sim$scena[i],".csv"),sep=";")
  
  
  d2 = tibble()
  S_seq = seq(0, 1, length.out = 30)
  c_seq = seq(param["cinter1"], 4 * param["cinter1"], length.out = 30)
  for (S in S_seq) {
    
    for (c1 in c_seq) {
      
      param["S"] = S
      param["cinter1"] = c1
      julia_assign("p", param)
      
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
      
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = S, cinter1 = c1))
    }
  }
  d2[d2 < 10^-3] = 0
  colnames(d2) = c("rho_1", "rho_2", "rho_d", "rho_0", "S", "cinter1")
  d2$rho_plus = d2$rho_1 + d2$rho_2
  
  
  Post_processing_MF(d2, for_sim$scena[i], C_seq = c_seq)
}





## 3) Hysteresis size as a function of competition coefficient ----


Run_dynamics_hysteresis(plot = T, N_seq_c = 3, N_seq_S = 300)

# Analysis of hysteresis size
param = Get_MF_parameters()
Nc = 120
N_seq_S = 500
color_rho = c("coexistence" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")

# d = Run_dynamics_hysteresis(plot = F, N_seq_c = Nc, N_seq_S = N_seq_S, write_table = T)


d = read.table(paste0("../Table/2_species/Hysteresis_size_length_C_", Nc, "_length_S_", N_seq_S, ".csv"), sep = ";")
d_hysteresis = Compute_hysteresis(d) %>% add_column(Rel_comp = rep(seq(param["cinter2"], 4 * param["cinter1"], length.out = Nc), each = 2))

p = ggplot(d_hysteresis %>% mutate(Species = as.character(Species)) %>%
             mutate(., Species = recode_factor(Species, "1" = "stress_tol", "2" = "competitive"))) +
  geom_line(aes(x = Rel_comp / param["cinter2"], y = Hysteresis, color = Species, group = Species), lwd = 1) +
  the_theme +
  scale_color_manual(values = color_rho[c(2, 4)]) +
  labs(x = "", y = "Hysteresis size")

p2 = ggplot(d_hysteresis %>% mutate(Species = as.character(Species)) %>%
              mutate(., Species = recode_factor(Species, "1" = "stress_tol", "2" = "competitive")) %>%
              melt(., measure.vars = c("Tipping_D", "Tipping_R")) %>%
              mutate(., variable = recode_factor(variable, "Tipping_D" = "Degradation", "Tipping_R" = "Restoration"))) +
  geom_line(aes(x = Rel_comp / param["cinter2"], y = value, color = Species, group = Species), lwd = 1) +
  the_theme +
  scale_color_manual(values = color_rho[c(2, 4)]) +
  facet_wrap(. ~ variable) +
  labs(x = TeX(r'(Relative competitive ability \ $\frac{c_{2,1}}{c_{1,2}})'), y = "Tipping point stress (S)")


p_tot = ggarrange(p, p2, nrow = 2, labels = LETTERS[1:2], common.legend = T, legend = "bottom")

ggsave("../Figures/2_species/MF/Evolution_hysteresis_competition.pdf", p_tot, width = 7, height = 7)










## 4) Interplay between facilitation & competition  ----

# interplay between facilitation and competition on the coexistence and type of transition
# For that, we use 3 values of stress for which we do heat maps of the density of species.

state = c(rho_1 = 0.4, rho_2 = 0.4, rho_d = 0.1, rho_0 = 0.1)
tspan = c(0, 5000)
t = seq(0, 5000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)

param = c(r = 0.05, d = 0.1, f = 0.9, beta = 0.8, m = 0.1, e = .1, emax = 1.2, cintra = .1, cinter1 = .1, cinter2 = .1, S = 0)

c_seq = seq(param["cinter1"], 4 * param["cinter1"], length.out = 30)
f_seq = seq(0, 1.5, length.out = 30)
S_seq = c(0, .25, .5, .75)
d2 = tibble()

for (S in S_seq) {
  for (f in f_seq) {
    for (c1 in c_seq) {
      param[3] = f
      param[9] = c1
      param[11] = S
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
      d2 = rbind(d2, d[nrow(d), ] %>% add_column(facil = f, cinter1 = c1, S = S))
    }
  }
}
d2[d2 < 10^-3] = 0
d2$rho_plus = d2$rho_1 + d2$rho_2


density_col = colorRampPalette(c("red", "white", "blue"))

p = ggplot(d2) +
  geom_tile(aes(x = facil, y = cinter1 / param["cinter2"], fill = rho_plus)) +
  theme_classic() +
  labs(x = "Facilitation strength (f)", y = "Relative competition strength", fill = TeX("$\\rho_{+}$")) +
  scale_fill_gradientn(colours = density_col(100)) +
  facet_wrap(. ~ S)

ggsave("../Figures/2_species/MF/Interplay_facilitation_competition.pdf", width = 7, height = 6)










## 5) Net effects between species ----

# we want to calculate the net effect between both species. We evaluate that by increasing slightly the recruitment rate of 1 at eq

state = Get_MF_initial_state()
tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)


length_seq = 500;epsilon=10^(-5)
S_seq = seq(0, 1, length.out = length_seq)

for (func_MF in c("normal","SGH")){
  for (scena in c("low_inter")){
    
    d2 = d3 = tibble()
    
    if (scena == "low_inter"){
      param = c(
        r = 0.05, d = 0.1, f = 0.9, beta1 = 0.8, beta2 = .8, m = 0.1, e = .1, emax = 1.2, cintra = .2,
        cinter1 = .05, cinter2 = .05, S = 0
      )
    } else {
      param = c(
        r = 0.05, d = 0.1, f = 0.9, beta1 = 0.8, beta2 = .8, m = 0.1, e = .1, emax = 1.2, cintra = .1,
        cinter1 = .1, cinter2 = .1, S = 0
      )
    }
    
    c_seq = seq(param["cinter1"], 5 * param["cinter1"], length.out = length_seq)
    C_for_analyse = c_seq[c(1,round((3/4)*length(c_seq)),length(c_seq))]
    
    for (c1 in C_for_analyse) {
      for (S in S_seq) {
        for (sp in 1:2) { #for both species
          
          param["S"] = S
          param["cinter1"] = c1 # update parameters
          
          # first without press
          
          julia_assign("p", param)
          if (func_MF=="normal"){
            prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
          } else {
            prob = julia_eval("ODEProblem(MF_two_species_SGH_press, state, tspan, p)")
          }
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
          colnames(mean_densities)="Densities"
          d2 = rbind(d2, mean_densities %>% add_column(S = S, cinter1 = c1, Type = "control",Species=c(1,2)[sp])) 
          d3 = rbind(d3, as_tibble(d[nrow(d),]) %>% add_column(S = S, cinter1 = c1))
          
          # with press
          
          param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
          julia_assign("p", param)
          if (func_MF=="SGH"){
            prob = julia_eval("ODEProblem(MF_two_species_press, state, tspan, p)")
          } else {
            prob = julia_eval("ODEProblem(MF_two_species_SGH_press, state, tspan, p)")
          }
          sol = de$solve(prob, de$Tsit5(), saveat = t)
          d = as.data.frame(t(sapply(sol$u, identity)))
          mean_densities=as_tibble(t(get_mean_densities(d)))[,c(1,2)[-sp]]
          colnames(mean_densities)="Densities"
          d2 = rbind(d2, mean_densities %>% add_column(S = S, cinter1 = c1, Type = "press",Species=c(1,2)[sp]))
          
          param[paste0("beta",sp)] = .8
        }
      }
    }
    
    d2[d2 < epsilon] = 0
    colnames(d2) = c("eq", "S", "cinter1", "Type","Species")
    
    
    net_effect =sapply(1:length(seq(1, nrow(d2) , by = 2)),function(x){
      i=seq(1, nrow(d2) , by = 2)[x]
      return((d2$eq[i+1] - d2$eq[i]) / epsilon)
    })
    
    d_net=d2%>%
      filter(., Type=="control")%>%
      select(.,-Type)
    d_net$value=net_effect
    
    #for latex format in facet
    appender <- function(string) {
      TeX(paste("$\\frac{c_{0,1}}{c_{1,0}} = $", string))  
    }
    
    p=ggplot(d_net %>%
               mutate(., cinter1=as.factor(round(cinter1/param["cinter2"],2)))) +
      geom_line(aes(x = S, y = value, color = as.factor(Species)),lwd=1,alpha=.5) +
      the_theme+scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
      geom_hline(yintercept = 0,linetype=9)+
      facet_wrap(.~cinter1)+
      labs(x="Stress (S)",alpha=TeX('$c_{2,1} / c_{1,2}$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_i}}{\partial b_j}$)'))+
      theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank())
    
    
    
    colnames(d3)=c("rho_1","rho_2","rho_d","rho_0","S","cinter1")
    
    p2=ggplot(d3%>%select(., -rho_d,-rho_0)%>%
                melt(., id.vars=c("S","cinter1"))%>%
                mutate(., cinter1=as.factor(round(cinter1/param["cinter2"],2))))+
      
      geom_line(aes(x=S,y=value,color=as.factor(variable),alpha=cinter1),lwd=1,alpha=.5)+
      the_theme+scale_color_manual(values=c("blue","green"),labels=c("rho_1"= "Stess_tol","rho_2"="Competitive"))+
      facet_wrap(.~cinter1,labeller = as_labeller(appender, default = label_parsed))+
      labs(x="Stress (S)",alpha=TeX('$c_{2,1}/c_{1,2}$'),color="Species",y="Densities")
    
    
    p_tot=ggarrange(p2,p,common.legend = T,legend = "bottom",nrow=2,align = "v")
    ggsave(paste0("../Figures/2_species/MF/Net_effects_",scena,"_",func_MF,".pdf"),p_tot,width = 8,height = 5)
  }
}


## 6) Comparing recruitment rate along stress S ----



state =Get_MF_initial_state()
tspan = c(0, 10000)
t = seq(0, 10000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)


for_sim = tibble(
  scena = c("low_inter", "inter=intra", "intra=0"),
  cinter2 = c(.05, .1, .15),
  cintra = c(.2, .1, 0)
)

d2 = tibble()

for (i in 1:nrow(for_sim)) {
  
  param = c(
    r = 0.05, d = 0.1, f = 0.9, beta = 0.8, m = 0.1, e = .1, emax = 1.2, cintra = for_sim$cintra[i],
    cinter1 = for_sim$cinter2[i], cinter2 = for_sim$cinter2[i], S = 0
  )
  
  
  S_seq = seq(0, 1, length.out = 100)
  c_seq = seq(param["cinter1"], 4 * param["cinter1"], length.out = 3)
  for (S in S_seq) {
    
    for (c1 in c_seq) {
      
      param["S"] = S
      param["cinter1"] = c1
      julia_assign("p", param)
      
      prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
      recruit_rate=Get_recruitment_rates(d[nrow(d), ],param)
      d2 = rbind(d2, as_tibble(t(recruit_rate)) %>% add_column(Scenario=for_sim$scena[i],S = S, cinter1 = round(c1/param["cinter2"],4),rho_1=d[nrow(d), 1],rho_2=d[nrow(d), 2]))
    }
  }
  d2[d2 < 10^-3] = 0
  colnames(d2) = c("r1","r2","Scenario", "S", "cinter1","rho_1","rho_2")
  
}

p=ggplot(d2 %>% 
           melt(.,measure.vars=c("rho_1","rho_2","r1","r2"))%>%
           mutate(type_data=sapply(1:nrow(.),function(x){
             if (.$variable[x] %in% c("rho_1","rho_2")){
               return("Density")
             } else {
               return("Recruitment_rate")
             }
           }))%>%
           mutate(species=sapply(1:nrow(.),function(x){
             if (.$variable[x] %in% c("rho_1","r1")){
               return("Stress_tol")
             } else {
               return("Competitive")
             }
           }))%>%
           
           mutate(.,variable=recode_factor(variable,"rho_1"="stress_tol","rho_2"="competitive")))+
  geom_line(aes(x=S,y=value,color=variable),lwd=1)+labs(x="Stress (S)",y="Density/Growth rate",color="")+
  the_theme+  theme(legend.text = element_text(size=12))+
  facet_grid(Scenario~cinter1)+
  scale_color_manual(values=c(color_rho[c(2,4)],"r1"="#EA7474",'r2'="#8E63AB"))

ggsave("../Figures/2_species/MF/Relative_recruitment_rate.pdf",p,width = 7,height = 5)
#at eq drho/dt = 0 ...


# Step 2) Pair approximation (PA) ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#

## 1) Exploration along competitive ability ----




# As I get oscillations, function to get the mean densities
param=Get_PA_parameters()
state=Get_PA_initial_state()

tspan = c(0, 6000)
t = seq(0, 6000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("p", param)
julia_assign("tspan", tspan)


N_rep = 50
d2 = tibble()
S_seq = seq(0, 1, length.out = N_rep)
alpha_seq = seq(param["alpha12"], 4 * param["alpha12"], length.out = N_rep)
for (S in S_seq) {
  
  for (alpha21 in alpha_seq) {
    
    param["S"] = S
    param["alpha21"] = alpha21
    julia_assign("p", param)
    
    prob = julia_eval("ODEProblem(PA_two_species, state, tspan, p)")
    sol = de$solve(prob, de$Tsit5(), saveat = t)
    d = as.data.frame(t(sapply(sol$u, identity)))
    
    colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
    
    d2 = rbind(d2, as_tibble(t(get_mean_densities(d))) %>% add_column(S = S, alpha21 = alpha21))
  }
}
d2[d2 < 10^-4] = 0
colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha21")
d2$rho_plus = d2$rho_1 + d2$rho_2

# postprocessing
post_processing_2species_PA(d2, alpha_seq, S_seq)


#exploring clustering coefficients (e.g., q_{1|1})







## 2) Testing the influence of global and local competition ----

param=Get_PA_parameters()
state=Get_PA_initial_state()
tspan = c(0, 6000)
t = seq(0, 6000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)

scenarios_sim=tibble(name=c("high_global","only_local","local=global"),
                     cglobal=c(.2,0,.1),
                     alpha_intra=c(0.01,0.1,.05),
                     alpha_inter=c(0.01,.1,.05))

N_rep=50
S_seq=seq(0,1,length.out=N_rep)

for (nsim in 1:nrow(scenarios_sim)){
  
  d2 = tibble()
  param["cg"]=scenarios_sim$cglobal[nsim]
  param["alpha11"]=param["alpha22"]=scenarios_sim$alpha_intra[nsim]
  param["alpha21"]=param["alpha12"]=scenarios_sim$alpha_inter[nsim]
  alpha_seq = seq(param["alpha12"], 4 * param["alpha12"], length.out = N_rep)
  
  for (S in S_seq) {
    
    for (alpha21 in alpha_seq) {
      
      param["S"] = S
      param["alpha21"] = alpha21
      
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(PA_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
      
      d2 = rbind(d2, as_tibble(t(get_mean_densities(d))) %>% add_column(S = S, alpha21 = alpha21,scenario=scenarios_sim$name[nsim]))
    }
  }
  d2[d2 < 10^-4] = 0
  colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha21","scenario")
  
  assign(paste0("p_tot_",nsim),post_processing_2species_PA(d2, alpha_seq, S_seq)+ggtitle(scenarios_sim$name[nsim]))
}


p_final=ggarrange(p_tot_1,p_tot_2,p_tot_3,nrow=3,common.legend = T,legend = "bottom")
ggsave("../Figures/2_species/PA/PA_scenarios_global_local.pdf",p_final,width = 14,height = 18)

## 3) Influence of the strength of competition (local + global) ----


param=Get_PA_parameters()
state=Get_PA_initial_state()
tspan = c(0, 4000)
t = seq(0, 4000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)

scenarios_sim=tibble(name=c("high","medium","normal"),
                     cglobal=c(.3,0.2,.1),
                     alpha_intra=c(0.3,0.2,.1),
                     alpha_inter=c(0.3,.2,.1))

N_rep=50
S_seq=seq(0,1,length.out=N_rep)

for (nsim in 1:nrow(scenarios_sim)){
  d2 = tibble()
  param["cg"]=scenarios_sim$cglobal[nsim]
  param["alpha11"]=param["alpha22"]=scenarios_sim$alpha_intra[nsim]
  param["alpha21"]=param["alpha12"]=scenarios_sim$alpha_inter[nsim]
  alpha_seq = seq(param["alpha12"], 4 * param["alpha12"], length.out = N_rep)
  
  for (S in S_seq) {
    
    for (alpha21 in alpha_seq) {
      
      param["S"] = S
      param["alpha21"] = alpha21
      
      julia_assign("p", param)
      prob = julia_eval("ODEProblem(PA_two_species, state, tspan, p)")
      sol = de$solve(prob, de$Tsit5(), saveat = t)
      d = as.data.frame(t(sapply(sol$u, identity)))
      colnames(d) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
      
      d2 = rbind(d2, as_tibble(t(get_mean_densities(d))) %>% add_column(S = S, alpha21 = alpha21,scenario=scenarios_sim$name[nsim]))
    }
  }
  d2[d2 < 10^-6] = 0
  colnames(d2) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm", "S", "alpha21","scenario")
  
  assign(paste0("p_tot_",nsim),post_processing_2species_PA(d2, alpha_seq, S_seq)+ggtitle(scenarios_sim$name[nsim]))
}

p_final=ggarrange(p_tot_1,p_tot_2,p_tot_3,nrow=3,common.legend = T,legend = "bottom")
ggsave("../Figures/2_species/PA/PA_strength_competition.pdf",p_final,width = 14,height = 18)


## 4) Net effects of species ----

param=Get_PA_parameters()
state=Get_PA_initial_state()
tspan = c(0, 8000)
t = seq(0, 8000, by = 1)
julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)


length_seq = 200;epsilon=10^(-4)
S_seq = seq(0, 1, length.out = length_seq)

for (scena in c("low_inter","high_inter")[-2]){ 
  
  d2 = d3 = tibble()
  
  if (scena == "low_inter"){
    param = c(    r = 0.05, d = 0.1, f = .9, beta1 = 0.8,beta2 = 0.8, m = 0.1, e = .9, emax = 1.2, cg = 0.05, alpha11 = .2,
                  alpha12 = .05, alpha21 = .05, alpha22 = .2, S = 0, delta = .1, z = 4)
  } else {
    param = c(    r = 0.05, d = 0.1, f = .9, beta1 = 0.8,beta2 = 0.8, m = 0.1, e = .9, emax = 1.2, cg = .05, alpha11 = .1,
                  alpha12 = .1, alpha21 = .1, alpha22 = .1, S = 0, delta = .1, z = 4)
  }
  
  alpha_seq = seq(param["alpha12"], 10 * param["alpha12"], length.out = length_seq)
  alpha_for_analyse = alpha_seq[c(1,round(length(alpha_seq)/2),length(alpha_seq))]
  
  for (alpha21 in alpha_for_analyse) {
    for (S in S_seq) {
      for (sp in 1:2) { #for both species
        
        param["S"] = S
        param["alpha21"] = alpha21 # update parameters
        
        # first without press
        
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(PA_two_species_press, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        mean_densities=as_tibble(t(get_mean_densities(d)))[,c(c(1,2)[-sp],4,7,8)] # we get the species for which we didn't change recruitment rate and the pair densities 
        colnames(mean_densities)="Densities"
        d2 = rbind(d2, mean_densities %>% add_column(S = S, alpha21 = alpha21, Type = "control",Species=c(1,2)[sp])) 
        d3 = rbind(d3, as_tibble(d[nrow(d),c(1,2,4,7,8)]) %>% add_column(S = S, alpha21 = alpha21))
        
        # with press
        
        param[paste0("beta",sp)] = param[paste0("beta",sp)] + epsilon # press
        julia_assign("p", param)
        prob = julia_eval("ODEProblem(PA_two_species_press, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        mean_densities=as_tibble(t(get_mean_densities(d)))[,c(c(1,2)[-sp],4,7,8)]
        colnames(mean_densities)="Densities"
        d2 = rbind(d2, mean_densities %>% add_column(S = S, alpha21 = alpha21, Type = "press",Species=c(1,2)[sp]))
        
        param[paste0("beta",sp)] = .8
      }
    }
  }
  
  d2[d2 < epsilon] = 0
  d3[d3 < epsilon] = 0
  colnames(d2) = c("eq","rho_12","rho_11",'rho_22', "S", "alpha21", "Type","Species")
  
  d_net=d2%>%
    filter(., Type=="control")%>%
    select(.,-Type)
  
  for (cols_rho in 1:4){ #evaluating net effect for each global density or each pair
    d_net[,cols_rho]=unlist(sapply(1:length(seq(1, nrow(d2) , by = 2)),function(x,y=cols_rho){
      i=seq(1, nrow(d2) , by = 2)[x]
      return((d2[i+1,y] - d2[i,y]) / epsilon)
    }))}
  
  
  #for latex format in facet
  appender <- function(string) {
    TeX(paste("$\\frac{c_{0,1}}{c_{1,0}} = $", string))  
  }
  
  #global densities plots
  p=ggplot(d_net %>%
             mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2)))) +
    geom_line(aes(x = S, y = eq, color = as.factor(Species)),lwd=1,alpha=.5) +
    geom_hline(yintercept = 0,linetype=9)+
    scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
    facet_wrap(.~alpha21)+
    labs(x="Stress (S)",alpha=TeX('$c_{2,1} / c_{1,2}$'),color="Species",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_i}{\partial b_j}$)'))+
    the_theme+theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank())
  
  
  
  colnames(d3)=c("rho_1","rho_2","rho_12","rho_11",'rho_22',"S","alpha21")
  
  pbis=ggplot(d3%>%
                melt(., measure.vars=c("rho_1","rho_2"))%>%
                mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2))))+
    
    geom_line(aes(x=S,y=value,color=as.factor(variable),alpha=alpha21),lwd=1,alpha=.5)+
    the_theme+scale_color_manual(values=c("blue","green"),labels=c("rho_1"= "Stess_tol","rho_2"="Competitive"))+
    facet_wrap(.~alpha21,labeller = as_labeller(appender, default = label_parsed))+
    labs(x="Stress (S)",alpha=TeX('$c_{2,1}/c_{1,2}$'),color="Species",y="Densities")
  
  
  Pglob=ggarrange(pbis,p,common.legend = T,legend = "bottom",nrow=2,align = "v",labels = letters[1:2])
  ggsave(paste0("../Figures/2_species/PA/Net_effects_PA_",scena,"_global.pdf"),Pglob,width = 8,height = 5)
  
  
  #local densities plots
  
  assign(paste0("plocal_",11),ggplot(d_net %>%
                                       mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2)))) +
           geom_line(aes(x = S, y = rho_11, color = as.factor(Species)),lwd=1,alpha=.5) +
           geom_hline(yintercept = 0,linetype=9)+
           scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
           facet_wrap(.~alpha21)+
           labs(x="Stress (S)",alpha=TeX('$c_{2,1} / c_{1,2}$'),color="j",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_1,\psi_1}}{\partial b_j}$)'))+
           the_theme+theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank()))
  
  assign(paste0("plocal_",22),ggplot(d_net %>%
                                       mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2)))) +
           geom_line(aes(x = S, y = rho_22, color = as.factor(Species)),lwd=1,alpha=.5) +
           geom_hline(yintercept = 0,linetype=9)+
           scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
           facet_wrap(.~alpha21)+
           labs(x="Stress (S)",alpha=TeX('$c_{2,1} / c_{1,2}$'),color="j",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_0,\psi_0}}{\partial b_j}$)'))+
           the_theme+theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank()))
  
  assign(paste0("plocal_",12),ggplot(d_net %>%
                                       mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2)))) +
           geom_line(aes(x = S, y = rho_12, color = as.factor(Species)),lwd=1,alpha=.5) +
           geom_hline(yintercept = 0,linetype=9)+
           scale_color_manual(values=c("blue","green"),labels=c("1"= "Stess_tol","2"="Competitive"))+
           facet_wrap(.~alpha21)+
           labs(x="Stress (S)",alpha=TeX('$c_{2,1} / c_{1,2}$'),color="j",y=TeX(r'(Net effect \ \ $\frac{\partial \rho_{\psi_1,\psi_0}}{\partial b_j}$)'))+
           the_theme+theme(panel.grid = element_blank())+ theme(strip.text.x = element_blank(),strip.background = element_blank()))
  
  
  p2=ggplot(d3%>%
              melt(., measure.vars=c("rho_11","rho_22","rho_12"))%>%
              mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2))))+
    
    geom_line(aes(x=S,y=value,color=as.factor(variable),alpha=alpha21),lwd=1,alpha=.8)+
    the_theme+
    facet_wrap(.~alpha21,labeller = as_labeller(appender, default = label_parsed))+
    labs(x="Stress (S)",alpha=TeX('$c_{2,1}/c_{1,2}$'),color="Species",y="Densities")
  
  
  Ploc=ggarrange(p2,ggarrange(plocal_11+xlab(""),plocal_12+xlab(""),plocal_22,nrow=3,align="v",common.legend = T,legend = "bottom"),
                 nrow=2,labels = LETTERS[1:2],heights = c(1.5,3),align = "v")
  ggsave(paste0("../Figures/2_species/PA/Net_effects_PA_",scena,"_local.pdf"),Ploc,width = 8,height = 9)
  
  
  #local conditional probabilities
  pq=ggplot(d3%>%
              mutate(., q12=rho_12/rho_2,q21=rho_12/rho_1,q11=rho_11/rho_1,q22=rho_22/rho_2)%>%
              melt(., measure.vars=c("q12","q21",'q11',"q22","rho_1","rho_2"))%>%
              mutate(., alpha21=as.factor(round(alpha21/param["alpha12"],2))))+
    
    geom_line(aes(x=S,y=value,color=variable),lwd=1,alpha=.8)+
    
    the_theme+scale_color_manual(
      values=c(hue_pal()(4),"blue","green"),
      labels=c("q12"=TeX("$q_{\\psi_1|\\psi_0}$"),"q11"=TeX("$q_{\\psi_1|\\psi_1}$"),
               "q22"=TeX("$q_{\\psi_0|\\psi_0}$"),"q21"=TeX("$q_{\\psi_0|\\psi_1}$"),
               "rho_1"=TeX("$\\rho_{\\psi_1}$"),"rho_2"=TeX("$\\rho_{\\psi_0}$")))+
    
    facet_wrap(.~alpha21,labeller = as_labeller(appender, default = label_parsed))+
    labs(x="Stress (S)",alpha=TeX('$c_{2,1}/c_{1,2}$'),color="",y="Densities")  
  ggsave("../Figures/2_species/PA/PA_qij.pdf",pq,width = 7,height = 4)
  
}



# Step 3) CA between both species : example and clustering ----
rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()
#

## 0) Comparing Gillespie simulations with classical IBM  ----
rm(list = ls())
source("./2_Species_analysis_functions.R")


t_seq = seq(1, 2000, 1)
param = c(
  r = 0.02, d = 0.1, f = .9, beta = 0.8, m = 0.05, e = 0, emax = 1,
  cintra = 0.1, cinter1 = .1, cinter2 = .1, S = 0, delta = .1, z = 4, leap = 1
)
Lattice_ini = Get_initial_lattice()

plot_dynamics(Run_CA_2_species(t_seq, param, ini = Lattice_ini)$d)
plot_dynamics(Gillespie_tau_leaping_R(Lattice_ini, t_seq, param)$state)


t_compare = tibble()
for (k in 1:10) {
  t1 = system.time(Run_CA_2_species(t_seq, param, ini = Lattice_ini))[1]
  t2 = system.time(Gillespie_tau_leaping_R(Lattice_ini, t_seq, param))[1]
  t_compare = rbind(t_compare, tibble(Comp_time = c(t1, t2), Method = c("Classic", "Gillespie")))
}
t_compare$Comp_time = as.numeric(t_compare$Comp_time)
t_compare %>%
  group_by(Method) %>%
  summarise(m_t = mean(Comp_time), s_t = sd(Comp_time))



d = tibble()
for (tau in c(1, 0.1, .01, .001)) {
  param["leap"] = tau
  d = rbind(d, (Gillespie_tau_leaping_R(Lattice_ini, t_seq, param)$state) %>% add_column(., tau_leap = tau))
}
color_rho = c("fertile" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")
d$time = rep(t_seq, 4)
p = ggplot(d %>% melt(., id.vars = c("time", "tau_leap")) %>%
             mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive", "rho_0" = "fertile", "rho_d" = "desert"))) +
  geom_line(aes(x = time, y = value, color = variable), lwd = 1) +
  theme_classic() +
  scale_color_manual(values = color_rho) +
  labs(x = "Time", y = "Densities", color = "") +
  theme(legend.text = element_text(size = 11), legend.position = "bottom") +
  facet_wrap(. ~ tau_leap, labeller = label_both)

ggsave("../Figures/2_species/Tau_values_dynamics_Gillespie.pdf", width = 7, height = 4)

## 1) An example along interspecific competition gradient ----
t_seq = seq(1, 3000, 1)
param = Get_CA_parameters()
Lattice_ini = Get_initial_lattice(size = 100)
test = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
plot_dynamics(test$state)
test2 = Gillespie_tau_leaping_R(time = t_seq, param = param, landscape = Lattice_ini)
plot_dynamics(test2$state)

d = tibble()
for (s in c(0, .25, .5, .75, 1)) {
  param["S"] = s
  Lattice_ini = Get_initial_lattice()
  CA = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
  d = rbind(d, as_tibble(melt(CA$landscape)) %>% add_column(., S = s))
}

color_CA = c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
p2 = ggplot(d) +
  geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
  theme_transparent() +
  scale_fill_manual(values = color_CA, labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
  facet_grid(. ~ S, labeller = label_both) +
  theme(panel.border = element_blank()) +
  theme(legend.position = "bottom") +
  labs(fill = "")


p_tot = ggarrange(p1, p2, ncol = 2, widths = c(1, 4))

ggsave("../Figures/2_species/CA/CA_example.pdf", p_tot, width = 16, height = 4, device = "pdf")



# We may be able to have coexistence at S=0 with no difference in competitive ability by changing delta
param = c(r = 0.02, d = 0.1, f = .9, beta = 0.8, m = 0.05, e = .1, emax = 1, cintra = 0.1, cinter1 = .1, cinter2 = .1, S = 0, delta = .1, z = 4, tau_leap = 1)

d = tibble()
for (delta in c(0.1, .25, .5, .75, 1)) {
  param["delta"] = delta
  Lattice_ini = Get_initial_lattice()
  CA = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
  d = rbind(d, as_tibble(melt(CA$landscape)) %>% add_column(., delta = delta))
}
p = ggplot(d) +
  geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
  theme_transparent() +
  scale_fill_manual(values = color_CA, labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
  facet_grid(. ~ delta, labeller = label_both) +
  theme(panel.border = element_blank()) +
  theme(legend.position = "bottom") +
  labs(fill = "")

ggsave("../Figures/2_species/CA/Varying_delta_CA.pdf", p, width = 13, height = 4, device = "pdf")


## 2) Clustering metrics along inter-specific competition gradient ----

# we use the Julia Code 2_species_clustering.jl
param = c(r = 0.02, d = 0.1, f = .9, beta = 0.8, m = 0.05, e = .1, emax = 1, cintra = 0.1, cinter1 = .3, cinter2 = .1, S = 0, delta = .1, dt = 1, z = 4)

d = read.table("../Table/2_species/Clustering_species_interspe_compet_gradient.csv", sep = ",")
d[is.na(d)] = 0 # for sake of simplicity
colnames(d) = c("rep", "q12", "c12", "cpp", "Relative_compet", "S")

# q12
p1 = d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    q12_m = median(q12), q12_q1 = quantile(q12, 0.25), q12_q3 = quantile(q12, 0.75), .groups = "keep"
  ) %>%
  ggplot(.) +
  geom_line(aes(x = Relative_compet / param["cinter2"], y = q12_m), color = "#F1B943", lwd = 1) +
  geom_ribbon(aes(x = Relative_compet / param["cinter2"], ymin = q12_q1, ymax = q12_q3), fill = "#F1B943", alpha = .5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ S, labeller = label_both) +
  labs(x = "", y = TeX("$\\q_{1|2}$"))

# c12
p2 = d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    c12_m = median(c12), c12_q1 = quantile(c12, 0.25), c12_q3 = quantile(c12, 0.75), .groups = "keep"
  ) %>%
  ggplot(.) +
  geom_line(aes(x = Relative_compet / param["cinter2"], y = c12_m), color = "#119680") +
  geom_ribbon(aes(x = Relative_compet / param["cinter2"], ymin = c12_q1, ymax = c12_q3), fill = "#119680", alpha = .5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ S) +
  labs(x = "", y = TeX("$\\c_{12}$")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )


# cpp
p3 = d %>%
  group_by(S, Relative_compet) %>%
  summarise(
    cpp_m = median(cpp), cpp_q1 = quantile(cpp, 0.25), cpp_q3 = quantile(cpp, 0.75), .groups = "keep"
  ) %>%
  ggplot(.) +
  geom_line(aes(x = Relative_compet / param["cinter2"], y = cpp_m), color = "#E63535") +
  geom_ribbon(aes(x = Relative_compet / param["cinter2"], ymin = cpp_q1, ymax = cpp_q3), fill = "#E63535", alpha = .5) +
  theme_classic() +
  theme(legend.position = "bottom") +
  facet_wrap(. ~ S) +
  labs(x = "Relative competition ability", y = TeX("$\\c_{++}$")) +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  )



p_tot = ggarrange(p1, p2, p3, nrow = 3)
ggsave("../Figures/2_species/CA/Clustering_figure_interspe.pdf", p_tot, width = 8, height = 6)




## 3) Same but along stress gradient  ----

d = read.table("../Table/2_species/Clustering_species_stress_gradient.csv", sep = ",")
colnames(d) = c("rep", "q12",  "q21",  "q11" , "q22" , "c12",  "c21",  "c11",  "c22",  "qpp" , "cpp" , "stress")


d %>%
  melt(., measure.vars=c("q12",'q21',"c21"))%>%
  ggplot(.) +
  geom_line(aes(x = stress, y = value,color=variable),lwd=1) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(x = "Stress (S)", y = "")


# # cpp
# p4 = d %>%
#     group_by(stress) %>%
#     summarise(
#         cpp_m = median(cpp), cpp_q1 = quantile(cpp, 0.25), cpp_q3 = quantile(cpp, 0.75), .groups = "keep"
#     ) %>%
#     ggplot(.) +
#     geom_line(aes(x = stress, y = cpp_m), color = "#E63535") +
#     geom_ribbon(aes(x = stress, ymin = cpp_q1, ymax = cpp_q3), fill = "#E63535", alpha = .5) +
#     theme_classic() +
#     theme(legend.position = "bottom") +
#     labs(x = "Relative competition ability", y = TeX("$\\c_{++}$"))


p_tot = ggarrange(p1, p2, p3, nrow = 3)
ggsave("../Figures/2_species/CA/Clustering_figure_stress_gradient.pdf", width = 5, height = 5)




## 4) Testing CA on specific parameters values ----

t_seq = seq(1, 5000, 1)
param = Get_CA_parameters()
param["alpha11"] = 0.4
param["alpha12"] = 0.01
param["alpha21"] = 0.1
param["alpha22"] = 0.01
param["cg"] = 0.01
param["S"] = 0.5
param["tau_leap"] = 1
Lattice_ini = Get_initial_lattice()

test = Run_CA_2_species(time = t_seq, param = param, landscape = Lattice_ini)
plot_dynamics(test$state)

test2 = Gillespie_tau_leaping_R(time = t_seq, param = param, landscape = Lattice_ini)

plot_dynamics(test2$state)
Plot_landscape(test2$landscape)

#
## 5) Patch size distribution ----

### a) Example  ----

mapping(100,100)
S_seq=seq(0,1,length.out=6)
N_rep = 3 

d_freq=d_size=tibble()
display_landscape=tibble()
for (stress in S_seq){
  for (rep in 1:N_rep){
    print(stress)
    if (stress == 0){
      stress="0.0"
    }
    if (stress == 1){
      stress="1.0"
    }
    
    landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Example_sim/Landscape_size_100_stress_",stress,"_",rep,".csv"),sep=","))
    
    if (rep==1){
      display_landscape = rbind(display_landscape, as_tibble(melt(landscape)) %>% add_column(., S = stress))
      
    }
    
    Freq_patches=Get_frequency_number_patches(landscape)
    d_freq=rbind(d_freq,Freq_patches$Patches_frequency%>%add_column(., S=stress,Nrep=rep))
    d_size=rbind(d_size,Freq_patches$Patches_size%>%add_column(., S=stress,Nrep=rep))
  }
}



p1=ggplot(d_freq)+geom_point(aes(x=log(Size),log(Frequency),color=Species))+
  the_theme+labs(x="Patch size (k)",y="Frequency patch > k")+
  facet_wrap(.~S)+scale_color_manual(values=c("gray40",as.character(color_rho[c(4,2)])))

p2=ggplot(display_landscape) +
  geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
  theme_transparent() +
  scale_fill_manual(values =  c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
                    , labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
  theme(panel.border = element_blank()) +
  theme(legend.position = "bottom") +
  labs(fill = "")+
  facet_wrap(.~S)

p_tot=ggarrange(p1,p2,nrow=2,labels=LETTERS[1:2])
ggsave("../Figures/2_species/CA/Example_PSD.pdf",p_tot,width = 7,height = 10)



### b) Analysing the simulations for different competition regimes ----


# the simulations were made on Julia (Step 4 : Patch size distribution : different competition scenario)

mapping(100,100)
S_seq=seq(0,1,length.out=6)
N_rep = 3 

for (scena in c("low_inter","intra=inter")){
  for (rela_comp in c("1.0", "2.5", "4.0")){
  
    d_freq=d_size=d_max_size=tibble()
    display_landscape=tibble()
  
    for (stress in S_seq){
      
      for (rep in 1:N_rep){
        
        print(stress)
        
        if (stress == 0){
          stress="0.0"
        }
        if (stress == 1){
          stress="1.0"
        }
        
        landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Competition_regime/Landscape_size_100_stress_",stress,"_",scena,"_",rep,"_relat_compet_",rela_comp,".csv"),sep=","))
        if (rep==1 & length(unique(as.numeric(landscape)))>2){ #i.e. we have vegetation
          display_landscape = rbind(display_landscape, as_tibble(melt(landscape)) %>% add_column(., S = stress))
          
        }
        
        Freq_patches=Get_frequency_number_patches(landscape)
        d_freq=rbind(d_freq,Freq_patches$Patches_frequency%>%add_column(., S=stress,Nrep=rep))
        d_size=rbind(d_size,Freq_patches$Patches_size%>%add_column(., S=stress,Nrep=rep))
        d_max_size=rbind(d_max_size,d_size%>%group_by(Species,S)%>%summarise(max_patch=max(patch_size),.groups = "keep")) #max across replicates
      }
    }
  
    p1=ggplot(d_freq)+geom_point(aes(x=log(Size),log(Frequency),color=Species))+
      the_theme+labs(x="Patch size (k)",y="Frequency patch > k")+
      facet_wrap(.~S,labeller = label_both)+scale_color_manual(values=c("gray40",as.character(color_rho[c(4,2)])))
    
    p2=ggplot(display_landscape) +
      geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
      theme_transparent() +
      scale_fill_manual(values =  c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
                        , labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
      theme(panel.border = element_blank()) +
      theme(legend.position = "bottom") +
      labs(fill = "")+
      facet_grid(.~S,labeller = label_both)
    
    p_tot=ggarrange(p1,p2,nrow=2,labels=LETTERS[1:2],heights = c(1,.5))
    ggsave(paste0("../Figures/2_species/CA/Patch_size_distribution_",scena,"_Relative_comp_",rela_comp,".pdf"),p_tot,width = 7,height = 8)
  }
}

landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Landscape_test.csv"),sep=","))
Plot_landscape(landscape)



#



### c) Fitting Power law ----
mapping(100,100)
S_seq = seq(0, 0.8, length.out=30)
n_save = 25
relativ_comp = c("1.0", "2.5", "4.0")

#doing the loop on simulations
for (scena in c(1)[-2]){
  
  for (stress in S_seq){
    d_PL=tibble()
    
    stress=round(stress, 2)
    
    if (stress == 0){
      stress="0.0"
    }
    if (stress == 1){
      stress="1.0"
    }
    
    for (rela_c in relativ_comp){
      
      for (n in 1:n_save){
        
        
        #get landscape
        landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Big_sim/Landscape_l_S_",stress,
                                    "_a21_", rela_c,"_nsave_", n , ".csv"),sep=","))
        

        
        #psd = patch size distribution
        Freq_patches=Get_frequency_number_patches(landscape)
        d_psd=Freq_patches$Patches_frequency
        colnames(d_psd) =  c("n","p","size","Species")
  
        #perform analysis for each species or for global vegetation
        for (sp in c("+","1","2")){
  
          best_model= NA
          psd_sp=filter(d_psd,Species==sp)
  
          # 1: check if desert for the focal species

          if(sum(landscape %in% c(1,2))/length(landscape) <0.01& sp=="+" |
             sum(landscape %in% c(1))/length(landscape) <0.01   & sp=="1"  |
             sum(landscape %in% c(2))/length(landscape) <0.01 & sp=="2" ){ best_model = 1 }
  
          # 2: check if vegetated for all vegetation
          if(is.na(best_model) & sum(landscape %in% c(2))/length(landscape) >= 0.70  & sp=="2" |
             is.na(best_model) & sum(landscape %in% c(1,2))/length(landscape) >= 0.70  & sp=="+" |
             is.na(best_model) & sum(landscape %in% c(1))/length(landscape) >= 0.70  & sp=="1")  {
            best_model = 5
          }
  
          # 3: fit power law models and compare via AIC via function
  
          if(is.na(best_model)) {
  
  
            p_spanning <- tail(psd_sp$p,1)
  
            result <- fitPL(psd_sp, p_spanning)
  
            best_model = result$best
          }
  
          class = c("DEG", "DOWN","PL", "UP", "COV")[best_model]
  
          if (class %in% c("DEG","COV")){alpha=NA}
          if (class %in% c("UP")){alpha=coefficients(result$TPLup)["alpha"]}
          if (class %in% c("DOWN")){alpha=coefficients(result$TPLdown)["alpha"]}
          if (class %in% c("PL")){alpha=coefficients(result$PL)["alpha"]}
  
          d_PL=rbind(d_PL,tibble(Class=class,Max_patch=max(psd_sp$size),Alpha_expo=alpha,N_rep=n,
                                 Species=sp,a21=rela_c,Stress=stress))
  
        } #end loop species
      
      } #end loop replicates
      
    } #end loop a12
    
    write.table(d_PL,paste0("../Table/2_species/PL_summary/PL_",scena,"_S_",stress,".csv"),sep=";")
    
  } # end loop stress

} #end loop scenarios competition


#merging dataframes
d_PL=tibble()
for ( i in round(seq(0,.8,length.out=30),2)){if (i==0) {i = "0.0"}
  d_PL=rbind(d_PL,read.table(paste0("../Table/2_species/PL_summary/PL_1_S_",i,".csv"),sep=";"))}

d_PL$Max_patch[is.infinite(d_PL$Max_patch)]=0

d_PL$Class[which(d_PL$Species=="2" & d_PL$a21==4 & d_PL$Class=="PL")]="DOWN"

#for graphical purposes
appender <- function(string) {
  TeX(paste("$\\frac{\\alpha_{0,1}}{\\alpha_{1,0}} = $", string))}


color_sp=c("+"="gray70","1"=as.character(color_rho[4]),"2"= as.character(color_rho[2]))
color_class_PL =  c("COV"="#68B15E", "UP"= "#ECD57A",   "DOWN" ="#D22D2D","PL"="#E88D35",   "DEG"="#ADA9A9")

#Max patch size colored by species
d_PL_maxsize_sp=d_PL%>%
  group_by(., Species,a21,Stress)%>%
  summarise(.groups = "keep",Max_patch=mean(Max_patch))


p=ggplot(NULL)+
  geom_point(data=d_PL,aes(x=Stress,Max_patch,color=Species),size=.5,alpha=.3)+
  geom_line(data=d_PL_maxsize_sp,aes(x=Stress,Max_patch,color=Species),lwd=.9)+
  the_theme+labs(x="Stress (S)", y="Max patch size",color="")+
  facet_wrap(.~a21,labeller = as_labeller(appender, default = label_parsed))+
  scale_color_manual(values=color_sp,
                     labels=c("+"="All vegetation","1"="Stress_tol","2"= "Competitive"))+
  scale_y_log10()

ggsave("../Figures/2_species/CA/Max_patch_size_by_species.pdf",p,width = 7,height = 4)


#Max patch size colored by PL class
d_PL2=transform(d_PL,
                a21 = factor(a21, levels=c(1,2.5,4), labels=c("frac(alpha[0,1],alpha[1,0]) : 1", "frac(alpha[0,1],alpha[1,0]) : 2.5",
                                                              "frac(alpha[0,1],alpha[1,0]) : 4")),
                Species=factor(Species,levels=c("+","1","2"),labels=c("psi[0] + psi[1]","psi[1]","psi[0]")))


d_PL_maxsize_class=d_PL%>%
  group_by(., Species,a21,Stress,Class)%>%
  summarise(.groups = "keep",Max_patch=mean(Max_patch))

d_PL_maxsize_class=transform(d_PL_maxsize_class,
                a21 = factor(a21, levels=c(1,2.5,4), labels=c("frac(alpha[0,1],alpha[1,0]) : 1", "frac(alpha[0,1],alpha[1,0]) : 2.5",
                                                              "frac(alpha[0,1],alpha[1,0]) : 4")),
                Species=factor(Species,levels=c("+","1","2"),labels=c("psi[0] + psi[1]","psi[1]","psi[0]")))

p2=ggplot(NULL)+
  geom_point(data=d_PL2,aes(x=Stress,Max_patch,color=Class),size=.5,alpha=.3)+
  geom_line(data=d_PL_maxsize_class,aes(x=Stress,Max_patch,color=Class),lwd=.9)+
  the_theme+labs(x="Stress (S)", y="Max patch size",color="")+
  facet_grid(Species~a21,labeller=label_parsed)+
  scale_color_manual(values=color_class_PL,
                     labels=c("COV"="Covered","UP"="Up-bent PL","DOWN"= "Down-bent PL", "PL"="PL","DEG"="Degraded"))+
  scale_y_log10()

ggsave("../Figures/2_species/CA/Max_patch_size_by_PL_class.pdf",p2,width = 7,height = 5)



#Exponent PL by species

d_PL_lambda_sp=d_PL%>%
  group_by(., Species,a21,Stress)%>%
  summarise(.groups = "keep",Alpha_expo=mean(Alpha_expo))

p3=ggplot(NULL)+
  geom_line(data=d_PL_lambda_sp,aes(x=Stress,Alpha_expo,color=Species),lwd=.9)+
  geom_point(data=d_PL,aes(x=Stress,Alpha_expo,color=Species),size=.5,alpha=.3)+
  the_theme+labs(x="Stress (S)",  y = TeX(r'(PL exponent \ \  $\lambda)'),color="")+
  facet_wrap(.~a21,labeller = as_labeller(appender, default = label_parsed))+
  scale_color_manual(values=color_sp,
                     labels=c("+"="All vegetation","1"="Stress_tol","2"= "Competitive"))

ggsave("../Figures/2_species/CA/PL_exponent_by_species.pdf",p3,width = 7,height = 4)


#Exponent PL by PL class
d_PL_lambda_class=d_PL%>%
  group_by(., Species,a21,Stress,Class)%>%
  summarise(.groups = "keep",Alpha_expo=mean(Alpha_expo))

d_PL_lambda_class=transform(d_PL_lambda_class,
                            a21 = factor(a21, levels=c(1,2.5,4), labels=c("frac(alpha[0,1],alpha[1,0]) : 1", "frac(alpha[0,1],alpha[1,0]) : 2.5",
                                                                          "frac(alpha[0,1],alpha[1,0]) : 4")),
                            Species=factor(Species,levels=c("+","1","2"),labels=c("psi[0] + psi[1]","psi[1]","psi[0]")))

d_PL2=transform(d_PL,
                a21 = factor(a21, levels=c(1,2.5,4), labels=c("frac(alpha[0,1],alpha[1,0]) : 1", "frac(alpha[0,1],alpha[1,0]) : 2.5",
                                                              "frac(alpha[0,1],alpha[1,0]) : 4")),
                Species=factor(Species,levels=c("+","1","2"),labels=c("psi[0] + psi[1]","psi[1]","psi[0]")))
p4=ggplot(NULL)+
  geom_line(data=d_PL_lambda_class,aes(x=Stress,Alpha_expo,color=Class,group=interaction(Class,Species,a21)),lwd=.9)+
  geom_point(data=d_PL2%>%mutate(., a21=as.character(a21)),aes(x=Stress,Alpha_expo,color=Class),alpha=.3,size=.5)+
  the_theme+labs(x="Stress (S)", y = TeX(r'(PL exponent \ \  $\lambda)'),color="")+
  facet_grid(Species~a21,labeller=label_parsed)+
  scale_color_manual(values=color_class_PL,
                     labels=c("COV"="Covered","UP"="Up-bent PL","DOWN"= "Down-bent PL", "PL"="PL","DEG"="Degraded"))

ggsave("../Figures/2_species/CA/PL_exponent_by_PL_class.pdf",p4,width = 7,height = 5)

  
#
### d) Analysing species dynamics ----

mapping(100,100)
S_seq = seq(0, 0.8, length.out=30)
n_save = 25
relativ_comp = c("1.0", "2.5", "4.0")
d_vege=tibble()
#doing the loop on simulations
for (scena in c("inter=intra","low_inter")[-2]){
  
  for (stress in S_seq){

    stress=round(stress, 2)
    
    if (stress == 0){
      stress="0.0"
    }
    if (stress == 1){
      stress="1.0"
    }
    
    for (rela_c in relativ_comp){
      
      for (n in 1:n_save){
        
        
        #get landscape
        landscape=as.matrix(read.table(paste0("../Table/2_species/Patch_size/Big_sim/Landscape_l_S_",stress,
                                              "_a21_", rela_c,"_nsave_", n , ".csv"),sep=","))
        
        d_vege=rbind(d_vege,tibble(Rho_1=sum(landscape==1)/length(landscape),
                               Rho_2=sum(landscape==2)/length(landscape),Rho_p=sum(landscape %in% c(1,2))/length(landscape),N_rep=n,
                               a21=rela_c,Stress=as.numeric(stress)))
      } #end loop replicates
    } #end loop a12
  } # end loop stress
} #end loop scenarios competition


d_vege_merged=d_vege%>%
  group_by(.,a21,Stress)%>%
  summarise(.,Rho_1=mean(Rho_1),Rho_2=mean(Rho_2),Rho_p=mean(Rho_p),.groups = "keep" )

color_sp=c("Rho_p"="gray70","Rho_1"=as.character(color_rho[4]),"Rho_2"= as.character(color_rho[2]))
appender <- function(string) {
  TeX(paste("$\\frac{\\alpha_{0,1}}{\\alpha_{1,0}} = $", string))}



p=ggplot(NULL)+
  geom_path(data=d_vege_merged%>% melt(., id.vars=c("a21","Stress")),
            aes(x=Stress,value,color=variable,group=variable),lwd=1)+
  geom_path(data=d_vege%>%melt(., id.vars=c("a21","Stress","N_rep")),
            aes(x=Stress,value,color=variable,group=interaction(variable,N_rep)),alpha=.17)+
  the_theme+labs(x="Stress (S)", y="Patch density",color="")+
  facet_wrap(.~a21,labeller = as_labeller(appender, default = label_parsed))+
  scale_color_manual(values=color_sp,
                     labels=c("Rho_p"="All vegetation","Rho_1"="Stress_tol","Rho_2"= "Competitive"))

ggsave("../Figures/2_species/CA/Species_dynamics_CA_inter=intra.pdf",p,width = 7,height = 4)




#
# Step 4) Testing EWS on dynamics ----

rm(list = ls())
source("./2_Species_analysis_functions.R")
julia_setup()
de = diffeq_setup()


the_theme = theme_classic() + theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#CCE8D8"),
    strip.text.y = element_text(size = 10, angle = -90),
    strip.text.x = element_text(size = 8), axis.text = element_text(size = 11), axis.title = element_text(size = 13),
    legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)

plot_dynamics = function(d) {
    color_rho = c("fertile" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")

    if ("time" %in% colnames(d) | "Time" %in% colnames(d)) {
        return(ggplot(d %>% melt(., id.vars = colnames(d)[ncol(d)]) %>%
            mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive", "rho_0" = "fertile", "rho_d" = "desert"))) +
            geom_line(aes(x = time, y = value, color = variable), lwd = 1) +
            the_theme +
            scale_color_manual(values = color_rho) +
            labs(x = "Time", y = "Densities", color = ""))
    } else {
        d$time = 1:nrow(d)
        return(ggplot(d %>% melt(., id.vars = colnames(d)[ncol(d)]) %>%
            mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive", "rho_0" = "fertile", "rho_d" = "desert"))) +
            geom_line(aes(x = time, y = value, color = variable), lwd = 1) +
            the_theme +
            scale_color_manual(values = color_rho) +
            labs(x = "Time", y = "Densities", color = ""))
    }
}




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

Noise_MF = julia_eval("

function Noise_MF_2_species(du,u,p,t)

  # sigma = 0.01

  rho_1,rho_2,rho_d,rho_0=u

  du[1] = 0.001*u[1]
  du[2] = 0.001*u[2]
  du[3] = 0.001*u[3]
  du[4] = 0.001*u[4]


end")


param = c(r = 0.05, d = 0.1, f = 0.9, beta = 0.8, m = 0.1, e = .1, emax = 1.2, cintra = .1, cinter1 = .1, cinter2 = .1, S = 0)
state = c(rho_1 = 0.4, rho_2 = 0.4, rho_d = 0.1, rho_0 = 0.1)
tspan = c(0, 6000)
t = seq(0, 6000, by = 1)

julia_library("DifferentialEquations")
julia_assign("state", state)
julia_assign("tspan", tspan)
julia_assign("p", param)

C_seq = seq(0, 1, length.out = 30)
d2 = tibble()

# pdf("../Figures/2_species/EWS_dynamics.pdf",width = 6,height = 3)
for (S in C_seq) {
    d3 = tibble()
    for (rep in 1:20) {
        param[11] = S
        julia_assign("p", param)

        prob = julia_eval("SDEProblem(MF_two_species,Noise_MF_2_species, state, tspan, p)")
        prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
        sol = de$solve(prob, de$Tsit5(), saveat = t)
        d = as.data.frame(t(sapply(sol$u, identity)))
        colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
        plot_dynamics(d)
        d$rho_plus = d$rho_1 + d$rho_2

        d = d[((nrow(d) - 100):nrow(d)), ]

        # print(plot_dynamics(d))

        # getting the variance of species abundances
        d3 = rbind(d3, tibble(
            S = S, acf_1 = acf(d$rho_1, lag = 1, pl = F)$acf[2],
            acf_2 = acf(d$rho_2, lag = 1, pl = F)$acf[2],
            acf_tot = acf(d$rho_plus, lag = 1, pl = F)$acf[2],
            var_1 = mean(d$rho_1) / var(d$rho_1),
            var_2 = mean(d$rho_2) / var(d$rho_2),
            var_plus = mean(d$rho_plus) / var(d$rho_plus)
        ))
    }
    d2 = rbind(d2, colMeans(d3))

    # Just to check if there is a shift
    #
    #   prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
    #   sol = de$solve(prob,de$Tsit5(),saveat=t)
    #   d=as.data.frame(t(sapply(sol$u,identity)))
    #   d$rho_plus=d[,1]+d[,2] # sum species
    #   colnames(d)=c("rho_1","rho_2","rho_d","rho_0","rho_plus")
    #   d2=rbind(d2,d[nrow(d),]%>%add_column(S=S))
}
dev.off()

ggplot(d2 %>% add_column(., S = C_seq) %>% melt(., id.vars = "S")) +
    geom_line(aes(x = S, y = value, color = variable)) +
    theme_classic()




