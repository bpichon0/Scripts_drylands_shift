x = c("tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", 
      "JuliaCall", "diffeqr", "simecol", "tseries","phaseR","ggpattern",
      "ggquiver", "scales","boot","spatialwarnings","hillR","RColorBrewer","ggnewscale")
lapply(x, require, character.only = TRUE)


the_theme = theme_classic() + theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#CCE8D8"),
    strip.text.y = element_text(size = 10, angle = -90),
    strip.text.x = element_text(size = 8), axis.text = element_text(size = 11), axis.title = element_text(size = 13),
    legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)


color_rho = c("coexistence" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")

color_rho = c("coexistence" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")
color_Nsp=colorRampPalette(c("#077D10",as.character(color_rho[2]),as.character(color_rho[4]),"#2A39EF"))

#creating folders
dir.create("../",showWarnings = F)
dir.create("../Table",showWarnings = F)
dir.create("../Table/2_species",showWarnings = F)
dir.create("../Figures",showWarnings = F)
dir.create("../Figures/2_species",showWarnings = F)
dir.create("../Figures/2_species/CA",showWarnings = F)
dir.create("../Figures/2_species/PA",showWarnings = F)
dir.create("../Figures/2_species/MF",showWarnings = F)
dir.create("../Table/2_species/CA/PL_summary",showWarnings = F)
dir.create("../Table/2_species/CA/EWS_spatial",showWarnings = F)



# A) Mean field analysis ----

Get_MF_parameters = function() {
    return(c(
        r = 0.01, d = .025, f = 0.9, beta = 1, m = 0.15, e = .1,  cintra = .3,
        alpha_0 = .3, S = 0,h=1
    ))
}


Get_MF_initial_state = function(type="equal") {
  
  if (is.character(type)){
    if (type=="equal"){
      state = c(rho_1 = 0.4, rho_2 = 0.4, rho_m = 0.1, rho_0 = .1)
      names(state) = c("rho_1", "rho_2","rho_m","rho_0")
    } else if (type=="random_D"){
      state1=runif(1)*.8
      state=c(state1,.8-state1,.1,.1)
    } else if (type=="random_R"){
      state1=runif(1)*.01
      state=c(state1,.01-state1,.49,.5)
      
    } 
    }else{
    state=c(rho_1 = type[1], rho_2 = type[2], rho_m = type[3], rho_0 = 1-sum(type))
    
  }
  return(state)
}



# normal function
MF_two_species_julia = julia_eval("

function MF_two_species(du,u,p,t)

  r,d,f,beta,m,e,cintra,alpha_0,S,h=p
  rho_1,rho_2,rho_d,rho_0=u
  
  du[1] = rho_0 * ( rho_1 *  (beta * ( 1 - S * (1 - e)) - (cintra*rho_1 + alpha_0 * (1+h*exp(-1)) * rho_2))) - rho_1 * m
  du[2] = rho_0 * ( rho_2 *  (beta *( 1 - S) - (cintra*rho_2 + alpha_0*rho_1 ))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1)) - rho_0 * ( rho_1 *  (beta * ( 1 - S * (1 - e)) - (cintra*rho_1 + alpha_0 * (1+h*exp(-1)) * rho_2 ))) -
      rho_0 * ( rho_2 *  (beta *( 1 - S) - (cintra*rho_2 + alpha_0*rho_1 ))) + rho_2 * m + rho_1 * m

end")


# function of press perturbation
MF_two_species_julia_press = julia_eval("

function MF_two_species_press(du,u,p,t)

  r,d,f,beta1,beta2,m,e,cintra,alpha_0,S,h=p
  rho_1,rho_2,rho_d,rho_0=u
  
  du[1] = rho_0 * ( rho_1 *  (beta1 * ( 1 - S * (1 - e)) - (cintra*rho_1 + alpha_0 * (1+h*exp(-1)) * rho_2))) - rho_1 * m
  du[2] = rho_0 * ( rho_2 *  (beta2 *( 1 - S) - (cintra*rho_2 + alpha_0*rho_1 ))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(rho_1))
  du[4] = -d*rho_0 + rho_d*(r+f*(rho_1)) - rho_0 * ( rho_1 *  (beta1 * ( 1 - S * (1 - e)) - (cintra*rho_1 + alpha_0 * (1+h*exp(-1)) * rho_2 ))) -
      rho_0 * ( rho_2 *  (beta2 *( 1 - S) - (cintra*rho_2 + alpha_0*rho_1 ))) + rho_2 * m + rho_1 * m

end")


#MF varying trait
MF_two_species_julia_varying_trait = julia_eval("

function MF_two_species_varying_trait(du,u,p,t)

  r,d,f,beta,m,e,cintra,alpha_0,S,h,psi_1,psi_2=p
  rho_1,rho_2,rho_d,rho_0=u
  
  du[1] = rho_0 * ( rho_1 *  (beta * ( 1 - S * (1 - psi_1 * e)) - (cintra*rho_1 + alpha_0 * (1 + psi_1 * h * exp( - abs(psi_1 - psi_2))) * rho_2))) - rho_1 * m
  du[2] = rho_0 * ( rho_2 *  (beta * ( 1 - S * (1 - psi_2 * e)) - (cintra*rho_2 + alpha_0 * (1 + psi_2 * h * exp( - abs(psi_1 - psi_2))) * rho_1))) - rho_2 * m
  du[3] = d*rho_0 - rho_d*(r+f*(psi_1*rho_1+psi_2*rho_2))
  du[4] = -d*rho_0 + rho_d*(r+f*(psi_1*rho_1+psi_2*rho_2)) - rho_0 * ( rho_1 *  (beta * ( 1 - S * (1 - psi_1 * e)) - (cintra*rho_1 + alpha_0 * (1 + psi_1 * h * exp( - abs(psi_1 - psi_2))) * rho_2 ))) -
      rho_0 * ( rho_2 *  (beta *( 1 - S * (1 - psi_2 * e)) - (cintra*rho_2 + alpha_0 * (1 + psi_2 * h * exp( - abs(psi_1 - psi_2))) * rho_1 ))) + rho_2 * m + rho_1 * m

end")


# Hysteresis functions


Run_dynamics_hysteresis = function(plot = F, N_seq_c = 100, N_seq_S = 100, write_table = T) {
    color_rho = c("coexistence" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")

    tspan = c(0, 5000)
    t = seq(0, 5000, by = 1)
    julia_library("DifferentialEquations")
    julia_assign("tspan", tspan)

    param=Get_MF_parameters()    
    N_seq_s = N_seq_S
    d2 = tibble()
    c_seq = seq(0,0.3, length.out = N_seq_c)
    param["cintra"]=.3


    for (comp in c_seq) { # for each competition coefficient
        param["alpha_0"] = comp

        for (branch in c("Degradation", "Restoration")) { # do the bifurcation plot

            if (branch == "Degradation") { # forward
                state = c(rho_1 = 0.4, rho_2 = 0.4, rho_d = 0.1, rho_0 = 0.1)
                julia_assign("state", state)
                S_seq = seq(0, 1, length.out = N_seq_s)
            } else { # i.e. backward
                state = c(rho_1 = 0.005, rho_2 = 0.005, rho_d = 0.49, rho_0 = 0.5)
                julia_assign("state", state)
                S_seq = rev(seq(0, 1, length.out = N_seq_s))
            }

            for (s in S_seq) { # each stress value

                param["S"] = s
                julia_assign("p", param)
                prob = julia_eval("ODEProblem(MF_two_species, state, tspan, p)")
                sol = de$solve(prob, de$Tsit5(), saveat = t)
                d = as.data.frame(t(sapply(sol$u, identity)))
                colnames(d) = c("rho_1", "rho_2", "rho_d", "rho_0")
                d2 = rbind(d2, d[nrow(d), ] %>% add_column(S = s, alpha_0 = comp, orientation = branch))
            }
        }
    }

    d2[d2 < 10^-3] = 0
    d2$rho_plus = d2$rho_1 + d2$rho_2

    if (plot) {

        # global vegetation
        p = ggplot(d2) +
            geom_line(aes(x = S, y = rho_plus, linetype = orientation), lwd = .6, alpha = .8) +
            the_theme +
            facet_wrap(. ~ alpha_0, labeller = label_bquote(cols=alpha[e]==.(alpha_0))) +
            scale_linetype_manual(values = c(1, 2)) +
            labs(x = "Stress (S)", y = TeX("$\\rho_{+}$"), color = "",linetype="")

        # by species
        p2 = ggplot(d2 %>% melt(., measure.vars = c("rho_1", "rho_2")) %>%
            mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive"))) +
            geom_line(aes(x = S, y = value, color = variable, linetype = orientation), lwd = .6) +
            the_theme +
          facet_wrap(. ~ alpha_0, labeller = label_bquote(cols=alpha[e]==.(alpha_0))) +
          scale_linetype_manual(values = c(1, 2)) +
            scale_color_manual(values = color_rho[c(2, 4)]) +
            labs(x = "Stress (S)", y = TeX("$\\rho_{+_1}, \\rho_{+_2}$"), color = "",linetype="")

        p_tot = ggarrange(p, p2, nrow = 2)
        ggsave("../Figures/2_species/MF/Example_hysteresis.pdf", width = 7, height = 6)
    }

    if (write_table) {
        write.table(d2, paste0("../Table/2_species/Hysteresis_size_length_C_", N_seq_c, "_length_S_", N_seq_S, ".csv"), sep = ";")
    }
    return(d2)
}



Compute_hysteresis = function(d, n_species = 2) {
    nb_S = 2 * length(unique(d$S)) # we want both the degradation & restoration orientation
    d_hysteresis = tibble()

    for (i in 1:(nrow(d) / (nb_S))) { # for each bifurcation plot

        interval = ((i - 1) * nb_S + 1):(i * nb_S)

        for (sp in 1:n_species) {
            d2 = d[interval, c(paste0("rho_", sp), "S", "orientation")]
            # hysteresis size
            biomass_D = filter(d2, orientation == "Degradation")
            biomass_R = filter(d2, orientation == "Restoration")

            if (sp == 1) { # we delete the left part of the bifurcation diagram before species become dominant
                max_S = biomass_D$S[which(biomass_D[, 1] == max(biomass_D[, 1]))]
                biomass_D = filter(biomass_D, S > max_S)
                biomass_R = filter(biomass_R, S > max_S)
            }

            hysteresis = sum(biomass_D[, 1] - rev(biomass_R[, 1]))
            if (max(abs(diff(biomass_D[, 1]))) > .1) {
                # #now tipping points
                Tipping_D = biomass_D$S[which(abs(diff(biomass_D[, 1])) == max(abs(diff(biomass_D[, 1]))))]
                Tipping_R = biomass_R$S[which(abs(diff(biomass_R[, 1])) == max(abs(diff(biomass_R[, 1]))))]
            } else {
                Tipping_D = 0
                Tipping_R = 0
            }

            d_hysteresis = rbind(d_hysteresis, tibble(Hysteresis = hysteresis, Species = sp, Tipping_D = Tipping_D, Tipping_R = Tipping_R))
        }
    }
    return(d_hysteresis)
}



plot_dynamics = function(d) {
    color_rho2 = c("fertile" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")
    colnames(d) = c("rho_1", "rho_2", "rho_0", "rho_d")
    if ("time" %in% colnames(d) | "Time" %in% colnames(d)) {
        return(ggplot(d %>% melt(., id.vars = colnames(d)[ncol(d)]) %>%
            mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive", "rho_0" = "fertile", "rho_d" = "desert"))) +
            geom_line(aes(x = time, y = value, color = variable), lwd = 1) +
            theme_classic() +
            scale_color_manual(values = color_rho2) +
            labs(x = "Time", y = "Densities", color = "") +
            theme(legend.text = element_text(size = 11), legend.position = "bottom"))
    } else {
        d$time = 1:nrow(d)
        return(ggplot(d %>% melt(., id.vars = colnames(d)[ncol(d)]) %>%
            mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive", "rho_0" = "fertile", "rho_d" = "desert"))) +
            geom_line(aes(x = time, y = value, color = variable), lwd = 1) +
            the_theme +
            scale_color_manual(values = color_rho2) +
            labs(x = "Time", y = "Densities", color = "") +
            theme(legend.text = element_text(size = 11), legend.position = "bottom"))
    }
}


# Recruitment rate

Get_recruitment_rates = function(eq, param) {
    recruit_rate = c(
        r1 =  (1 - eq$rho_1 - eq$rho_2 - eq$rho_d) * param["beta"] * eq$rho_1 * ( (1 - param["S"] * (1 - param["e"])) - (param["cintra"] * eq$rho_1 + param["cinter1"] * eq$rho_2)) - eq$rho_1 * param["m"],
        r2 =  (1 - eq$rho_1 - eq$rho_2 - eq$rho_d) * param["beta"] * eq$rho_2 * ( (1 - param["S"]) - (param["cintra"] * eq$rho_2 + param["cinter2"] * eq$rho_1)) - eq$rho_2 * param["m"]
    )

    names(recruit_rate) = c("r1", "r2")
    return(recruit_rate)
}


# B) Pair approximation analysis ----

Get_PA_parameters = function() {
  param=c(
    r = 0.01, d = 0.025, f = .9, beta = 1, m = 0.15, e = .1, cintra=.3,alpha_0=.3, S = 0, delta = .1, z = 4,h=1
  )       
  
  return(param)
}

Get_PA_initial_state = function(ini=c(.4,.4,.1)) {
  
    if (length(ini==4)){ini=ini[-4]}
  
    state_pair = c(
        rho_12 = ini[1] * ini[2], rho_1m = ini[1] * ini[3], rho_2m = ini[2] * ini[3],
        rho_11 = ini[1] * ini[1], rho_22 = ini[2] * ini[2], rho_mm = ini[3] * ini[3]
    )
    ini = c(ini, state_pair)
    names(ini) = c("rho_1", "rho_2", "rho_m", "rho_12", "rho_1m", "rho_2m", "rho_11", "rho_22", "rho_mm")
    return(ini)
}

get_mean_densities = function(d) {
    length_transient = 3000
    d = d[(length_transient + 1):nrow(d), ]
    return(colMeans(d))
}


## a) local facilitation ----

PA_two_species_local_C_local_F = julia_eval("

function PA_two_species_local_C_local_F(du,u,p,t)


r,d,f,beta,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0=(1- rho_1-rho_2-rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m 



#rho_1
du[1] =  rho_0 *  (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
( (cintra * (rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) ) ) - rho_1 * m

#rho_2
du[2] =  rho_0 *  (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S)  - 
( (cintra * (rho_20/rho_0) + alpha_0 * (rho_10/rho_0))  )  ) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * (rho_1m / rho_m))

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * 
        (beta * (1 - S ) - ( (cintra * ((z - 1) / z)*(rho_20/rho_0) + (alpha_0/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) )   )  +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * 
        (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        ( (cintra * ((z - 1) / z) * (rho_10/rho_0) + ((alpha_0 * (1+h*exp(-1)))/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0)) )  ) -
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        (cintra * ((z - 1) / z) * (rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0))     ) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_1m / rho_m) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S  ) - 
        ( (cintra * ((z - 1) / z) * (rho_20/rho_0) + alpha_0 *  ((z - 1) / z) *(rho_10/rho_0)) )) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_1m / rho_m )))

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta *  (1 - S * (1-e )) - (cintra * ((z - 1) / z) * (rho_10/rho_0) + (cintra/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0))  ) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta *  
        (1 - S ) - ( (cintra * ((z - 1) / z)*(rho_20/rho_0) + (cintra/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) )    ) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_1m / rho_m))


end


")

PA_two_species_global_C_local_F = julia_eval("

function PA_two_species_global_C_local_F(du,u,p,t)

r,d,f,beta,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0=(1- rho_1-rho_2-rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m 


#rho_1
du[1] =  rho_0 *  (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
(cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) ) - rho_1 * m

#rho_2
du[2] =  rho_0 *  (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S)  - 
(cintra *rho_2 + alpha_0 * rho_1)) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * (rho_1m / rho_m))

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * 
        (beta * (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1))  +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * 
        (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) )-
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_1m / rho_m) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S  ) - 
        ((cintra *rho_2 + alpha_0 *rho_1) )) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_1m / rho_m )))

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta *  (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta *  
        (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_1m / rho_m))

end


")





PA_two_species_local_C_local_F_press = julia_eval("

function PA_two_species_local_C_local_F_press(du,u,p,t)

r,d,f,beta1,beta2,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0  = (1- rho_1-rho_2-rho_m)
rho_10 = (rho_1 - rho_11 - rho_12 - rho_1m)
rho_20 = (rho_2 - rho_22 - rho_12 - rho_2m)
rho_0m = (rho_m - rho_mm - rho_1m - rho_2m)

#rho_1
du[1] =  rho_0 * (delta * rho_1 + (1 - delta) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) )) - rho_1 * m

#rho_2
du[2] =  rho_0 * (delta * rho_2 + (1 - delta) * ((rho_20) / rho_0)) * (beta2 * (1 - S)  - ( (cintra *(rho_20/rho_0) + alpha_0 *(rho_10/rho_0)) )) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * (rho_1m / rho_m))

#rho_12
du[4] = (rho_10) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0)) * (beta2 * (1 - S ) - ( (cintra *((z - 1) / z)*(rho_20/rho_0) + (alpha_0/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) ))  +
        (rho_20) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 * ((z - 1) / z) * (rho_10/rho_0) + (alpha_0 * (1+h*exp(-1))/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0)) ))-
        2 * rho_12 * m

#rho_1m
du[5] = (rho_10) * d + rho_0m * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) )) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_1m / rho_m) ))

#rho_2m
du[6] = (rho_20) * d + rho_0m * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0) * (beta2 * (1 - S  ) - ( (cintra *(rho_20/rho_0) + alpha_0 *(rho_10/rho_0)) ))) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_1m / rho_m )))

#rho_11
du[7] = 2* (rho_10) * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0 )) * ( beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) * ((z - 1) / z) + (alpha_0/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) * (rho_20/rho_0)) )) -
        2 * rho_11 * m

#rho_22
du[8] = 2* (rho_20) * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0 )) * (beta2 * (1 - S ) - ( (cintra *(rho_20/rho_0) * ((z - 1) / z) + (cintra/z) + alpha_0 * ((z - 1) / z) * (rho_10/rho_0)) )) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * rho_0m * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_1m / rho_m))

end


")




PA_two_species_global_C_local_F_press = julia_eval("

function PA_two_species_global_C_local_F_press(du,u,p,t)

r,d,f,beta1,beta2,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0=(1- rho_1-rho_2-rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m 


#rho_1
du[1] =  rho_0 *  (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta1 * (1 - S * (1-e )) - 
(cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) ) - rho_1 * m

#rho_2
du[2] =  rho_0 *  (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta2 * (1 - S)  - 
(cintra *rho_2 + alpha_0 * rho_1)) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * (rho_1m / rho_m))

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * 
        (beta2 * (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1))  +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * 
        (rho_10 / rho_0)) * (beta1 * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) )-
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta1 * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_1m / rho_m) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta2 * (1 - S  ) - 
        ((cintra *rho_2 + alpha_0 *rho_1) )) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_1m / rho_m )))

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta1 *  (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta2 *  
        (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_1m / rho_m))

end


")


## b) global facilitation ----

PA_two_species_local_C_global_F = julia_eval("

function PA_two_species_local_C_global_F(du,u,p,t)


r,d,f,beta,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0=(1- rho_1-rho_2-rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m 



#rho_1
du[1] =  rho_0 *  (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
( (cintra * (rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) ) ) - rho_1 * m

#rho_2
du[2] =  rho_0 *  (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S)  - 
( (cintra * (rho_20/rho_0) + alpha_0 *(rho_10/rho_0))  )  ) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * rho_1)

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * 
        (beta * (1 - S ) - ( (cintra * ((z - 1) / z)*(rho_20/rho_0) + (alpha_0/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) )   )  +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * 
        (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        ( (cintra * ((z - 1) / z) * (rho_10/rho_0) + ((alpha_0 * (1+h*exp(-1)))/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0)) )  ) -
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        (cintra * ((z - 1) / z) * (rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0))     ) -
        rho_1m * m - rho_1m * (r + f*rho_1)

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S  ) - 
        ( (cintra * ((z - 1) / z) * (rho_20/rho_0) + alpha_0 *  ((z - 1) / z) *(rho_10/rho_0)) )) -
        rho_2m * m - rho_2m * (r + f*rho_1)

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta *  (1 - S * (1-e )) - (cintra * ((z - 1) / z) * (rho_10/rho_0) + (cintra/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0))  ) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta *  
        (1 - S ) - ( (cintra * ((z - 1) / z)*(rho_20/rho_0) + (cintra/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) )    ) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * rho_1)


end


")

PA_two_species_global_C_global_F = julia_eval("

function PA_two_species_global_C_global_F(du,u,p,t)

r,d,f,beta,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0=(1- rho_1-rho_2-rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m 


#rho_1
du[1] =  rho_0 *  (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) ) - rho_1 * m

#rho_2
du[2] =  rho_0 *  (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S)  - 
(cintra *rho_2 + alpha_0 * rho_1)) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * rho_1)

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * 
        (beta * (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1))  +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * 
        (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) )-
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        rho_1m * m - rho_1m * (r + f*rho_1)

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S  ) - 
        ((cintra *rho_2 + alpha_0 *rho_1) )) -
        rho_2m * m - rho_2m * (r + f*rho_1)

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta *  (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta *  
        (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * rho_1)

end


")




PA_two_species_local_C_global_F_press = julia_eval("

function PA_two_species_local_C_global_F_press(du,u,p,t)

r,d,f,beta1,beta2,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

#rho_1
du[1] =  rho_0 * (delta * rho_1 + (1 - delta) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) )) - rho_1 * m

#rho_2
du[2] =  rho_0 * (delta * rho_2 + (1 - delta) * ((rho_20) / rho_0)) * (beta2 * (1 - S)  - ( (cintra *(rho_20/rho_0) + alpha_0 *(rho_10/rho_0)) )) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * rho_1)

#rho_12
du[4] = (rho_10) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0)) * (beta2 * (1 - S ) - ( (cintra *((z - 1) / z)*(rho_20/rho_0) + (alpha_0/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) ))  +
        (rho_20) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 * ((z - 1) / z) * (rho_10/rho_0) + (alpha_0 * (1+h*exp(-1))/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0)) ))-
        2 * rho_12 * m

#rho_1m
du[5] = (rho_10) * d + rho_0m * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) )) -
        rho_1m * m - rho_1m * (r + f*rho_1)

#rho_2m
du[6] = (rho_20) * d + rho_0m * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0) * (beta2 * (1 - S  ) - ( (cintra *(rho_20/rho_0) + alpha_0 *(rho_10/rho_0)) ))) -
        rho_2m * m - rho_2m * (r + f*rho_1)

#rho_11
du[7] = 2* (rho_10) *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0 )) * (beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) * ((z - 1) / z) + (alpha_0/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) * (rho_20/rho_0)) )) -
        2 * rho_11 * m

#rho_22
du[8] = 2* (rho_20) * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0 )) * (beta2 * (1 - S ) - ( (cintra *(rho_20/rho_0) * ((z - 1) / z) + (cintra/z) + alpha_0 * ((z - 1) / z) * (rho_10/rho_0)) )) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * rho_0m * d - 2* rho_mm * (r + f * rho_1)

end


")



PA_two_species_global_C_global_F_press = julia_eval("

function PA_two_species_global_C_global_F_press(du,u,p,t)

r,d,f,beta1,beta2,m,e,cintra,alpha_0,S,delta,z,h=p
rho_1,rho_2,rho_m,rho_12,rho_1m,rho_2m,rho_11,rho_22,rho_mm=u

rho_0=(1- rho_1-rho_2-rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m 


#rho_1
du[1] =  rho_0 *  (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta1 * (1 - S * (1-e )) - 
(cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) ) - rho_1 * m

#rho_2
du[2] =  rho_0 *  (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta2 * (1 - S)  - 
(cintra *rho_2 + alpha_0 * rho_1)) - rho_2 * m

#rho_m
du[3] =  rho_0 * d - rho_m * (r + f * rho_1)

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * 
        (beta2 * (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1))  +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * 
        (rho_10 / rho_0)) * (beta1 * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2) )-
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta1 * (1 - S * (1-e )) - 
        (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        rho_1m * m - rho_1m * (r + f*rho_1)

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0) * (beta2 * (1 - S  ) - 
        ((cintra *rho_2 + alpha_0 *rho_1) ))) -
        rho_2m * m - rho_2m * (r + f*rho_1)

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta1 *  (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta2 *  
        (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * rho_1)

end


")




## c) varying trait ----

PA_two_species_varying_trait = julia_eval("

function PA_two_species_varying_trait(du,u,p,t)

r, d, f, beta, m, e, cintra, alpha_0, S, delta, z, h, psi_1,psi_2 = p
rho_1, rho_2, rho_m, rho_12, rho_1m, rho_2m, rho_11, rho_22, rho_mm = u

rho_0 = (1 - rho_1 - rho_2 - rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m


#rho_1
du[1] = rho_0 * (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1)) -
                                                                        (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) - rho_1 * m

#rho_2
du[2] = rho_0 * (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2)) -
                                                                        (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) - rho_2 * m

#rho_m
du[3] = rho_0 * d - rho_m * (r + f * ( psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m)))

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_2)) - (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) *
                                      (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1)) -
                                                               (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1)) -
                                                                                                            (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        rho_1m * m - rho_1m * (r + f * ((1 / z) * psi_1 + ((z - 1) / z) * (psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m)) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2)) -
                                                                                                           (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        rho_2m * m - rho_2m * (r + f * ((1 / z) * psi_2 + ((z - 1) / z) * (psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m))))

#rho_11
du[7] = 2 * rho_10 * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_1 )) - (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2 * rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta *
                                                                                                                     (1 - S * (1 - e * psi_2)) - (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2 * rho_mm * (r + f * ((z - 1) / z) * (psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m)))

end


")

PA_two_species_varying_trait_trade_off = julia_eval("

function PA_two_species_varying_trait_trade_off(du,u,p,t)

r, d, f, beta, m, e, cintra, alpha_0, S, delta, z, h, psi_1,psi_2,shape = p
rho_1, rho_2, rho_m, rho_12, rho_1m, rho_2m, rho_11, rho_22, rho_mm = u

rho_0 = (1 - rho_1 - rho_2 - rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m


#rho_1
du[1] = rho_0 * (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1^shape)) -
                                                                        (cintra * rho_1 + alpha_0 * (1 + psi_1^shape * h * exp(-abs(psi_1 - psi_2))) * rho_2)) - rho_1 * m

#rho_2
du[2] = rho_0 * (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2^shape)) -
                                                                        (cintra * rho_2 + alpha_0 * (1 + psi_2^shape * h * exp(-abs(psi_1 - psi_2))) * rho_1)) - rho_2 * m

#rho_m
du[3] = rho_0 * d - rho_m * (r + f * ( psi_1^shape * (rho_1m / rho_m) + psi_2^shape * (rho_2m / rho_m)))

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_2^shape)) - (cintra * rho_2 + alpha_0 * (1 + psi_2^shape * h * exp(-abs(psi_1 - psi_2))) * rho_1)) +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) *
                                      (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1^shape)) -
                                                               (cintra * rho_1 + alpha_0 * (1 + psi_1^shape * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1^shape)) -
                                                                                                            (cintra * rho_1 + alpha_0 * (1 + psi_1^shape * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        rho_1m * m - rho_1m * (r + f * ((1 / z) * psi_1^shape + ((z - 1) / z) * (psi_1^shape * (rho_1m / rho_m) + psi_2^shape * (rho_2m / rho_m)) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2^shape)) -
                                                                                                           (cintra * rho_2 + alpha_0 * (1 + psi_2^shape * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        rho_2m * m - rho_2m * (r + f * ((1 / z) * psi_2^shape + ((z - 1) / z) * (psi_1^shape * (rho_1m / rho_m) + psi_2^shape * (rho_2m / rho_m))))

#rho_11
du[7] = 2 * rho_10 * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_1^shape )) - (cintra * rho_1 + alpha_0 * (1 + psi_1^shape * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2 * rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta *
                                                                                                                     (1 - S * (1 - e * psi_2^shape)) - (cintra * rho_2 + alpha_0 * (1 + psi_2^shape * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2 * rho_mm * (r + f * ((z - 1) / z) * (psi_1^shape * (rho_1m / rho_m) + psi_2^shape * (rho_2m / rho_m)))

end


")


PA_two_species_varying_trait_press = julia_eval("

function PA_two_species_varying_trait_press(du,u,p,t)

r, d, f, beta, m, e, cintra, alpha_0, S, delta, z, h, psi_1,psi_2 = p
rho_1, rho_2, rho_m, rho_12, rho_1m, rho_2m, rho_11, rho_22, rho_mm = u

rho_0 = (1 - rho_1 - rho_2 - rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m


#rho_1
du[1] = rho_0 * (delta * rho_1 + (1 - delta) * (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1)) -
                                                                        (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) - rho_1 * m

#rho_2
du[2] = rho_0 * (delta * rho_2 + (1 - delta) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2)) -
                                                                        (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) - rho_2 * m

#rho_m
du[3] = rho_0 * d - rho_m * (r + f * ( psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m)))

#rho_12
du[4] = rho_10 * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_2)) - (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) +
        rho_20 * (delta * rho_1 + (1 - delta) * ((z - 1) / z) *
                                      (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1)) -
                                                               (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_12 * m

#rho_1m
du[5] = rho_10 * d + (rho_0m) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) * (beta * (1 - S * (1 - e * psi_1)) -
                                                                                                            (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        rho_1m * m - rho_1m * (r + f * ((1 / z) * psi_1 + ((z - 1) / z) * (psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m)) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2)) -
                                                                                                           (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        rho_2m * m - rho_2m * (r + f * ((1 / z) * psi_2 + ((z - 1) / z) * (psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m))))

#rho_11
du[7] = 2 * rho_10 * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_1 )) - (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2 * rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta *
                                                                                                                     (1 - S * (1 - e * psi_2)) - (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2 * rho_mm * (r + f * ((z - 1) / z) * (psi_1 * (rho_1m / rho_m) + psi_2 * (rho_2m / rho_m)))

end


")
