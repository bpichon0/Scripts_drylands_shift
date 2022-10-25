x = c("tidyverse", "ggpubr", "latex2exp", "deSolve", "reshape2", 
      "JuliaCall", "diffeqr", "simecol", "tseries","phaseR","ggpattern",
      "ggquiver", "scales","boot","spatialwarnings","hillR","RColorBrewer")
lapply(x, require, character.only = TRUE)


the_theme = theme_classic() + theme(
    legend.position = "bottom",
    strip.background = element_rect(fill = "#CCE8D8"),
    strip.text.y = element_text(size = 10, angle = -90),
    strip.text.x = element_text(size = 8), axis.text = element_text(size = 11), axis.title = element_text(size = 13),
    legend.text = element_text(size = 10), text = element_text(family = "NewCenturySchoolbook")
)


color_rho = c("coexistence" = "#D8CC7B", "competitive" = "#ACD87B", "desert" = "#696969", "stress_tol" = "#7BD8D3")


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
    } else{
      state1=runif(1)*.8
      state=c(state1,.8-state1,.1,.1)
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
            facet_wrap(. ~ alpha_0, labeller = label_both) +
            scale_linetype_manual(values = c(1, 2)) +
            labs(x = "Stress (S)", y = TeX("$\\rho_{+}$"), color = "",linetype="")

        # by species
        p2 = ggplot(d2 %>% melt(., measure.vars = c("rho_1", "rho_2")) %>%
            mutate(., variable = recode_factor(variable, "rho_1" = "stress_tol", "rho_2" = "competitive"))) +
            geom_line(aes(x = S, y = value, color = variable, linetype = orientation), lwd = .6) +
            the_theme +
            facet_wrap(. ~ alpha_0, labeller = label_both) +
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
du[3] =  rho_0 * d - rho_m * (r + f * (rho_10 / rho_0))

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
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_10 / rho_0) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S  ) - 
        ( (cintra * ((z - 1) / z) * (rho_20/rho_0) + alpha_0 *  ((z - 1) / z) *(rho_10/rho_0)) )) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_10 / rho_0 )))

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta *  (1 - S * (1-e )) - (cintra * ((z - 1) / z) * (rho_10/rho_0) + (cintra/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0))  ) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta *  
        (1 - S ) - ( (cintra * ((z - 1) / z)*(rho_20/rho_0) + (cintra/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) )    ) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_10 / rho_0))


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
du[3] =  rho_0 * d - rho_m * (r + f * (rho_10 / rho_0))

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
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_10 / rho_0) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S  ) - 
        ((cintra *rho_2 + alpha_0 *rho_1) )) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_10 / rho_0 )))

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta *  (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta *  
        (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_10 / rho_0))

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
du[3] =  rho_0 * d - rho_m * (r + f * ((rho_10) / rho_0))

#rho_12
du[4] = (rho_10) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0)) * (beta2 * (1 - S ) - ( (cintra *((z - 1) / z)*(rho_20/rho_0) + (alpha_0/z) + alpha_0 *((z - 1) / z)*(rho_10/rho_0)) ))  +
        (rho_20) * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 * ((z - 1) / z) * (rho_10/rho_0) + (alpha_0 * (1+h*exp(-1))/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) *(rho_20/rho_0)) ))-
        2 * rho_12 * m

#rho_1m
du[5] = (rho_10) * d + rho_0m * (delta * rho_1 + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0)) * (beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) + alpha_0 * (1+h*exp(-1)) *(rho_20/rho_0)) )) -
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * ((rho_10) / rho_0) ))

#rho_2m
du[6] = (rho_20) * d + rho_0m * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0) * (beta2 * (1 - S  ) - ( (cintra *(rho_20/rho_0) + alpha_0 *(rho_10/rho_0)) ))) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * ((rho_10) / rho_0 )))

#rho_11
du[7] = 2* (rho_10) * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_10) / rho_0 )) * ( beta1 * (1 - S * e ) - ( (alpha_0 *(rho_10/rho_0) * ((z - 1) / z) + (alpha_0/z) + alpha_0 * (1+h*exp(-1)) * ((z - 1) / z) * (rho_20/rho_0)) )) -
        2 * rho_11 * m

#rho_22
du[8] = 2* (rho_20) * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * ((rho_20) / rho_0 )) * (beta2 * (1 - S ) - ( (cintra *(rho_20/rho_0) * ((z - 1) / z) + (cintra/z) + alpha_0 * ((z - 1) / z) * (rho_10/rho_0)) )) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * rho_0m * d - 2* rho_mm * (r + f * ((z - 1) / z) * ((rho_10) / rho_0))

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
du[3] =  rho_0 * d - rho_m * (r + f * (rho_10 / rho_0))

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
        rho_1m * m - rho_1m * (r + f*( (1 / z) + ((z - 1) / z) * (rho_10 / rho_0) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta2 * (1 - S  ) - 
        ((cintra *rho_2 + alpha_0 *rho_1) )) -
        rho_2m * m - rho_2m * (r + f*( ((z - 1) / z) * (rho_10 / rho_0 )))

#rho_11
du[7] = 2* rho_10 *  (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0 )) * 
        (beta1 *  (1 - S * (1-e )) - (cintra *rho_1 + alpha_0 * (1+h*exp(-1))*rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2* rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0 )) * (beta2 *  
        (1 - S ) - (cintra *rho_2 + alpha_0 *rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2* rho_mm * (r + f * ((z - 1) / z) * (rho_10 / rho_0))

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
du[3] = rho_0 * d - rho_m * (r + f * ( psi_1 * (rho_10 / rho_0) + psi_2 * (rho_20 / rho_0)))

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
        rho_1m * m - rho_1m * (r + f * ((1 / z) * psi_1 + ((z - 1) / z) * (psi_1 * (rho_10 / rho_0) + psi_2 * (rho_20 / rho_0)) ))

#rho_2m
du[6] = rho_20 * d + (rho_0m) * (delta * rho_2 + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta * (1 - S * (1 - e * psi_2)) -
                                                                                                           (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        rho_2m * m - rho_2m * (r + f * ((1 / z) * psi_2 + ((z - 1) / z) * (psi_1 * (rho_10 / rho_0) + psi_2 * (rho_20 / rho_0))))

#rho_11
du[7] = 2 * rho_10 * (delta * rho_1 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_10 / rho_0)) *
        (beta * (1 - S * (1 - e * psi_1 )) - (cintra * rho_1 + alpha_0 * (1 + psi_1 * h * exp(-abs(psi_1 - psi_2))) * rho_2)) -
        2 * rho_11 * m

#rho_22
du[8] = 2 * rho_20 * (delta * rho_2 + ((1 - delta) / z) + (1 - delta) * ((z - 1) / z) * (rho_20 / rho_0)) * (beta *
                                                                                                                     (1 - S * (1 - e * psi_2)) - (cintra * rho_2 + alpha_0 * (1 + psi_2 * h * exp(-abs(psi_1 - psi_2))) * rho_1)) -
        2 * rho_22 * m

#rho_mm
du[9] = 2 * (rho_0m) * d - 2 * rho_mm * (r + f * ((z - 1) / z) * (psi_1 * (rho_10 / rho_0) + psi_2 * (rho_20 / rho_0)))

end


")

# C) CA analysis ----

Get_CA_parameters = function() {
    return(c(
        r = 0.05, d = 0.1, f = .9, beta = 0.8, m = 0.1, e = .05,  alpha11 = .1,
        alpha12 = .1, alpha21 = .1, alpha22 = .1, S = 0, delta = .1, z = 4, leap=0.5
    ))
}


CA_2_species = function(landscape, param) {

    # Variables : 1 = stress_tol, 2 = competitive, 0 = fertile, -1 = degraded
    rho_1 = length(which((landscape == 1))) / length(landscape) # fraction stress_tol
    rho_2 = length(which((landscape == 2))) / length(landscape) # fraction competitive
    rho_0 = length(which((landscape == 0))) / length(landscape) # fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_0 # fraction degraded

    # Neighbors :

    # using simcol package
    neigh_1 = simecol::neighbors(x = landscape, state = 1, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)
    neigh_2 = simecol::neighbors(x = landscape, state = 2, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)
    neigh_f = simecol::neighbors(x = landscape, state = 0, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)
    neigh_d = simecol::neighbors(x = landscape, state = -1, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)


    r = param[1]
    d = param[2]
    f = param[3]
    beta = param[4]
    m = param[5]
    e = param[6]
    emax = param[7]
    cg = param[8]
    alpha11 = param[9]
    alpha12 = param[10]
    alpha21 = param[11]
    alpha22 = param[12]
    S = param[13]
    delta = param[14]
    z = param[15]

    colonization_1 = beta * (delta * rho_1 + (1 - delta) * neigh_1 / z) * ((1 - S * (1 - e)) - (cg*(rho_1+rho_2)+(alpha11*(neigh_1/z) + alpha21*(neigh_2/z))))
    colonization_2 = beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * ((1 - S)           - (cg*(rho_1+rho_2)+(alpha22*(neigh_2/z) + alpha12*(neigh_1/z))))


    # calculate regeneration, degradation & mortality rate
    death = m
    regeneration = (r + f * (neigh_1) / z)
    degradation = d

    # Apply rules
    rnum = runif(length(landscape)) # one random number between 0 and 1 for each cell
    landscape_update = landscape

    ## New stress_tol
    landscape_update[which(landscape == 0 & rnum <= colonization_1)] = 1

    ## New competitive
    landscape_update[which(landscape == 0 & rnum > colonization_1 & rnum <= colonization_2 + colonization_1)] = 2

    ## New fertile
    landscape_update[which(landscape == 1 & rnum <= death)] = 0
    landscape_update[which(landscape == 2 & rnum <= death)] = 0
    landscape_update[which(landscape == -1 & rnum <= regeneration)] = 0

    ## New degraded
    landscape_update[which(landscape == 0 & rnum > colonization_2 + colonization_1 & rnum <= colonization_2 + colonization_1 + degradation)] = -1

    rho_1 = length(which((landscape == 1))) / length(landscape) # fraction stress_tol
    rho_2 = length(which((landscape == 2))) / length(landscape) # fraction competitive
    rho_0 = length(which((landscape == 0))) / length(landscape) # fraction fertile
    rho_d = 1 - rho_1 - rho_2 - rho_0 # fraction degraded

    # print(max(colonization_1))
    # print(max(colonization_2))
    # print(max(regeneration))


    if (any(c(
        colonization_1 + degradation,
        colonization_2 + degradation,
        colonization_2 + colonization_1 + degradation,
        death, regeneration
    ) > 1)) {
        warning("a set probability is exceeding 1 in run! decrease delta!!!")
    }
    if (any(c(
        colonization_2, colonization_1,
        degradation, regeneration,
        death
    ) < 0)) {
        warning("a set probability falls below 0 in run balance parameters!!!")
    }


    return(list(State = c(rho_1, rho_2, rho_0, rho_d), Landscape = landscape_update))
}



Run_CA_2_species = function(time = seq(1, 1000, 1), param, landscape) {
    d = tibble(
        rho_1 = sum(landscape == 1) / length(landscape), rho_2 = sum(landscape == 2) / length(landscape),
        rho_0 = sum(landscape == 0) / length(landscape), rho_d = sum(landscape == -1) / length(landscape), time = 1
    )

    for (k in 2:length(time)) {
        param["dt"] = time[k] - time[k - 1]
        output_CA = CA_2_species(landscape, param = param)
        landscape = output_CA$Landscape
        d = rbind(d, tibble(
            rho_1 = output_CA$State[1], rho_2 = output_CA$State[2],
            rho_0 = output_CA$State[3], rho_d = output_CA$State[4],
            time = k
        ))
    }

    return(list(state = d, landscape = landscape))
}


Get_initial_lattice = function(frac = c(.4, .4, .1, .1), size = 25) {
    return(matrix(sample(c(2, 1, 0, -1), replace = T, size = size * size, prob = frac), ncol = size, nrow = size))
}


Plot_landscape = function(landscape) {
    color_CA = c("1" = "#7BD8D3", "2" = "#ACD87B", "0" = "#D8CC7B", "-1" = "#696969")
    colnames(landscape)=rownames(landscape)=1:nrow(landscape)
    ggplot(melt(as.matrix(landscape))) +
        geom_tile(aes(x = Var1, y = Var2, fill = as.character(value))) +
        theme_transparent() +
        scale_fill_manual(values = color_CA, labels = c("Stress tol", "Competitive", "Fertile", "Desert")) +
        theme(panel.border = element_blank()) +
        theme(legend.position = "bottom") +
        labs(fill = "")
}




Gillespie_tau_leaping_R = function(landscape, time, param) {
    d2 = tibble()
    # param
    r = param[1]
    d = param[2]
    f = param[3]
    beta = param[4]
    m = param[5]
    e = param[6]
    emax = param[7]
    cg = param[8]
    alpha11 = param[9]
    alpha12 = param[10]
    alpha21 = param[11]
    alpha22 = param[12]
    S = param[13]
    delta = param[14]
    z = param[15]
    leap = param[16]

    Rate_landscape = array(0, c(nrow(landscape), ncol(landscape), 6))

    rules_change = matrix(c(0, 1, 0, 2, 0, -1, 1, 0, 2, 0, -1, 0), ncol = 2, nrow = 6, byrow = T)


    for (dt in 1:length(time)) { # for each time step

        # calculate all rates in the lattice

        # number of neighbors
        neigh_1 = simecol::neighbors(x = landscape, state = 1, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)
        neigh_2 = simecol::neighbors(x = landscape, state = 2, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)
        neigh_f = simecol::neighbors(x = landscape, state = 0, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)
        neigh_d = simecol::neighbors(x = landscape, state = -1, wdist = matrix(c(0, 1, 0, 1, 0, 1, 0, 1, 0), nrow = 3), bounds = 1)

        # global densities
        rho_1 = length(which((landscape == 1))) / length(landscape) # fraction stress_tol
        rho_2 = length(which((landscape == 2))) / length(landscape) # fraction competitive
        rho_0 = length(which((landscape == 0))) / length(landscape) # fraction fertile
        rho_d = 1 - rho_1 - rho_2 - rho_0 # fraction degraded

        # calculate the rates
        Rate_landscape[, , 1] = beta * (delta * rho_1 + (1 - delta) * neigh_1 / z) * ((1 - S * (1 - e)) - (cg*(rho_1+rho_2)+(alpha11*(neigh_1/z) + alpha21*(neigh_2/z))))
        Rate_landscape[, , 2] = beta * (delta * rho_2 + (1 - delta) * neigh_2 / z) * ((1 - S ) - (cg*(rho_1+rho_2)+(alpha22*(neigh_2/z) + alpha12*(neigh_1/z))))
        Rate_landscape[, , 3] = d
        Rate_landscape[, , 4] = m
        Rate_landscape[, , 5] = m
        Rate_landscape[, , 6] = (r + f * neigh_1 / z)
        
        Rate_landscape[,,1][Rate_landscape[,,1]<0]=0 #to avoid problems at high stress
        Rate_landscape[,,2][Rate_landscape[,,2]<0]=0 #to avoid problems at high stress
        
        
        # calculate their propensity (i.e. sum of type of events)
        Propensity = c(
            sum(Rate_landscape[, , 1][which(landscape == 0)]), # colonization stress tol
            sum(Rate_landscape[, , 2][which(landscape == 0)]), # colonization competitive
            sum(Rate_landscape[, , 3][which(landscape == 0)]), # degradation
            sum(Rate_landscape[, , 4][which(landscape == 1)]), # death stress_tol
            sum(Rate_landscape[, , 5][which(landscape == 2)]), # death competitive
            sum(Rate_landscape[, , 6][which(landscape == -1)]) # regeneration
        )
        
        nb_events = sapply(1:length(Propensity), function(x) {
            rpois(1, lambda = Propensity[x] * leap)
        }) # number of events per event type

        for (ev in 1:length(nb_events)) { # for each type of events

            patches = which(landscape == rules_change[ev, 1])

            if (nb_events[ev] != 0 & length(patches) >= nb_events[ev]) {
              
              proba_sample = Rate_landscape[, , ev][patches]
                
              
              if (length(patches)==1){
                landscape[patches] = rules_change[ev, 2]
              } else{
                landscape[sample(patches,prob = proba_sample,nb_events[ev], replace = T)] = rules_change[ev, 2]
              }
            }
        }


        d2 = rbind(d2, tibble(
            rho_1 = length(which((landscape == 1))) / length(landscape), # fraction stress_tol
            rho_2 = length(which((landscape == 2))) / length(landscape), # fraction competitive
            rho_0 = length(which((landscape == 0))) / length(landscape), # fraction fertile
            rho_d = 1 - rho_1 - rho_2 - rho_0 # fraction degraded
        ))
    }

    return(list(state = d2, landscape = landscape))
}


# D) Spatial analysis ----

#From Schneider et al., 2016 TE
Get_patches = function(landscape, state) {
  
  if (state=="+"){
    landscape[landscape %in% c(1,2)] ="+"
  }
  
  pattern = landscape
  pattern = pattern %in% state #keeping the state of interest
  map = rep(NA, times = length(landscape))
  old = rep(99, times = length(landscape)) #to compare
  
  while(!identical(old[pattern], map[pattern])) { 
    
    old = map
    count = as.integer(1)
    
    for(i in which(pattern)) { #for each cell of interest

      neighbors = map[x_with_border][x_to_evaluate[i]+interact] #get its neighbors
      if(all(is.na(neighbors)) ) { 
        map[i] = count #then no patch -> = 1
      } else {
        map[i] = min(neighbors, na.rm = TRUE) 
      }
      count = count +1
    }
    
  }
  
  map = as.factor(map)
  patchvec = as.vector(sapply(levels(map), function(i) length(which(map == i) ) )) 
  
  out = vector()
  if(length(patchvec) > 0) out = sort(patchvec) else out = NA
  return(out)
  
} 


mapping = function(width, height, boundary = "periodic", i_matrix = matrix(c(0,1,0,1,NA,1,0,1,0), ncol = 3, byrow = TRUE)) {
  
  X = matrix(as.integer(1:(width*height)), ncol = width, byrow =TRUE)
  X = cbind(X[,width], X, X[,1] )  
  X = rbind(X[height,], X, X[1,] ) 
  x_with_border = as.integer(t(X))
  
  assign("x_with_border", as.integer(t(X))  , envir = .GlobalEnv )
  
  assign("x_to_evaluate", sort(matrix(1:prod(dim(X)), ncol = dim(X)[2], byrow =TRUE)[-c(1, dim(X)[1]), -c(1,dim(X)[2])]  )	, envir = .GlobalEnv )
  
  
  I = i_matrix	

  neighbours_in_I = which(is.finite(abs(I)/abs(I)), arr.in = TRUE)

  relrow = neighbours_in_I[,1]-which(is.na(I), arr.ind = TRUE)[1]
  relcol = neighbours_in_I[,2]-which(is.na(I), arr.ind = TRUE)[2]
  
  assign("interact", relrow * dim(X)[2] + relcol, envir = .GlobalEnv )
  
}

count  = function(x, neighbor) {
  
  neighbors = numeric(length = prod(x$dim))
  x_logical_with_border = (x$cells %in% neighbor)[x_with_border]
  for(k in interact) {
    neighbors = neighbors + x_logical_with_border[x_to_evaluate+k]
  }
  return(neighbors)  
}



Get_frequency_number_patches=function(landscape){
  
  if (is.matrix(landscape)){
    landscape=as.vector(landscape)
  }
  
  d_patch=tibble()
  for (i in c(1:2)){ #for the two species
    d_patch=rbind(d_patch,tibble(patch_size=as.numeric(Get_patches(landscape,i)))%>% add_column(., Species=i))
  }
  d_patch=rbind(d_patch,tibble(patch_size=as.numeric(Get_patches(landscape,"+")))%>% add_column(., Species="+")) #and for global vegetation
  
  d_freq_patch=d_size=tibble()
  for (k in unique(d_patch$Species)){  #for each species and for total vegetation
    size_patches_k=filter(d_patch,Species==k)$patch_size
    cumbins=sort(unique(size_patches_k)) 
    d_freq_patch=rbind(d_freq_patch,tibble(Number=as.numeric(sapply(cumbins, function(k) length(which(size_patches_k >= k)) )),
                       Frequency = as.numeric(sapply(cumbins, function(k) length(which(size_patches_k >= k))/length(size_patches_k) )),
                       Size=cumbins,
                       Species=k))
    
    d_size=rbind(d_size,tibble(Max_size=max(cumbins),
                               Species=k))
    
  }
  
  return(list(Patches_size=d_patch,Patches_frequency=d_freq_patch,Max_size=d_size))

}



# fitting Power Laws



fitpoly =  function(data , indices, modelout = FALSE) {
  model = lm(log(p) ~  - 1 + log(size) + I(log(size)^2), data = data[indices,] )
  if(modelout) {return(model)} else {return(coefficients(model)[2])} 
} 
fitlm =  function(data , indices, modelout = FALSE) {
  model =lm(log(p) ~  - 1 + log(size), data = data[indices,] )
  if(modelout) {return(model)} else {return(coefficients(model)[2])} 
} 



fitPL = function(psd, p_spanning) {
  
  "
  psd = patch size distribution 
  p_spanning = lower limitnàfnthe up-bent power law
  
  "
  
  # code of fitted classes
  

  out = list()
  out$best = NA
  out$AIC = vector("numeric", length = 3)
  out$dAIC = vector("numeric", length = 3)
  
  # criteria for vegetated state & desert state
  
  ##### linear power law model for parameter estimation
  PLlm = lm(I(log(p)) ~  1 - I(log(size)) , data = psd) 
  
  ###########
  
  try( {out$TPLdown = nls(I(log(p)) ~ alpha * log(size) + Sx * (1 - size) , 
                           data = psd,
                           start = list(alpha =  PLlm$coefficients, Sx = 1/1000),
                           #algorithm = "port",
                           trace = FALSE
  )}, silent = TRUE
  )    
  
  if(!is.null(out$TPLdown) & !coefficients(out$TPLdown)["Sx"] <= 0) {
    out$AIC[1] = AIC(out$TPLdown) 
  } else {
    out$TPLdown = list(NA)
    out$AIC[1] = NA
  }
  
  #####
  
  try({out$PL = nls(I(log(p)) ~ alpha * log(size), 
                     data = psd,
                     start = list( alpha =  PLlm$coefficients ),
                     trace = FALSE,
                     nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  if(!is.null(out$PL)) {
    out$AIC[2] = AIC(out$PL)
  } else {
    out$PL  = list(NA)
    out$AIC[2] = NA
  }
  
  ###########
    
    
  try({out$TPLup = nls(I(log(p)) ~  log(b) + log(1+(size^(alpha))/b ) , 
                        data = psd,
                        start = list( alpha =  PLlm$coefficients, b = p_spanning ) , 
                        nls.control(maxiter = 50)
  )}, silent = TRUE
  )
  
  
  if(!is.null(out$TPLup)) {
    out$AIC[3] = AIC(out$TPLup) 
  } else { 
    #result$fit$summary$TPLup  = list(NA)
    out$TPLup  = list(NA)
    out$AIC[3] = NA
  }
  
  ###########
  
  out$dAIC =   out$AIC -min(out$AIC, na.rm = TRUE)
  
  out$best = which.min(out$AIC)+1
  
  return(out)
} 



# E) Temporal analysis ----


Get_temporal_EWS=function(d,burning_phase=3000){
  
  colnames(d)=c("Time","Rho_1","Rho_2","Rho_0","Rho_d")
  d$Rho_tot = d$Rho_1+d$Rho_2
  
  Temp_sd = c(var(d$Rho_1[-(1:burning_phase)]) , var(d$Rho_2[-(1:burning_phase)]),var(d$Rho_tot[-(1:burning_phase)]))
  
  AR1 = c(arima(d$Rho_1[-(1:burning_phase)], order=c(1,0,0))$coef[1],
          arima(d$Rho_2[-(1:burning_phase)], order=c(1,0,0))$coef[1],
          arima(d$Rho_tot[-(1:burning_phase)], order=c(1,0,0))$coef[1])
  
  
  return(tibble(Species=c("1","2","Tot"),Temp_sp=Temp_sd,AR1=AR1))
  
}



# F) Functional diversity ----


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


