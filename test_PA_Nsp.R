source("./Dryland_shift_functions.R")

p=Get_PA_parameters();p["psi_1"]=1;p["psi_2"]=0;p["S"]=.1
p["alpha_0"]=.1

for (i in 1:length(p)) assign(paste0(names(p)[i]),p[i])

u=Get_PA_initial_state()
for (i in 1:length(u)) assign(paste0(names(u)[i]),u[i])


rho_0 = (1 - rho_1 - rho_2 - rho_m)
rho_20 = rho_2 - rho_22 - rho_12 - rho_2m
rho_10 = rho_1 - rho_11 - rho_12 - rho_1m
rho_0m = rho_m - rho_mm - rho_1m - rho_2m

du=rep(0,9)
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


names(du)=names(u)
# for (i in 1:length(du)) assign(paste0(names(du)[i]),du[i])

d_rho_00=-sum(du[1:3]) -(du[2]-du["rho_22"]-du["rho_12"]-du["rho_2m"])-
  (du["rho_1"]-du["rho_11"]-du["rho_12"]-du["rho_1m"])-
  (du["rho_m"]-du["rho_mm"]-du["rho_1m"]-du["rho_2m"]) #d_rho_00

d_rho_00

#d_rho_m0
d_rho_m0 = (du["rho_m"]-du["rho_mm"]-du["rho_1m"]-du["rho_2m"])
d_rho_m0

#d_rho_i0
d_rho_10 = (du["rho_1"]-du["rho_11"]-du["rho_12"]-du["rho_1m"])
d_rho_20 = (du[2]-du["rho_22"]-du["rho_12"]-du["rho_2m"])
d_rho_10;d_rho_20


