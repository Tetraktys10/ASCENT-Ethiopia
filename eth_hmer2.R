#ASCENT model Ethiopia using hmer package

library(hmer)
library(deSolve)
library(ggplot2)
library(reshape2)
library(purrr)
library(tidyverse)
library(lhs)
set.seed(123)

# Define the helper function `ode_results`, to obtain the solution of the ODEs.
ode_results <- function(parms, end_time = 10) {  #11 years from 2011 to 2021
  forcer = matrix(c(0, parms['beta1'], 3285, parms['beta2']), #two values of beta 2011-2019 and 2020-2021 to take into account the possible effect of covid restrictions on TB transmission
                  ncol = 2, byrow = TRUE)
  force_func = approxfun(x = forcer[,1], y = forcer[,2], method = "linear", rule = 2)
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      #children
      dSc       <-  pi - (force_func(time)*Sc*(Dc+Da))/N - (mu*1000/N)*Sc - Sc*0.066    
      dLec      <-  force_func(time)*(Sc+alpha*(Llc+Rnc+Rc))*(Dc+Da)/N - (omega+(0.7*0.47+0.3)*p+mu*1000/N)*Lec - Lec*0.066
      dLlc      <-  omega*Lec - ((alpha*force_func(time)*(Dc+Da))/N+(0.7*0.47+0.3)*v+mu*1000/N)*Llc - Llc*0.066
      dDc       <-  (0.7*0.47+0.3)*p*Lec + (0.7*0.47+0.3)*v*Llc + f*Dtc + rho_n*Rnc + rho*Rc - (sigma*q+tau_n+mu_n)*Dc - Dc*0.066 -(mu*1000/N)*Dc
      #dDt_starc <-  sigma*g*q*Dc - (f_star+tau_star+mut_star)*Dt_starc - Dt_starc*0.066
      dDtc      <-  sigma*q*Dc - (f+tau+mut)*Dtc - Dtc*0.066 -(mu*1000/N)*Dtc
      dRnc      <-  tau_n*Dc - alpha*force_func(time)*Rnc*(Dc+Da)/N - (rho_n+mu*1000/N)*Rnc - Rnc*0.066
      #dR_starc  <-  tau_star*Dt_starc - alpha*force_func(time)*R_starc*(Dc+Da)/N - (rho_star+mu*1000/N)*R_starc - R_starc*0.066
      dRc       <-  tau*Dtc - alpha*force_func(time)*Rc*(Dc+Da)/N - (rho+mu*1000/N)*Rc - Rc*0.066
      
      #adults
      dSa       <-  Sc*0.066 - (force_func(time)*Sa*(Dc+Da))/N - (mu*1000/N)*Sa
      dLea      <-  Lec*0.066 + force_func(time)*(Sa+alpha*(Lla+Rna+Ra))*(Dc+Da)/N - (omega+p+mu*1000/N)*Lea
      dLla      <-  Llc*0.066 + omega*Lea - ((alpha*force_func(time)*(Dc+Da))/N+v+mu*1000/N)*Lla
      dDa       <-  Dc*0.066 + p*Lea + v*Lla + f*Dta + rho_n*Rna + rho*Ra - (sigma*q+tau_n+mu_n)*Da -(mu*1000/N)*Da
      #dDt_stara <-  Dt_starc*0.066 + sigma*g*q*Da - (f_star+tau_star+mut_star)*Dt_stara
      dDta      <-  Dtc*0.066 + sigma*q*Da - (f+tau+mut)*Dta -(mu*1000/N)*Dta
      dRna      <-  Rnc*0.066 + tau_n*Da - alpha*force_func(time)*Rna*(Dc+Da)/N - (rho_n+mu*1000/N)*Rna
      #dR_stara  <-  R_starc*0.066 + tau_star*Dt_stara - alpha*force_func(time)*R_stara*(Dc+Da)/N - (rho_star+mu*1000/N)*R_stara
      dRa       <-  Rc*0.066 + tau*Dta - alpha*force_func(time)*Ra*(Dc+Da)/N - (rho+mu*1000/N)*Ra
      
      dN        <- Sc+Lec+Llc+Dc+Dtc+Rc+Rnc+Sa+Lea+Lla+Da+Dta+Ra+Rna
      dDead     <- (mu_n+mu*1000/N)*Dc+(mut+mu*1000/N)*Dtc+(mu_n+mu*1000/N)*Da+(mut+mu*1000/N)*Dta-Dead  #+mu*1000/N*Rc+mu*1000/N*Rnc+mu*1000/N*Ra+mu*1000/N*Rna
      dInc      <- (0.7*0.47+0.3)*p*Lec + (0.7*0.47+0.3)*v*Llc + rho_n*Rnc + rho*Rc + p*Lea + v*Lla + rho_n*Rna + rho*Ra-Inc
      
      return(list(c(dSc,dLec,dLlc,dDc,dDtc,dRc,dRnc,dSa,dLea,dLla,dDa,dDta,dRa,dRna,dN,dDead,dInc)))
      
    })
  }
  yini = c(Sc=38655496, Lec=82510, Llc=2147490, Dc=4768, Dtc=33936, Rc=30671, Rnc=129,  
           Sa=26110378, Lea=867000, Lla=23400000, Da=33378, Dta=237549, Ra=901, Rna=213794, N=91818000, Dead=56000, Inc=378000)
  times = seq(0, end_time, by = 1)
  out = deSolve::ode(yini, times, des, parms)
  return(out)
}

# Define the helper function `get_results` that acts as `ode_results`, but has the additional feature
# of allowing us to decide which outputs and times should be returned
get_results <- function(params, times, outputs) {
  t_max <- max(times)
  all_res <- ode_results(params, t_max)
  actual_res <- all_res[all_res[,'time'] %in% times, c('time', outputs)]
  shaped <- reshape2::melt(actual_res[,outputs])
  return(setNames(shaped$value, paste0(shaped$Var2, actual_res[,'time'], sep = "")))
}

# Define the helper function `plot_runs` which takes a data frame containing parameter sets and plots the trajectories for
# S, E, I, R for them.
plot_runs <- function(points){
  sol <- list()
  for (i in 1:nrow(points)) sol[[i]]<-ode_results(points[i,])
  par(mar = c(3, 3, 3, 3))
  do.call(plot, sol)
}

# parameters values:[17]
parameters_values <- c(
  pi = 30.9*1000,
  mu = 6.5,
  beta1 = 75, beta2 = 75,
  omega = 0.872,
  p = 0.0826,
  v = 0.000594,
  tau_n = 0.144, 
  tau = 0.923,
  mu_n = 0.077,
  mut = 0.04,
  sigma = 0.71,
  q = 1,
  f = 0.012,
  alpha = 0.21,
  rho_n = 0.02,
  rho = 0.01
)

solution <- ode_results(parameters_values)
par(mar = c(2, 2, 2, 2))
plot(solution)

############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################  3. WAVE0 - PARAMETER RANGES, TARGETS AND DESIGN POINTS  ###########################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Define the parameter ranges 
ranges = list(
  pi = c(22000, 36500), # birth rate
  mu = c(5.4, 9), # rate of death from other causes
  beta1 = c(0, 100), # infection rate from time t=0
  beta2 = c(0, 100), # infection rate from time t=9
  omega = c(0.871, 0.873), ##########################
  p = c(0.0825, 0.0827), ############################
  v = c(0.000593, 0.000595), ########################
  tau_n = c(0.088, 0.220), 
  tau = c(0.911,0.936), #recovery rate standard care
  mu_n = c(0.062, 0.094), 
  mut = c(0.031,0.048),
  sigma = c(0.53,1),
  q = c(0.9,1), #######################################
  f = c(0.008,0.016),
  alpha = c(0.14,0.30),
  rho_n = c(0.015,0.025),################################
  rho = c(0,0.02) ####################################
)

#uncertainty around targets
targets <- list(
  Dead1 = c(27000,	58000),
  Dead2 = c(29000,	62000),
  Dead3 = c(28000,	62000),
  Dead4 = c(22000,	43000),
  Dead5 = c(21000,	43000),
  Dead6 = c(20000,	43000),
  Dead7 = c(18000,	38000),
  Dead8 = c(16000,	35000),
  Dead9 = c(15000,	33000),
  Dead10 = c(13000,	30000),
  Inc1 = c(101000,	399000),
  Inc2 = c(103000,	374000),
  Inc3 = c(105000,	342000),
  Inc4 = c(128000,	280000),
  Inc5 = c(118000,	271000),
  Inc6 = c(117000,	250000),
  Inc7 = c(115000,	232000),
  Inc8 = c(112000,	215000),
  Inc9 = c(105000,	206000),
  Inc10 = c(96000,	199000)
)

# Define two Latin hypercube designs through the function `maximinLHS`. This function assumes that each parameter 
# is distributed on [0,1]
initial_LHS_training <- maximinLHS(170, 17)
initial_LHS_validation <- maximinLHS(170, 17)
initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)

# Rescale the parameter ranges from [0,1] to the correct ranges, and add columns names to identify the parameters
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                                              function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                                                unlist(lapply(ranges, function(x) x[1]))))), names(ranges))

# Run the model on `initial_points` and add column names to identify the different targets 
initial_results <- setNames(data.frame(t(apply(initial_points, 1, get_results, 
                                               c(1,2,3,4,5,6,7,8,9,10), c('Dead','Inc')))), names(targets))

# Bind `initial_points` and the corresponding model outputs `initial_results` by column
wave0 <- cbind(initial_points, initial_results)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#################################################  4. EMULATORS  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Split the dataframe `wave0` into a training and a validation set
training <- wave0[1:170,]
validation <- wave0[171:340,]


# # # # # # # # # # # # # # # # # # # #  4.2 TRAINING EMULATORS  # # # # # # # # # # # # # # # # # # # # # # #

# Train the first set of emulators using the function `emulator_from_data`
ems_wave1 <- emulator_from_data(training, names(targets), ranges)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###########################################  6. EMULATOR DIAGNOSTICS  #########################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

#vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets, plt=TRUE)

#Validate emulators iteratively
for (j in 1:length(ems_wave1)) {
  misclass <- nrow(classification_diag(ems_wave1[[j]], targets, validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave1[[j]] <- ems_wave1[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave1[[j]], targets, validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave1)) {
  bad.model <- nrow(comparison_diag(ems_wave1[[j]], targets, validation, plt = FALSE))
  if (bad.model > floor(nrow(validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave1 <- ems_wave1[!seq_along(ems_wave1) %in% bad.ems]
vd <- validation_diagnostics(ems_wave1, validation = validation, targets = targets, plt=TRUE)


new_points <- generate_new_runs(ems_wave1, 340, targets, verbose = TRUE)
# plot_wrap(new_points, ranges)

R_squared1 <- list()
for (i in 1:length(ems_wave1)) {
  R_squared1[[i]] <- summary(ems_wave1[[i]]$model)$adj.r.squared
}
names(R_squared1) <- names(ems_wave1)
unlist(R_squared1)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  9. SECOND WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Start by evaluating the function `get_results` on `new_points` 
new_initial_results <- setNames(data.frame(t(apply(new_points, 1, get_results,
                                                   c(1,2,3,4,5,6,7,8,9,10), c('Dead','Inc')))), names(targets))

# Bind by columns `new_points` to the model output `new_initial_results` to create the data.frame `wave1`
wave1 <- cbind(new_points, new_initial_results)

# Split `wave1` into training and validation sets
new_t_sample <- sample(1:nrow(wave1), 170)
new_training <- wave1[new_t_sample,]
new_validation <- wave1[-new_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave2 <- emulator_from_data(new_training, names(targets), ranges, check.ranges = TRUE)

#Validate emulators iteratively
for (j in 1:length(ems_wave2)) {
  misclass <- nrow(classification_diag(ems_wave2[[j]], targets, new_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave2[[j]] <- ems_wave2[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave2[[j]], targets, new_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave2)) {
  bad.model <- nrow(comparison_diag(ems_wave2[[j]], targets, new_validation, plt = FALSE))
  if (bad.model > floor(nrow(new_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave2 <- ems_wave2[!seq_along(ems_wave2) %in% bad.ems]
vd <- validation_diagnostics(ems_wave2, validation = new_validation, targets = targets, plt=TRUE)

new_new_points <- generate_new_runs(c(ems_wave2, ems_wave1), 340, targets, verbose=TRUE)
#plot_wrap(new_new_points, ranges)

R_squared_new <- list()
for (i in 1:length(ems_wave2)) {
  R_squared_new[[i]] <- summary(ems_wave2[[i]]$model)$adj.r.squared
}
names(R_squared_new) <- names(ems_wave2)
unlist(R_squared_new)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#############################  10. VISUALISATIONS OF NON-IMPLAUSIBLE SPACE BY WAVE  ###########################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Show the distribution of the non-implausible space before the wave 1, at the end of wave 1 and at the end of
# wave 2 using the function `wave_points`
#wave_points(list(initial_points, new_points, new_new_points), input_names = names(ranges), p_size = 1)

# Create a dataframe `wave2` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new_new_initial_results <- setNames(data.frame(t(apply(new_new_points, 1, get_results,c(1,2,3,4,5,6,7,8,9,10), c('Dead','Inc')))), names(targets))


wave2 <- cbind(new_new_points, new_new_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  THIRD WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new2_t_sample <- sample(1:nrow(wave2), 170)
new2_training <- wave2[new2_t_sample,]
new2_validation <- wave2[-new2_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave3 <- emulator_from_data(new2_training, names(targets), ranges, check.ranges = TRUE)

#Validate emulators iteratively
for (j in 1:length(ems_wave3)) {
  misclass <- nrow(classification_diag(ems_wave3[[j]], targets, new2_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave3[[j]] <- ems_wave3[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave3[[j]], targets, new2_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave3)) {
  bad.model <- nrow(comparison_diag(ems_wave3[[j]], targets, new2_validation, plt = FALSE))
  if (bad.model > floor(nrow(new2_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave3 <- ems_wave3[!seq_along(ems_wave2) %in% bad.ems]
vd <- validation_diagnostics(ems_wave3, validation = new2_validation, targets = targets, plt=TRUE)

new3_points <- generate_new_runs(c(ems_wave3, ems_wave2), 340, targets, verbose=TRUE)
#plot_wrap(new_new_points, ranges)


