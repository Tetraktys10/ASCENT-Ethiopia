#ASCENT model Ethiopia using hmer package

# library(hmer)
# library(deSolve)
# library(ggplot2)
# library(reshape2)
# library(purrr)
# library(tidyverse)
# library(lhs)
# set.seed(123)
# library("writexl")

eventdat <- data.frame(var=c("Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb"),
                       time=c(0,1,2,3,4,5,6,7,8,9,10,0,1,2,3,4,5,6,7,8,9,10),#c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1),
                       value=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       method=c("rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep"))


# Define the helper function `ode_results`, to obtain the solution of the ODEs.
ode_results <- function(parms, end_time = 10) {  #11 years from 2011 to 2021
  forcer1 = matrix(c(0, parms['pi1'], 10, parms['pi2']), 
                   ncol = 2, byrow = TRUE)
  force_func1 = approxfun(x = forcer1[,1], y = forcer1[,2], method = "linear", rule = 2)
  forcer2 = matrix(c(0, parms['mu1'], 10, parms['mu2']), 
                   ncol = 2, byrow = TRUE)
  force_func2 = approxfun(x = forcer2[,1], y = forcer2[,2], method = "linear", rule = 2)
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      #children
      dSc       <-  force_func1(time)*N - (betac*Sc*(Dc+Da))/N - force_func2(time)*Sc - Sc*0.059
      dLec      <-  betac*(Sc+alpha*(Llc+Rnc+Rc))*(Dc+Da)/N - (omega+(0.7*0.29+0.3)*p+force_func2(time))*Lec - Lec*0.059
      dLlc      <-  omega*Lec - ((alpha*betac*(Dc+Da))/N+(0.7*0.29+0.3)*v+force_func2(time))*Llc - Llc*0.059 + 0.5*(Rnc+Rc)
      dDc       <-  (0.7*0.29+0.3)*p*Lec + (0.7*0.29+0.3)*v*Llc + f*Dtc + rho_n*Rnc + rho*Rc - (sigma*q+tau_n+mu_n)*Dc - Dc*0.059
      dDtc      <-  sigma*q*Dc - (f+tau+mut)*Dtc - Dtc*0.059
      dRnc      <-  tau_n*Dc - alpha*betac*Rnc*(Dc+Da)/N - (rho_n+force_func2(time))*Rnc - Rnc*0.059
      dRc       <-  tau*Dtc - (rho+force_func2(time)+0.5)*Rc - Rc*0.059
      
      #adults
      dSa       <-  Sc*0.059 - (betaa*Sa*(Dc+Da))/N - force_func2(time)*Sa
      dLea      <-  Lec*0.059 + betaa*(Sa+alpha*(Lla+Rna+Ra))*(Dc+Da)/N - (omega+p+force_func2(time))*Lea
      dLla      <-  Llc*0.059 + omega*Lea - ((alpha*betaa*(Dc+Da))/N+v+force_func2(time))*Lla + 0.5*(Rna+Ra)
      dDa       <-  Dc*0.059 + p*Lea + v*Lla + f*Dta + rho_n*Rna + rho*Ra - (sigma*q+tau_n+mu_n)*Da
      dDta      <-  Dtc*0.059 + sigma*q*Da - (f+tau+mut)*Dta
      dRna      <-  Rnc*0.059 + tau_n*Da - alpha*betaa*Rna*(Dc+Da)/N - (rho_n+force_func2(time))*Rna
      dRa       <-  Rc*0.059 + tau*Dta - (rho+force_func2(time)+0.5)*Ra
      
      dN        <- dSc+dLec+dLlc+dDc+dDtc+dRc+dRnc+dSa+dLea+dLla+dDa+dDta+dRa+dRna
      dInc      <- (((0.7*0.29+0.3)*p*Lec + (0.7*0.29+0.3)*v*Llc + rho_n*Rnc + rho*Rc + p*Lea + v*Lla + rho_n*Rna + rho*Ra)/N)*1000
      dDtb      <- ((mu_n*Da + mut*Dta + mu_n*Dc + mut*Dtc)/N)*1000   # + (Dc+Dtc+Da+Dta)*mu
      
      return(list(c(dSc,dLec,dLlc,dDc,dDtc,dRc,dRnc,dSa,dLea,dLla,dDa,dDta,dRa,dRna,dN,dInc,dDtb)))
      
    })
  }
  # yini = c(Sc=40242416, Lec=193531, Llc=687913, Dc=8061, Dtc=6979, Rc=10000, Rnc=1200,
  #         Sa=26947265, Lea=2527913, Lla=23076526, Da=64613, Dta=63583, Ra=100000, Rna=1000, N=91818000, Inc=(140000/91818000)*1000, Dtb = (21000/91818000)*1000) #, Inc=159000, Dtb = 30000
  #yini = c(Sc=40242416, Lec=53531, Llc=687913, Dc=8061, Dtc=6979, Rc=10000, Rnc=1200,
  #         Sa=22947265, Lea=2127913, Lla=27076526, Da=64613, Dta=63583, Ra=100000, Rna=1000, N=91818000, Inc=(160000/91818000)*1000, Dtb = (21000/91818000)*1000) #, Inc=159000, Dtb = 30000
  yini = c(Sc=40242416, Lec=53531, Llc=687913, Dc=11061, Dtc=6979, Rc=10000, Rnc=1200,
           Sa=21947265, Lea=1827913, Lla=28076526, Da=168000, Dta=63583, Ra=100000, Rna=1000, N=91818000, Inc=(160000/91818000)*1000, Dtb = (21000/91818000)*1000) #, Inc=159000, Dtb = 30000
  
  
  times = seq(0, end_time, by = 1)
  out = deSolve::ode(yini, times, des, parms, events = list(data = eventdat))
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

# parameters values:[19]
parameters_values <- c(
  pi1 = 0.03580115,pi2 = 0.032398037,#pi3 = 0.0252,#pi = 0.029,#birth 
  mu1 = 0.008503,mu2 = 0.006799,#mu3 = 0.0054,#mu = 0.007,#mortality
  betac = 14.28772075,
  betaa = 15.72088596,
  omega = 0.872002,
  p = 0.0826,
  v = 0.000594,
  tau_n = 0.182631,#0.144, 
  tau = 0.923838,
  mu_n = 0.198601,#0.094,
  mut = 0.039813,#0.04,   #0.019 literature #0.04 facilities data
  sigma = 0.653369,
  q = 0.947627,
  f = 0.012095,
  alpha = 0.220558,
  rho_n = 0.020017,
  rho = 0.010175
)

solution <- ode_results(parameters_values)
par(mar = c(2, 2, 2, 2))
plot(solution)

############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
############################  3. WAVE0 - PARAMETER RANGES, TARGETS AND DESIGN POINTS  ###########################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 

# Define the parameter ranges (19)
ranges = list(
  pi1 = c(0.0357,0.0359),pi2 = c(0.0323,0.0325),#pi = 0.029,#birth 
  mu1 = c(0.0084,0.0086),mu2 = c(0.0067,0.0069),#mu = 0.007,#mortality 
  betac = c(8, 15), # infection rate from time t=0 children
  betaa = c(8, 12), # infection rate from time t=0 adults
  omega = c(0.871, 0.873), ##########################
  p = c(0.0825, 0.0827), ############################
  v = c(0.000593, 0.000595), ########################
  tau_n = c(0.123, 0.247), 
  tau = c(0.911,0.936), #recovery rate standard care
  mu_n = c(0.171, 0.236), 
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
  Dtb1 = c( 0.20, 0.43),
  Dtb2 = c( 0.22, 0.46),
  Dtb3 = c( 0.21, 0.46),
  Dtb4 = c( 0.16, 0.32),
  Dtb5 = c( 0.16, 0.32),
  Dtb6 = c( 0.15, 0.32),
  Dtb7 = c( 0.13, 0.28),
  Dtb8 = c( 0.12, 0.26),
  Inc1 = c(0.751,	2.968),
  Inc2 = c(0.766,	2.782),
  Inc3 = c(0.781,	2.544),
  Inc4 = c(0.952,	2.083),
  Inc5 = c(0.878,	2.016),
  Inc6 = c(0.870,	1.860),
  Inc7 = c(0.855,	1.726),
  Inc8 = c(0.833,	1.599)
)

# Define two Latin hypercube designs through the function `maximinLHS`. This function assumes that each parameter 
# is distributed on [0,1]
initial_LHS_training <- maximinLHS(190, 19)
initial_LHS_validation <- maximinLHS(190, 19)
initial_LHS <- rbind(initial_LHS_training, initial_LHS_validation)

# Rescale the parameter ranges from [0,1] to the correct ranges, and add columns names to identify the parameters
initial_points <- setNames(data.frame(t(apply(initial_LHS, 1, 
                                              function(x) x*unlist(lapply(ranges, function(x) x[2]-x[1])) + 
                                                unlist(lapply(ranges, function(x) x[1]))))), names(ranges))

# Run the model on `initial_points` and add column names to identify the different targets 
initial_results <- setNames(data.frame(t(apply(initial_points, 1, get_results, 
                                               c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))

# Bind `initial_points` and the corresponding model outputs `initial_results` by column
wave0 <- cbind(initial_points, initial_results)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 
#################################################  4. EMULATORS  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

# Split the dataframe `wave0` into a training and a validation set
training <- wave0[1:190,]
validation <- wave0[191:380,]


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


new_points <- generate_new_runs(ems_wave1, 380, targets, verbose = TRUE)
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
                                                   c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))

# Bind by columns `new_points` to the model output `new_initial_results` to create the data.frame `wave1`
wave1 <- cbind(new_points, new_initial_results)

# Split `wave1` into training and validation sets
new_t_sample <- sample(1:nrow(wave1), 190)
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

new_new_points <- generate_new_runs(c(ems_wave2, ems_wave1), 380, targets, verbose=TRUE)
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

new_new_initial_results <- setNames(data.frame(t(apply(new_new_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave2 <- cbind(new_new_points, new_new_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  THIRD WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new2_t_sample <- sample(1:nrow(wave2), 190)
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
ems_wave3 <- ems_wave3[!seq_along(ems_wave3) %in% bad.ems]
vd <- validation_diagnostics(ems_wave3, validation = new2_validation, targets = targets, plt=TRUE)

new3_points <- generate_new_runs(c(ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)
#plot_wrap(new_new_points, ranges)
R_squared_new3 <- list()
for (i in 1:length(ems_wave3)) {
  R_squared_new3[[i]] <- summary(ems_wave3[[i]]$model)$adj.r.squared
}
names(R_squared_new3) <- names(ems_wave3)
unlist(R_squared_new3)

# Create a dataframe `wave3` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new3_initial_results <- setNames(data.frame(t(apply(new3_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))

wave3 <- cbind(new3_points, new3_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3)
simulator_plot(all_points, targets)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  FOURTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new3_t_sample <- sample(1:nrow(wave3), 190)
new3_training <- wave3[new3_t_sample,]
new3_validation <- wave3[-new3_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave4 <- emulator_from_data(new3_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
#Validate emulators iteratively
for (j in 1:length(ems_wave4)) {
  misclass <- nrow(classification_diag(ems_wave4[[j]], targets, new3_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave4[[j]] <- ems_wave4[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave4[[j]], targets, new3_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave4)) {
  bad.model <- nrow(comparison_diag(ems_wave4[[j]], targets, new3_validation, plt = FALSE))
  if (bad.model > floor(nrow(new3_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave4 <- ems_wave4[!seq_along(ems_wave4) %in% bad.ems]

vd <- validation_diagnostics(ems_wave4, validation = new3_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new4_points <- generate_new_runs(c(ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)
#plot_wrap(new4_points, ranges)

R_squared_new4 <- list()
for (i in 1:length(ems_wave4)) {
  R_squared_new4[[i]] <- summary(ems_wave4[[i]]$model)$adj.r.squared
}
names(R_squared_new4) <- names(ems_wave4)
unlist(R_squared_new4)

# Create a dataframe `wave4` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new4_initial_results <- setNames(data.frame(t(apply(new4_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave4 <- cbind(new4_points, new4_initial_results)



# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4)
simulator_plot(all_points, targets)

all_points <- list(wave3, wave4)
simulator_plot(all_points, targets)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  FIFTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split the new wave into training and validation sets
new4_t_sample <- sample(1:nrow(wave4), 190)
new4_training <- wave4[new4_t_sample,]
new4_validation <- wave4[-new4_t_sample,]

# Train the new wave emulators using `emulator_from_data`, passing the new ranges to it
ems_wave5 <- emulator_from_data(new4_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all emulators
for (j in 1:length(ems_wave5)) {
  misclass <- nrow(classification_diag(ems_wave5[[j]], targets, new4_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave5[[j]] <- ems_wave5[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave5[[j]], targets, new4_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave5)) {
  bad.model <- nrow(comparison_diag(ems_wave5[[j]], targets, new4_validation, plt = FALSE))
  if (bad.model > floor(nrow(new4_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave5 <- ems_wave5[!seq_along(ems_wave5) %in% bad.ems]

vd <- validation_diagnostics(ems_wave5, validation = new4_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new5_points <- generate_new_runs(c(ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new5 <- list()
for (i in 1:length(ems_wave5)) {
  R_squared_new5[[i]] <- summary(ems_wave5[[i]]$model)$adj.r.squared
}
names(R_squared_new5) <- names(ems_wave5)
unlist(R_squared_new5)

# Create a dataframe `wave5` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new5_initial_results <- setNames(data.frame(t(apply(new5_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave5 <- cbind(new5_points, new5_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5)
simulator_plot(all_points, targets)

all_points <- list(wave4, wave5)
simulator_plot(all_points, targets)

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  SIXTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new5_t_sample <- sample(1:nrow(wave5), 190)
new5_training <- wave5[new5_t_sample,]
new5_validation <- wave5[-new5_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave6 <- emulator_from_data(new5_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave6)) {
  misclass <- nrow(classification_diag(ems_wave6[[j]], targets, new5_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave6[[j]] <- ems_wave6[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave6[[j]], targets, new5_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave6)) {
  bad.model <- nrow(comparison_diag(ems_wave6[[j]], targets, new5_validation, plt = FALSE))
  if (bad.model > floor(nrow(new5_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave6 <- ems_wave6[!seq_along(ems_wave6) %in% bad.ems]

vd <- validation_diagnostics(ems_wave6, validation = new5_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new6_points <- generate_new_runs(c(ems_wave6, ems_wave5,ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new6 <- list()
for (i in 1:length(ems_wave6)) {
  R_squared_new6[[i]] <- summary(ems_wave6[[i]]$model)$adj.r.squared
}
names(R_squared_new6) <- names(ems_wave6)
unlist(R_squared_new6)

# Create a dataframe `wave5` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new6_initial_results <- setNames(data.frame(t(apply(new6_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave6 <- cbind(new6_points, new6_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6)
simulator_plot(all_points, targets)

all_points <- list(wave5, wave6)
simulator_plot(all_points, targets)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  SEVENTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new6_t_sample <- sample(1:nrow(wave6), 190)
new6_training <- wave6[new6_t_sample,]
new6_validation <- wave6[-new6_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave7 <- emulator_from_data(new6_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave7)) {
  misclass <- nrow(classification_diag(ems_wave7[[j]], targets, new6_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave7[[j]] <- ems_wave7[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave7[[j]], targets, new6_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave7)) {
  bad.model <- nrow(comparison_diag(ems_wave7[[j]], targets, new6_validation, plt = FALSE))
  if (bad.model > floor(nrow(new6_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave7 <- ems_wave7[!seq_along(ems_wave7) %in% bad.ems]

vd <- validation_diagnostics(ems_wave7, validation = new6_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new7_points <- generate_new_runs(c(ems_wave7, ems_wave6, ems_wave5,ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new7 <- list()
for (i in 1:length(ems_wave7)) {
  R_squared_new7[[i]] <- summary(ems_wave7[[i]]$model)$adj.r.squared
}
names(R_squared_new7) <- names(ems_wave7)
unlist(R_squared_new7)

# Create a dataframe `wave5` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new7_initial_results <- setNames(data.frame(t(apply(new7_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave7 <- cbind(new7_points, new7_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7)
simulator_plot(all_points, targets)

all_points <- list(wave6, wave7)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  EIGTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new7_t_sample <- sample(1:nrow(wave7), 190)
new7_training <- wave7[new7_t_sample,]
new7_validation <- wave7[-new7_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave8 <- emulator_from_data(new7_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave8)) {
  misclass <- nrow(classification_diag(ems_wave8[[j]], targets, new7_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave8[[j]] <- ems_wave8[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave8[[j]], targets, new7_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave8)) {
  bad.model <- nrow(comparison_diag(ems_wave8[[j]], targets, new7_validation, plt = FALSE))
  if (bad.model > floor(nrow(new7_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave8 <- ems_wave8[!seq_along(ems_wave8) %in% bad.ems]

vd <- validation_diagnostics(ems_wave8, validation = new7_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new8_points <- generate_new_runs(c(ems_wave8, ems_wave7, ems_wave6, ems_wave5,ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new8 <- list()
for (i in 1:length(ems_wave8)) {
  R_squared_new8[[i]] <- summary(ems_wave8[[i]]$model)$adj.r.squared
}
names(R_squared_new8) <- names(ems_wave8)
unlist(R_squared_new8)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new8_initial_results <- setNames(data.frame(t(apply(new8_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave8 <- cbind(new8_points, new8_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8)
simulator_plot(all_points, targets)

all_points <- list(wave7, wave8)
simulator_plot(all_points, targets)



# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  NINTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new8_t_sample <- sample(1:nrow(wave8), 190)
new8_training <- wave8[new8_t_sample,]
new8_validation <- wave8[-new8_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave9 <- emulator_from_data(new8_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave9)) {
  misclass <- nrow(classification_diag(ems_wave9[[j]], targets, new8_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave9[[j]] <- ems_wave9[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave9[[j]], targets, new8_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave9)) {
  bad.model <- nrow(comparison_diag(ems_wave9[[j]], targets, new8_validation, plt = FALSE))
  if (bad.model > floor(nrow(new8_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave9 <- ems_wave9[!seq_along(ems_wave9) %in% bad.ems]

vd <- validation_diagnostics(ems_wave9, validation = new8_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new9_points <- generate_new_runs(c(ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5,ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new9 <- list()
for (i in 1:length(ems_wave9)) {
  R_squared_new9[[i]] <- summary(ems_wave9[[i]]$model)$adj.r.squared
}
names(R_squared_new9) <- names(ems_wave9)
unlist(R_squared_new9)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new9_initial_results <- setNames(data.frame(t(apply(new9_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave9 <- cbind(new9_points, new9_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8, wave9)
simulator_plot(all_points, targets)

all_points <- list(wave8, wave9)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  TENTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new9_t_sample <- sample(1:nrow(wave9), 190)
new9_training <- wave9[new9_t_sample,]
new9_validation <- wave9[-new9_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave10 <- emulator_from_data(new9_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave10)) {
  misclass <- nrow(classification_diag(ems_wave10[[j]], targets, new9_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave10[[j]] <- ems_wave10[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave10[[j]], targets, new9_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave10)) {
  bad.model <- nrow(comparison_diag(ems_wave10[[j]], targets, new9_validation, plt = FALSE))
  if (bad.model > floor(nrow(new9_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave10 <- ems_wave10[!seq_along(ems_wave10) %in% bad.ems]

vd <- validation_diagnostics(ems_wave10, validation = new9_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new10_points <- generate_new_runs(c(ems_wave10, ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new10 <- list()
for (i in 1:length(ems_wave10)) {
  R_squared_new10[[i]] <- summary(ems_wave10[[i]]$model)$adj.r.squared
}
names(R_squared_new10) <- names(ems_wave10)
unlist(R_squared_new10)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new10_initial_results <- setNames(data.frame(t(apply(new10_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave10 <- cbind(new10_points, new10_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8, wave9, wave10)
simulator_plot(all_points, targets)

all_points <- list(wave9, wave10)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  ELEVENTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new10_t_sample <- sample(1:nrow(wave10), 190)
new10_training <- wave10[new10_t_sample,]
new10_validation <- wave10[-new10_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave11 <- emulator_from_data(new10_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave11)) {
  misclass <- nrow(classification_diag(ems_wave11[[j]], targets, new10_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave11[[j]] <- ems_wave11[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave11[[j]], targets, new10_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave11)) {
  bad.model <- nrow(comparison_diag(ems_wave11[[j]], targets, new10_validation, plt = FALSE))
  if (bad.model > floor(nrow(new10_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave11 <- ems_wave11[!seq_along(ems_wave11) %in% bad.ems]

vd <- validation_diagnostics(ems_wave11, validation = new10_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new11_points <- generate_new_runs(c(ems_wave11, ems_wave10, ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new11 <- list()
for (i in 1:length(ems_wave11)) {
  R_squared_new11[[i]] <- summary(ems_wave11[[i]]$model)$adj.r.squared
}
names(R_squared_new11) <- names(ems_wave11)
unlist(R_squared_new11)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new11_initial_results <- setNames(data.frame(t(apply(new11_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave11 <- cbind(new11_points, new11_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8, wave9, wave10, wave11)
simulator_plot(all_points, targets)

all_points <- list(wave10, wave11)
simulator_plot(all_points, targets)
#write_xlsx(wave11,"N:\\Documents\\wave11_v3.xlsx")

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  TWELVETH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new11_t_sample <- sample(1:nrow(wave11), 190)
new11_training <- wave11[new11_t_sample,]
new11_validation <- wave11[-new11_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave12 <- emulator_from_data(new11_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave12)) {
  misclass <- nrow(classification_diag(ems_wave12[[j]], targets, new11_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave12[[j]] <- ems_wave12[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave12[[j]], targets, new11_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave12)) {
  bad.model <- nrow(comparison_diag(ems_wave12[[j]], targets, new11_validation, plt = FALSE))
  if (bad.model > floor(nrow(new11_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave12 <- ems_wave12[!seq_along(ems_wave12) %in% bad.ems]

vd <- validation_diagnostics(ems_wave12, validation = new11_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new12_points <- generate_new_runs(c(ems_wave12, ems_wave11, ems_wave10, ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new12 <- list()
for (i in 1:length(ems_wave12)) {
  R_squared_new12[[i]] <- summary(ems_wave12[[i]]$model)$adj.r.squared
}
names(R_squared_new12) <- names(ems_wave12)
unlist(R_squared_new12)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new12_initial_results <- setNames(data.frame(t(apply(new12_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave12 <- cbind(new12_points, new12_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8, wave9, wave10, wave11, wave12)
simulator_plot(all_points, targets)

all_points <- list(wave11, wave12)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  THIRTEENTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new12_t_sample <- sample(1:nrow(wave12), 190)
new12_training <- wave12[new12_t_sample,]
new12_validation <- wave12[-new12_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave13 <- emulator_from_data(new12_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave13)) {
  misclass <- nrow(classification_diag(ems_wave13[[j]], targets, new12_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave13[[j]] <- ems_wave13[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave13[[j]], targets, new12_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave13)) {
  bad.model <- nrow(comparison_diag(ems_wave13[[j]], targets, new12_validation, plt = FALSE))
  if (bad.model > floor(nrow(new12_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave13 <- ems_wave13[!seq_along(ems_wave13) %in% bad.ems]

vd <- validation_diagnostics(ems_wave13, validation = new12_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new13_points <- generate_new_runs(c(ems_wave13, ems_wave12, ems_wave11, ems_wave10, ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new13 <- list()
for (i in 1:length(ems_wave13)) {
  R_squared_new13[[i]] <- summary(ems_wave13[[i]]$model)$adj.r.squared
}
names(R_squared_new13) <- names(ems_wave13)
unlist(R_squared_new13)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new13_initial_results <- setNames(data.frame(t(apply(new13_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave13 <- cbind(new13_points, new13_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8, wave9, wave10, wave11, wave12, wave13)
simulator_plot(all_points, targets)

all_points <- list(wave12, wave13)
simulator_plot(all_points, targets)


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
###############################################  FOURTEENTH WAVE  ##############################################
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #


# Split `wave1` into training and validation sets
new13_t_sample <- sample(1:nrow(wave13), 190)
new13_training <- wave13[new13_t_sample,]
new13_validation <- wave13[-new13_t_sample,]

# Train wave 2 emulators using `emulator_from_data`, passing the new ranges to it
ems_wave14 <- emulator_from_data(new13_training, names(targets), ranges, check.ranges = TRUE)
# Produce diagnostics for all wave 2 emulators
for (j in 1:length(ems_wave14)) {
  misclass <- nrow(classification_diag(ems_wave14[[j]], targets, new13_validation, plt = FALSE))
  while(misclass > 0) {
    ems_wave14[[j]] <- ems_wave13[[j]]$mult_sigma(1.1)
    misclass <- nrow(classification_diag(ems_wave14[[j]], targets, new13_validation, plt = FALSE))
  }
}

bad.ems <- c()
for (j in 1:length(ems_wave14)) {
  bad.model <- nrow(comparison_diag(ems_wave14[[j]], targets, new13_validation, plt = FALSE))
  if (bad.model > floor(nrow(new13_validation)/4)) {
    bad.ems <- c(bad.ems, j)
  }
}
ems_wave14 <- ems_wave14[!seq_along(ems_wave14) %in% bad.ems]

vd <- validation_diagnostics(ems_wave14, validation = new13_validation, targets = targets, plt=TRUE)

# Generate 380 new parameter sets and plot them
new14_points <- generate_new_runs(c(ems_wave14, ems_wave13, ems_wave12, ems_wave11, ems_wave10, ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 380, targets, verbose=TRUE)

R_squared_new14 <- list()
for (i in 1:length(ems_wave14)) {
  R_squared_new14[[i]] <- summary(ems_wave14[[i]]$model)$adj.r.squared
}
names(R_squared_new14) <- names(ems_wave14)
unlist(R_squared_new14)

# Create a dataframe `wave8` binding the parameters sets generated at the end of wave 2 with the corresponding 
# model outputs

new14_initial_results <- setNames(data.frame(t(apply(new14_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


wave14 <- cbind(new14_points, new14_initial_results)

# Assess how much better parameter sets at later waves perform compared to the original `initial_points` through 
#`simulator_plot`
all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, wave6, wave7, wave8, wave9, wave10, wave11, wave12, wave13, wave14)
simulator_plot(all_points, targets)

all_points <- list(wave13, wave14)
simulator_plot(all_points, targets)


#########
final_points <- generate_new_runs(c(ems_wave14, ems_wave13, ems_wave12, ems_wave11, ems_wave10, ems_wave9, ems_wave8, ems_wave7, ems_wave6, ems_wave5, ems_wave4, ems_wave3, ems_wave2, ems_wave1), 1000, targets, verbose=TRUE)

final_results <- setNames(data.frame(t(apply(final_points, 1, get_results,c(1,2,3,4,5,6,7,8), c('Dtb','Inc')))), names(targets))


final_wave <- cbind(final_points, final_results)

write_xlsx(final_wave,"N:\\Documents\\final_wave.xlsx")

all_points <- list(final_wave)
simulator_plot(all_points, targets)

#########
final_r <- data.frame()
row_to_keep <- 1:1000
for (i in 1:1000) {
  if (final_wave[i,20]>=targets$Dtb1[1]&&final_wave[i,20]<=targets$Dtb1[2]&&
      final_wave[i,21]>=targets$Dtb2[1]&&final_wave[i,21]<=targets$Dtb2[2]&&
      final_wave[i,22]>=targets$Dtb3[1]&&final_wave[i,22]<=targets$Dtb3[2]&&
      final_wave[i,23]>=targets$Dtb4[1]&&final_wave[i,23]<=targets$Dtb4[2]&&
      final_wave[i,24]>=targets$Dtb5[1]&&final_wave[i,24]<=targets$Dtb5[2]&&
      final_wave[i,25]>=targets$Dtb6[1]&&final_wave[i,25]<=targets$Dtb6[2]&&
      final_wave[i,26]>=targets$Dtb7[1]&&final_wave[i,26]<=targets$Dtb7[2]&&
      final_wave[i,27]>=targets$Dtb8[1]&&final_wave[i,27]<=targets$Dtb8[2]&&
      final_wave[i,28]>=targets$Inc1[1]&&final_wave[i,28]<=targets$Inc1[2]&&
      final_wave[i,29]>=targets$Inc2[1]&&final_wave[i,29]<=targets$Inc2[2]&&
      final_wave[i,30]>=targets$Inc3[1]&&final_wave[i,30]<=targets$Inc3[2]&&
      final_wave[i,31]>=targets$Inc4[1]&&final_wave[i,31]<=targets$Inc4[2]&&
      final_wave[i,32]>=targets$Inc5[1]&&final_wave[i,32]<=targets$Inc5[2]&&
      final_wave[i,33]>=targets$Inc6[1]&&final_wave[i,33]<=targets$Inc6[2]&&
      final_wave[i,34]>=targets$Inc7[1]&&final_wave[i,34]<=targets$Inc7[2]&&
      final_wave[i,35]>=targets$Inc8[1]&&final_wave[i,35]<=targets$Inc8[2]) {
    #final_r[i,] <- final_wave[i,]#rbind(final_r, final_wave[i,]) #final_wave[i,]#
    # new_element <- final_wave[i,]
    # final_r[[length(final_r) + 1]] <- new_element
    row_to_keep[i] = TRUE
  } else {
    row_to_keep[i] = FALSE
  }
}
final_r <- final_wave
keep <- which(row_to_keep==TRUE)
final_r2 <- final_r[keep,]
write_xlsx(final_r2,"N:\\Documents\\final_results.xlsx")

all_points <- list(final_r2)
simulator_plot(all_points, targets)

all_points <- list(wave0, wave1, wave2, wave3, wave4, wave5, final_r2) #wave6, wave7, wave8, wave9, wave10, wave11, wave12, wave13, 
simulator_plot(all_points, targets)