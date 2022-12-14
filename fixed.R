# library(hmer)
# library(deSolve)
# library(ggplot2)
# library(reshape2)
# library(purrr)
# library(tidyverse)
# library(lhs)
# set.seed(123)

# Define the helper function `ode_results`, to obtain the solution of the ODEs.
ode_results <- function(parms, end_time = 10) {  #11 years from 2011 to 2021
  forcer3 = matrix(c(0, parms['betac1'], 10, parms['betac2']), 
                   ncol = 2, byrow = TRUE)
  force_func3 = approxfun(x = forcer3[,1], y = forcer3[,2], method = "linear", rule = 2)
  forcer4 = matrix(c(0, parms['betaa1'], 10, parms['betaa2']), 
                   ncol = 2, byrow = TRUE)
  force_func4 = approxfun(x = forcer4[,1], y = forcer4[,2], method = "linear", rule = 2)
  des = function(time, state, parms) {
    with(as.list(c(state, parms)), {
      #children
      dSc       <-  pi*N - (force_func3(time)*Sc*(Dc+Da))/N - mu*Sc - Sc*0.059
      dLec      <-  force_func3(time)*(Sc+alpha*(Llc+Rnc+R_starc+Rc))*(Dc+Da)/N - (omega+(0.7*0.47+0.3)*p+mu)*Lec - Lec*0.059
      dLlc      <-  omega*Lec - ((alpha*force_func3(time)*(Dc+Da))/N+(0.7*0.47+0.3)*v+mu)*Llc - Llc*0.059
      dDc       <-  (0.7*0.47+0.3)*p*Lec + (0.7*0.47+0.3)*v*Llc + f*Dt_starc + f*Dtc + rho_n*Rnc + rho*Rc + rho_star*R_starc - (sigma*q+tau_n+mu_n)*Dc - Dc*0.059
      dDt_starc <-  sigma*g*q*Dc - (f_star+tau_star+mut_star)*Dt_starc - Dt_starc*0.059
      dDtc      <-  sigma*(1-g)*q*Dc - (f+tau+mut)*Dtc - Dtc*0.059
      dRnc      <-  tau_n*Dc - alpha*force_func3(time)*Rnc*(Dc+Da)/N - (rho_n+mu)*Rnc - Rnc*0.059
      dR_starc  <-  tau_star*Dt_starc - alpha*force_func3(time)*R_starc*(Dc+Da)/N - (rho_star+mu)*R_starc - R_starc*0.059
      dRc       <-  tau*Dtc - alpha*force_func3(time)*Rc*(Dc+Da)/N - (rho+mu)*Rc - Rc*0.059
      
      #adults
      dSa       <-  Sc*0.059 - (force_func4(time)*Sa*(Dc+Da))/N - mu*Sa
      dLea      <-  Lec*0.05 + force_func4(time)*(Sa+alpha*(Lla+Rna+R_stara+Ra))*(Dc+Da)/N - (omega+p+mu)*Lea
      dLla      <-  Llc*0.059 + omega*Lea - ((alpha*force_func4(time)*(Dc+Da))/N+v+mu)*Lla
      dDa       <-  Dc*0.059 + p*Lea + v*Lla + f*Dt_stara + f*Dta + rho_n*Rna + rho*Ra + rho_star*R_stara - (sigma*q+tau_n+mu_n)*Da
      dDt_stara <-  Dt_starc*0.059 + sigma*g*q*Da - (f_star+tau_star+mut_star)*Dt_stara
      dDta      <-  Dtc*0.059 + sigma*(1-g)*q*Da - (f+tau+mut)*Dta
      dRna      <-  Rnc*0.059 + tau_n*Da - alpha*force_func4(time)*Rna*(Dc+Da)/N - (rho_n+mu)*Rna
      dR_stara  <-  R_starc*0.059 + tau_star*Dt_stara - alpha*force_func4(time)*R_stara*(Dc+Da)/N - (rho_star+mu)*R_stara
      dRa       <-  Rc*0.059 + tau*Dta - alpha*force_func4(time)*Ra*(Dc+Da)/N - (rho+mu)*Ra
      
      
      
      dN        <- dSc+dLec+dLlc+dDc+dDtc+dDt_starc+dRc+dRnc+dR_starc+dSa+dLea+dLla+dDa+dDta+dDt_stara+dRa+dRna+dR_stara
      dInc      <- (0.7*0.47+0.3)*p*Lec + (0.7*0.47+0.3)*v*Llc + rho_n*Rnc + rho*Rc + p*Lea + v*Lla + rho_n*Rna + rho*Ra-Inc
      dDtb      <- mu_n*Da + mut_star*Dt_stara + mut*Dta + mu_n*Dc + mut_star*Dt_starc + mut*Dtc - Dtb   # + (Dc+Dtc+Da+Dta)*mu
      
      return(list(c(dSc,dLec,dLlc,dDc,dDtc,dDt_starc,dRc,dRnc,dR_starc,dSa,dLea,dLla,dDa,dDta,dDt_stara,dRa,dRna,dR_stara,dN,dInc,dDtb)))
      
    })
  }
  #yini = c(Sc=38911307, Lec=300700, Llc=2269300, Dc=6761, Dtc=8710, Dt_starc=0, Rc=29420, Rnc=1522, R_starc=0,
  #Sa=26755282, Lea=426053, Lla=22470700, Da=52596, Dta=75169, Dt_stara=0, Ra=1713042, R_stara=0, Rna=651920, N=91818000, Dead=42000, Inc=159000) #*16   Dead=23300, Inc=133000  Dead=42000, Inc=159139
  
  # yini = c(Sc=38655496, Lec=82510, Llc=2147490, Dc=4768, Dtc=33936, Dt_starc=0, Rc=30671, Rnc=129, R_starc=0,  
  #          Sa=26110378, Lea=867000, Lla=23400000, Da=33378, Dta=237549, Dt_stara=0, Ra=901, Rna=213794, R_stara=0, N=91818000, Inc=85000, Dtb = 30000) #, Inc=159000, Dtb = 30000
  
  yini = c(Sc=38655496, Lec=93531, Llc=887913, Dc=8061, Dtc=6979, Dt_starc=0, Rc=75000, Rnc=4200, R_starc=0,  
           Sa=26110378, Lea=887913, Lla=23076526, Da=64613, Dta=63583, Dt_stara=0, Ra=400000, Rna=60000, R_stara=0, N=91818000, Inc=100000, Dtb = 20490) #, Inc=159000, Dtb = 30000
  
  
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

# parameters values:
parameters_values <- c(
  pi = 0.0336,#29405.03,
  mu = 0.007,
  betac1 = 13.5, betac2 = 9,
  betaa1 = 13.5, betaa2 = 9,
  omega = 0.8720559,
  p = 0.08259749,
  v = 0.0005940438,
  tau_n = 0.179,#0.144, 
  tau = 0.9231366,
  mu_n = 0.202,#0.094,
  mut = 0.04,#0.04,   #0.019 literature #0.04 facilities data
  sigma = 0.7656967,
  g = 0,#1,              #intervention switch
  q = 0.9489609,
  f = 0.012,
  alpha = 0.2199592,
  rho_n = 0.01995215,
  rho = 0.01100774,
  #MERM
  tau_star = 0.965,
  mut_star = 0.019,
  f_star = 0.018,
  rho_star = 0.01
  #99DOTS
  # tau_star = 0.963,
  # mut_star = 0.034,
  # f_star = 0.003,
  # rho_star = 0.01
)

solution <- ode_results(parameters_values)
par(mar = c(2, 2, 2, 2))
plot(solution)
