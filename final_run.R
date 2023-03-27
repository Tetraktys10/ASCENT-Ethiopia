eventdat <- data.frame(var=c("Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb"),
                       time=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),#c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1),
                       value=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
                       method=c("rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep"))

# Define the helper function `ode_results`, to obtain the solution of the ODEs.
ode_results <- function(parms, end_time = 24) {  #11 years from 2011 to 2021
  forcer1 = matrix(c(0, parms['pi1'], 10, parms['pi2'], 24, parms['pi3']), 
                   ncol = 2, byrow = TRUE)
  force_func1 = approxfun(x = forcer1[,1], y = forcer1[,2], method = "linear", rule = 2)
  forcer2 = matrix(c(0, parms['mu1'], 10, parms['mu2'], 24, parms['mu3']), 
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
  # yini = c(Sc=40242416, Lec=93531, Llc=887913, Dc=8061, Dtc=6979, Rc=75000, Rnc=4200,  
  #          Sa=25947265, Lea=887913, Lla=23076526, Da=64613, Dta=63583, Ra=400000, Rna=60000, N=91818000, Inc=90000, Dtb = 20490) #, Inc=159000, Dtb = 30000
  yini = c(Sc=40242416, Lec=53531, Llc=687913, Dc=11061, Dtc=6979, Rc=10000, Rnc=1200,
           Sa=21947265, Lea=1827913, Lla=28076526, Da=168000, Dta=63583, Ra=100000, Rna=1000, N=91818000, Inc=(160000/91818000)*1000, Dtb = (40000/91818000)*1000) #, Inc=159000, Dtb = 30000
  
  
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
  pi1 = 0.035798,pi2 = 0.0324,pi3 = 0.0282,#pi3 = 0.0252,#pi = 0.029,#birth
  mu1 = 0.008496738,mu2 = 0.006798442,mu3 = 0.0054,#mu3 = 0.0054,#mu = 0.007,#mortality
  betac = 13.71359921,
  betaa = 11.19208761,
  omega = 0.87200138,
  p = 0.082605,
  v = 0.000594,
  tau_n = 0.163345,#0.144,
  tau = 0.923375,
  mu_n = 0.18517,#0.094,
  mut = 0.03873,#0.04,   #0.019 literature #0.04 facilities data
  sigma = 0.569559,
  q = 0.942588,
  f = 0.012055,
  alpha = 0.239097,
  rho_n = 0.020554,
  rho = 0.011527
)

# # parameters values post calibration:
# parameters_values <- c(
#   pi1 = 0.0358,pi2 = 0.0324,pi3 = 0.0262,#pi = 0.029,#birth 
#   mu1 = 0.0085,mu2 = 0.0068,mu3 = 0.0054,#mu = 0.007,#mortality 
#   betac = 14,  #transmission in children
#   betaa = 14,  #transmission in adults
#   #betac1 = 13.5, betac2 = 13.0, betac3 = 12.5, betac4 = 12.0, #transmission in children
#   #betaa1 = 14.0, betaa2 = 13.5, betaa3 = 13, betaa4 = 12.5, #transmission in adults
#   omega = 0.872, #progression from early latent to late latent
#   p = 0.0826, #probability of active TB from early latent
#   v = 0.000594, #probability of active TB from late latent
#   tau_n = 0.178, #self-recovery
#   tau = 0.9233, #recovery through treatment (standard care)
#   mu_n = 0.2024, #TB mortality no treatment
#   mut = 0.03945, #TB mortality on treatment   #0.019 literature #0.04 facilities data
#   sigma = 0.65,  #detection rate
#   g = 0, #intervention coverage, 0 = no intervention
#   q = 0.9484, #treatment initiation
#   f = 0.012, #LTFU
#   alpha = 0.2274, #protection from reinfection in recovered individuals
#   rho_n = 0.02, #relapse following self-cure
#   rho = 0.0122, #relapse following treatment
#   #MERM
#   tau_star = 0.965,
#   mut_star = 0.019,
#   f_star = 0.018,
#   rho_star = 0.01
#   #99DOTS
#   # tau_star = 0.963,
#   # mut_star = 0.034,
#   # f_star = 0.003,
#   # rho_star = 0.01
# )


solution <- ode_results(parameters_values)
par(mar = c(2, 2, 2, 2))
plot(solution)

#library("writexl")
#solutionasdata<- data.frame(solution)
#write_xlsx(solutionasdata,"N:\\Documents\\solution_events.xlsx")