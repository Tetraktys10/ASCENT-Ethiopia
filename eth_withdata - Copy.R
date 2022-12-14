#ASCENT model Ethiopia
incid <- source("inc_data.r")
eth_2 <- function(pi, mu, beta, omega, p, v, tau_n, tau, tau_star, mu_n, mut, mut_star, sigma, g, q, f, f_star, alpha, rho_n, rho, rho_star, 
                  S0c, Le0c, Ll0c, D0c, Dt0c, Dt_star0c, R0c, Rn0c, R_star0c,
                  S0a, Le0a, Ll0a, D0a, Dt0a, Dt_star0a, R0a, Rn0a, R_star0a, N0, Inc0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      #children
      dSc       <-  pi - (beta*Sc*(Dc+Da))/N - (mu*1000/N)*Sc - Sc*0.066
      dLec      <-  beta*(Sc+alpha*(Llc+Rnc+R_starc+Rc))*(Dc+Da)/N - (omega+(0.7*0.47+0.3)*p+mu*1000/N)*Lec - Lec*0.066
      dLlc      <-  omega*Lec - ((alpha*beta*(Dc+Da))/N+(0.7*0.47+0.3)*v+mu*1000/N)*Llc - Llc*0.066
      dDc       <-  (0.7*0.47+0.3)*p*Lec + (0.7*0.47+0.3)*v*Llc + f*Dt_starc + f*Dtc + rho_n*Rnc + rho*Rc - (sigma*q+tau_n+mu_n)*Dc - Dc*0.066
      dDt_starc <-  sigma*g*q*Dc - (f_star+tau_star+mut_star)*Dt_starc - Dt_starc*0.066
      dDtc      <-  sigma*(1-g)*q*Dc - (f+tau+mut)*Dtc - Dtc*0.066
      dRnc      <-  tau_n*Dc - alpha*beta*Rnc*(Dc+Da)/N - (rho_n+mu*1000/N)*Rnc - Rnc*0.066
      dR_starc  <-  tau_star*Dt_starc - alpha*beta*R_starc*(Dc+Da)/N - (rho_star+mu*1000/N)*R_starc - R_starc*0.066
      dRc       <-  tau*Dtc - alpha*beta*Rc*(Dc+Da)/N - (rho+mu*1000/N)*Rc - Rc*0.066
      
      #adults
      dSa       <-  Sc*0.066 - (beta*Sa*(Dc+Da))/N - (mu*1000/N)*Sa
      dLea      <-  Lec*0.066 + beta*(Sa+alpha*(Lla+Rna+R_stara+Ra))*(Dc+Da)/N - (omega+p+mu*1000/N)*Lea
      dLla      <-  Llc*0.066 + omega*Lea - ((alpha*beta*(Dc+Da))/N+v+mu*1000/N)*Lla
      dDa       <-  Dc*0.066 + p*Lea + v*Lla + f*Dt_stara + f*Dta + rho_n*Rna + rho*Ra - (sigma*q+tau_n+mu_n)*Da
      dDt_stara <-  Dt_starc*0.066 + sigma*g*q*Da - (f_star+tau_star+mut_star)*Dt_stara
      dDta      <-  Dtc*0.066 + sigma*(1-g)*q*Da - (f+tau+mut)*Dta
      dRna      <-  Rnc*0.066 + tau_n*Da - alpha*beta*Rna*(Dc+Da)/N - (rho_n+mu*1000/N)*Rna
      dR_stara  <-  R_starc + tau_star*Dt_stara - alpha*beta*R_stara*(Dc+Da)/N - (rho_star+mu*1000/N)*R_stara
      dRa       <-  Rc*0.066 + tau*Dta - alpha*beta*Ra*(Dc+Da)/N - (rho+mu*1000/N)*Ra
      
      dN        <- Sc+Lec+Llc+Dc+Dtc+Dt_starc+Rc+Rnc+R_starc+Sa+Lea+Lla+Da+Dta+Dt_stara+Ra+Rna+R_stara
      dInc      <- (0.7*0.47+0.3)*p*Lec + (0.7*0.47+0.3)*v*Llc + rho_n*Rnc + rho*Rc + p*Lea + v*Lla + rho_n*Rna + rho*Ra
      
      return(list(c(dSc,dLec,dLlc,dDc,dDtc,dDt_starc,dRc,dRnc,dR_starc,dSa,dLea,dLla,dDa,dDta,dDt_stara,dRa,dRna,dR_stara,dN,dInc)))
      #return(list(c(dSc,dLec,dLlc,dDc,dDt_starc,dDtc,dRnc,dR_starc,dRc,dSa,dLea,dLla,dDta,dDt_stara,dDa,dRna,dR_stara,dRa,dN,dInc)))
      })
  }
  
  # parameters values:
  parameters_values <- c(
    pi = 30.9*1000,
    mu = 6.5,
    beta = 0.5,
    omega = 0.872,
    p = 0.0826,
    v = 0.000594,
    tau_n = 0.027, 
    tau = 0.94,
    mu_n = 0.047,
    mut = 0.036,
    sigma = 0.71,
    g = 0,
    q = 1,
    f = 0.013,
    alpha = 0.21,
    rho_n = 0.02,
    rho = 0.01,
    tau_star = 0.96,
    mut_star = 0.02,
    f_star = 0.011,
    rho_star = 0.01
  )
  
  # the initial values of variables:
  initial_values <- c(Sc=S0c, Lec=Le0c, Llc=Ll0c, Dc=D0c, Dtc=Dt0c, Dt_starc=Dt_star0c, Rc=R0c, Rnc=Rn0c, R_starc=R_star0c,
                      Sa=S0a, Lea=Le0a, Lla=Ll0a, Da=D0a, Dta=Dt0a, Dt_stara=Dt_star0a, Ra=R0a, Rna=Rn0a, R_stara=R_star0a, N=N0, Inc=Inc0)
  
  # solving
  out <- ode(initial_values, times, equations, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}



#2011 Ethiopian population: 40955 children (0-14), 50863 adults (15+)
#2011 Ethiopian TB prevalence (WHO): 237 (191-288) x100,000  (217,609)    
#2011 Ethiopian TB prevalence (survey): Prevalence pulmonary TB( smear negative and smear positive TB) : 156/100, 000 population (if ETB is included: 240 (182-298)/100,000)) 10.5% of total new cases were in children  ; 143,237 prevalent cases in the whole population: 17,905 in children and 125,339 in adults
#2011 notification rate (survey): The NTP's notification rate of all forms of TB was 175/100,000 in 2010/11. Given that the average duration of treatment is 8 months, it is estimated that 242/100,000 adults in this survey populations received TB treatment in a year. 
#2011 number of PTB cases notified (WHO): 105,091 with a total detection rate of 72-73% (divide by 12 and multiply by 7 to consider people still on treatment since the previous year):271,485
#2011 Ethiopian LTBI prevalence (Houben et al.): 2,230,000 (2,050,000 -  2,370,000) in 0-15
#                                                867,000 (692,000 -  1,060,000) 2 years infection
#                                                23,400,000 (21,500,000 - 25,800,000) all LTBI
predictions <- eth_2(parameters_values, S0c=38655496, Le0c=82510, Ll0c=2147490, D0c=4768, Dt0c=33936, Dt_star0c=0, R0c=30671, Rn0c=129, R_star0c=0,  
      S0a=26110378, Le0a=867000, Ll0a=23400000, D0a=33378, Dt0a=237549, Dt_star0a=0, R0a=901, Rn0a=213794, R_star0a=0, N0=91818000, Inc0=233000, times = seq(0, 10))

with(predictions, {
  # plotting the time series of incident cases:
  plot(time, Inc, type = "l", col = "blue",
       xlab = "years", ylab = "incidence")
})

# adding a legend:
legend("right", c("incidence"),
       col = c("blue"), lty = 1, bty = "n")