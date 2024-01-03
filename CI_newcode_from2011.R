# library(hmer)
# library(deSolve)
# library(ggplot2)
# library(reshape2)
# library(purrr)
# library(tidyverse)
# library(lhs)
# set.seed(123)
# library("writexl")
# library(openxlsx)

excel_path <- "C:\\Users\\LaraGosce\\Documents\\ASCENT\\Ethiopia\\Ethiopia code\\results\\calibration results\\final_results7.xlsx"  #insert pathway to the calibration output
pi1wave <- read.xlsx(excel_path, sheet = 1, cols = 1)
#pi2wave <- read.xlsx(excel_path, sheet = 1, cols = 2)
mu1wave <- read.xlsx(excel_path, sheet = 1, cols = 3)
#mu2wave <- read.xlsx(excel_path, sheet = 1, cols = 4)
betacwave <- read.xlsx(excel_path, sheet = 1, cols = 5)
betaawave <- read.xlsx(excel_path, sheet = 1, cols = 6)
omegawave <- read.xlsx(excel_path, sheet = 1, cols = 7)
pwave <- read.xlsx(excel_path, sheet = 1, cols = 8)
vwave <- read.xlsx(excel_path, sheet = 1, cols = 9)
tau_nwave <- read.xlsx(excel_path, sheet = 1, cols = 10)
tauwave <- read.xlsx(excel_path, sheet = 1, cols = 11)
mu_nwave <- read.xlsx(excel_path, sheet = 1, cols = 12)
mutwave <- read.xlsx(excel_path, sheet = 1, cols = 13)
sigmawave <- read.xlsx(excel_path, sheet = 1, cols = 14)
sigmalfuwave <- read.xlsx(excel_path, sheet = 1, cols = 15)
qwave <- read.xlsx(excel_path, sheet = 1, cols = 16)
fwave <- read.xlsx(excel_path, sheet = 1, cols = 17)
alphawave <- read.xlsx(excel_path, sheet = 1, cols = 18)
rho_nwave <- read.xlsx(excel_path, sheet = 1, cols = 19)
rhowave <- read.xlsx(excel_path, sheet = 1, cols = 20)
Wave <- data.frame(pi1wave,mu1wave,betacwave,betaawave,omegawave,pwave,vwave,tau_nwave,tauwave,mu_nwave,mutwave,sigmawave,sigmalfuwave,qwave,fwave,alphawave,rho_nwave,rhowave)
n<-nrow(Wave)

#define empty lists
my_list <- list()
incidence_list <- list()
deaths_list <- list()

#define empty dataframes
#solution <- data.frame(NA_row = rep(NA, n))
incidence <- data.frame(NA_row = rep(NA, 25))
deaths <- data.frame(NA_row = rep(NA, 25))
total_pop <- data.frame(NA_row = rep(NA, 25))

#run the model iteratively for each set of parameters resulting from the calibration
for (i in 1:n){
  # eventdat <- data.frame(var=c("Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc","Inc",
  #                              "Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb","Dtb"),
  #                        time=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,
  #                               0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24),#c(1.1,2.1,3.1,4.1,5.1,6.1,7.1,8.1,9.1,10.1),
  #                        value=c(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,
  #                                0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0),
  #                        method=c("rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep",
  #                                 "rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep","rep"))
  # 
  eventdat <- data.frame(var=c(rep("Inc", 25), rep("Dtb", 25)),
                         time=c(rep(0:24, 2)),
                         value=c(rep(0, 50)),
                         method=c(rep("rep", 50)))
  
  # Define the helper function `ode_results`, to obtain the solution of the ODEs.
  ode_results <- function(parms, end_time = 24) {  #25 years from 2011 to 2035
    #interpolate births and deaths to capture change in demographics between 2020 and 2035
    forcer1 = matrix(c(0, parms['pi1'], 24, parms['pi3']), 
                     ncol = 2, byrow = TRUE)
    force_func1 = approxfun(x = forcer1[,1], y = forcer1[,2], method = "linear", rule = 2)
    forcer2 = matrix(c(0, parms['mu1'], 24, parms['mu3']), 
                     ncol = 2, byrow = TRUE)
    force_func2 = approxfun(x = forcer2[,1], y = forcer2[,2], method = "linear", rule = 2)
    forcer3 = matrix(c(0, parms['g1'], 11, parms['g2'], 12, parms['g3'], 24, parms['g3']), 
                     ncol = 2, byrow = TRUE)
    force_func3 = approxfun(x = forcer3[,1], y = forcer3[,2], method = "linear", rule = 2)
    des = function(time, state, parms) {
      with(as.list(c(state, parms)), {
        #children
        dSc       <-  force_func1(time)*N - (betac*Sc*(Dc+Da+LFUc+LFUa))/N - force_func2(time)*Sc - Sc*0.059
        dLec      <-  betac*(Sc+alpha*Llc)*(Dc+Da+LFUc+LFUa)/N - (omega+(0.7*0.29+0.3)*p+force_func2(time))*Lec - Lec*0.059
        dLlc      <-  omega*Lec - ((alpha*betac*(Dc+Da+LFUc+LFUa))/N+(0.7*0.29+0.3)*v+force_func2(time))*Llc - Llc*0.059 + 0.5*(Rnc+Rc)
        dDc       <-  (0.7*0.29+0.3)*p*Lec + (0.7*0.29+0.3)*v*Llc + rho_n*Rnc + rho*Rc - (sigma*q+tau_n+mu_n)*Dc - Dc*0.059
        dDtc      <-  sigma*q*Dc + sigma_lfu*q*LFUc - (f+tau+mut)*Dtc - Dtc*0.059
        dRnc      <-  tau_n*Dc - (rho_n+force_func2(time))*Rnc - Rnc*0.059
        dRc       <-  tau*Dtc - (rho+force_func2(time)+0.5)*Rc - Rc*0.059
        dLFUc     <-  f*Dtc - (sigma_lfu*q+mu_n)*LFUc - LFUc*0.059
        
        #adults
        dSa       <-  Sc*0.059 - (betaa*Sa*(Dc+Da+LFUc+LFUa))/N - force_func2(time)*Sa
        dLea      <-  Lec*0.059 + betaa*(Sa+alpha*Lla)*(Dc+Da+LFUc+LFUa)/N - (omega+p+force_func2(time))*Lea
        dLla      <-  Llc*0.059 + omega*Lea + 0.5*(Rna+R_stara+Ra) - ((alpha*betaa*(Dc+Da+LFUc+LFUa))/N+v+force_func2(time))*Lla 
        dDa       <-  Dc*0.059 + p*Lea + v*Lla + rho_n*Rna + rho*Ra + rho_star*R_stara - (sigma*q+tau_n+mu_n)*Da
        dDta      <-  Dtc*0.059 + (1-force_func3(time))*sigma*q*Da + (1-force_func3(time))*sigma_lfu*q*LFUa - (f+tau+mut)*Dta
        dDt_stara <-  force_func3(time)*sigma*q*Da + force_func3(time)*sigma_lfu*q*LFUa - (f_star+tau_star+mut_star)*Dt_stara
        dRna      <-  Rnc*0.059 + tau_n*Da - (rho_n+force_func2(time))*Rna
        dR_stara  <-  tau_star*Dt_stara - (rho_star+force_func2(time)+0.5)*R_stara #- alpha*force_func4(time)*R_stara*(Dc+Da)/N 
        dRa       <-  Rc*0.059 + tau*Dta - (rho+force_func2(time)+0.5)*Ra
        dLFUa     <-  LFUc*0.059 + f*Dta +f_star*Dt_stara - (sigma_lfu*q+mu_n)*LFUa
        
        dN        <- dSc+dLec+dLlc+dDc+dDtc+dRnc+dRc+dLFUc+dSa+dLea+dLla+dDa+dDta+dDt_stara+dRna+dR_stara+dRa+dLFUa
        dInc      <- (((0.7*0.29+0.3)*p*Lec + (0.7*0.29+0.3)*v*Llc + rho_n*Rnc + rho*Rc + p*Lea + v*Lla + rho_n*Rna + rho*Ra + rho_star*R_stara)/N)*100000
        dDtb      <- ((mu_n*Da + mut*Dta + mut_star*Dt_stara + mu_n*Dc + mut*Dtc + mu_n*LFUc + mu_n*LFUa)/N)*100000
        
        return(list(c(dSc,dLec,dLlc,dDc,dDtc,dRnc,dRc,dLFUc,dSa,dLea,dLla,dDa,dDta,dDt_stara,dRna,dR_stara,dRa,dLFUa,dN,dInc,dDtb)))
        
      })
    }
    #initial values 2011
    yini = c(Sc=40242416, Lec=53531, Llc=687913, Dc=10061, Dtc=6979, Rnc=1200, Rc=10000, LFUc=1000,
             Sa=21947265, Lea=1827913, Lla=28076526, Da=160000, Dta=63583, Dt_stara=0, Rna=1000, R_stara=0, Ra=100000, LFUa=8000, N=91818000, Inc=(160000/91818000)*100000, Dtb = (40000/91818000)*100000) #, Inc=159000, Dtb = 30000
    #Sa=21947265, Lea=1827913
    #Sc=40242416, Lec=53531
    
    #initial values 2020
    #yini = c(Sc=41342044, Lec=604691, Llc=4698855, Dc=29993, Dtc=16867, Rnc=34138, Rc=27853, LFUc=769,
    #         Sa=37149336, Lea=557008, Lla=33263136, Da=80740, Dta=48352, Dt_stara=0, Rna=152859, R_stara=0, Ra=95640, LFUa=2935, N=116716829, Inc=0.8607523, Dtb = 0.2024454) #, Inc=159000, Dtb = 30000
    
    #initial values 2023 
    # yini = c(Sc=43170598, Lec=551450, Llc=5260727, Dc=27224, Dtc=15040, Rnc=39178, Rc=25176, LFUc=701,
    #          Sa=42591267, Lea=547686, Lla=34974799, Da=79983, Dta=46695, Dt_stara=0, Rna=186615, R_stara=0, Ra=87585, LFUa=2585, N=126218924, Inc=0.7796431, Dtb = 0.1791367) 
    
    
    times = seq(0, end_time, by = 1)
    out = deSolve::ode(yini, times, des, parms, events = list(data = eventdat))
    return(out)
  }
  
  
  # parameters values:[27]
  parameters_values <- c(
    pi1 = pi1wave[i,1],pi3 = 0.0282,#natality
    mu1 = mu1wave[i,1],mu3 = 0.0054,#mortality
    betac = betacwave[i,1],#transmission in children
    betaa = betaawave[i,1],#transmission in adults
    omega = omegawave[i,1],#progression from early to late latency
    p = pwave[i,1],#probability of active TB from early latency
    v = vwave[i,1],#progression from late latent to active TB
    tau_n = tau_nwave[i,1],#0.144,   #recovery in untreated people
    tau = tauwave[i,1],  #0.944 (0.850-1 CI) standard care   #treatment completion (control)
    tau_star = 0.940, #0.943 (0.835-1 CI) pillbox      0.940 (0.832-1 CI) label      #treatment completion (intervention)
    mu_n = mu_nwave[i,1],#0.094,     #mortality in untreated people
    mut = mutwave[i,1],  #0.031 (0.023-0.041 CI) stand care  #mortality in people under treatment (control)
    mut_star = 0.039, #0.031 (0.022-0.043 CI) pillbox      0.039 (0.031-0.051) label      #mortality in people under treatment (intervention)
    sigma = sigmawave[i,1],   #detection 
    sigma_lfu = sigmalfuwave[i,1], #new detection of lost to follow up
    g1 = 0, #proportion of people on DAT intervention from 2011 to 2022, always 0
    g2 = 0, #proportion of people on DAT intervention from 2022 to 2023  (0 if everyone remains in SoC)
    g3 = 0, #1 #proportion of people on DAT intervention from 2023 to 2035  (0 if everyone remains in SoC)
    g3 = 0, #1 #proportion of people on DAT intervention from 2035 onward  (0 if everyone remains in SoC)
    q = qwave[i,1],     #treatment initiation 
    f = fwave[i,1],    #0.016 (0.09-0.029 CI) stand care  #LTFU (control)
    f_star = 0.006, #0.01 (0.006-0.018 CI) pillbox       0.006 (0.003-0.013 CI) label      #LTFU (intervention)
    alpha = alphawave[i,1], #protection from reinfection of of recovered individuals
    rho_n = rho_nwave[i,1], #relapse following self-cure
    rho = rhowave[i,1],    #0.16 (0.01-0.027 CI) standard care  #relapse (control)
    rho_star = 0.019       #0.019 (0.011-0.032 CI) pillbox   0.019 (0.012-0.029 CI) label   #relapse (intervention)
  )
  
  
  
  #SAVE ONLY INCIDENCE AND DEATHS AS DATAFRAMES
  new_incidence <- ode_results(parameters_values)[,21]
  incidence[,i] <- new_incidence
  new_deaths <- ode_results(parameters_values)[,22]
  deaths[,i] <- new_deaths
  totalpop <- ode_results(parameters_values)[,20]
  total_pop[,i] <- totalpop
  
}

#TRANSPOSE
incid<-t(incidence)
dead<-t(deaths)
totpop<-t(total_pop)

#CALCULATE MEAN AND CONFIDENCE INTERVALS OF INCIDENCE AND DEATHS
incide <- data.frame(NA_row = rep(NA, 5))
dea <- data.frame(NA_row = rep(NA, 5))
tpop <- data.frame(NA_row = rep(NA, 5))

for (i in 1:25){ 
  incide[i,1]<-mean(incid[,i])   #mean
  standard_error <- sd(incid[,i])/sqrt(n)    #standard error
  t_score = qt(p=0.05/2, df=n-1,lower.tail=F)  
  margin_error <- t_score * standard_error
  incide[i,2] <- incide[i,1] - margin_error  #lower boundary
  incide[i,3] <- incide[i,1] + margin_error  #upper boundary
  #incide[i,4] <- min(incid[,i])              #minimum value
  #incide[i,5] <- max(incid[,i])              #maximum value
  incide[i,4]<-quantile(incid[,i], probs=c(.025))
  incide[i,5]<-quantile(incid[,i], probs=c(.975))
  
  dea[i,1]<-mean(dead[,i])   #mean
  standard_error_dead <- sd(dead[,i])/sqrt(n)    #standard error
  t_score = qt(p=0.05/2, df=n-1,lower.tail=F)
  margin_error_dead <- t_score * standard_error_dead
  dea[i,2] <- dea[i,1] - margin_error_dead  #lower boundary
  dea[i,3] <- dea[i,1] + margin_error_dead
  #dea[i,4] <- min(dead[,i])              #minimum value
  #dea[i,5] <- max(dead[,i])              #maximum value
  dea[i,4]<-quantile(dead[,i], probs=c(.025))
  dea[i,5]<-quantile(dead[,i], probs=c(.975))
  
  tpop[i,1]<-mean(totpop[,i])   #mean
}

#calib_years <- c(2,3,4,5,6,7,8,9)
targets_incide_low <- c(75.1,76.6,78.1,95.2,87.8,87.0,85.5,83.3)
targets_incide_high <- c(296.8,278.2,254.4,208.3,201.6,186.0,172.6,159.9)
targets_dea_low  <- c(20,22,21,16,16,15,13,12)
targets_dea_high <- c(43,46,46,32,32,32,28,26)

#PLOT MEANs AND CIs
#incidence_list <- list(incide)
par(mar = c(2, 2, 2, 2))
plot(seq(2011,2035), incide[,1],type = "l",
     col = c('black'),
     xlab = "Year",
     ylab = "Incidence", lty = 1,
     xlim = c(2013, 2035),
     ylim = c(47.8, 297.0))
title(main = "TB Incidence", cex.main = 0.8,   font.main= 2, )
# lines(incide[,2],type = "l",
#       col = c('red'),     xlab = "Year",
#       ylab = "Incidence",  lty = 1)
# lines(incide[,3],type = "l",
#       col = c('red'),     xlab = "Year",
#       ylab = "Incidence",  lty = 1)
lines(seq(2011,2035), incide[,4],type = "l",
      col = c('blue'),     xlab = "Year",
      ylab = "Incidence",  lty = 1)
lines(seq(2011,2035), incide[,5],type = "l",
      col = c('blue'),     xlab = "Year",
      ylab = "Incidence",  lty = 1)
arrows(seq(2012,2019), targets_incide_low, seq(2012,2019), targets_incide_high, angle=90, code=3, length=0.06, col="blue")



par(mar = c(2, 2, 2, 2))
plot(seq(2011,2035), dea[,1],type = "l",
     col = c('black'), lty = 1,
     xlab = "Year",
     ylab = "Deaths",
     xlim = c(2013, 2035),
     ylim = c(8, 46))
title(main = "TB Mortality", cex.main = 0.8,   font.main= 2, )
# lines(seq(2011,2035), dea[,2],type = "l",
#       col = c('red'), lty = 1)
# lines(seq(2011,2035), dea[,3],type = "l",
#       col = c('red'), lty = 1)
lines(seq(2011,2035), dea[,4],type = "l",
      col = c('blue'), lty = 1)
lines(seq(2011,2035), dea[,5],type = "l",
      col = c('blue'), lty = 1)
#arrows(x, y+errors, x, y-errors, angle=90, code=3, length=0.06, col="blue")
arrows(seq(2012,2019), targets_dea_low, seq(2012,2019), targets_dea_high, angle=90, code=3, length=0.06, col="blue")


#\ASCENT\Ethiopia\Ethiopia%20code\results

incid_SoC <- incid
dead_SoC <- dead
incide_SoC <- incide
dea_SoC <- dea

#########################################################################################################################################
############################ ESTIMATE DIFFERENCE BETWEEN SoC RESULTS AND INTERVENTION RESULTS ##########################################
########################################################################################################################################

#this part of the code should be run after the intervention runs have been saved, file name CI_newcode_from2011_intervention.R

incidSoC <- data.frame(NA_row = rep(NA, nrow(incid_label)), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)#, NA_col = rep(NA, 25))
for (i in 1:nrow(incid_SoC)){ 
  for (j in 1:10){
    incidSoC[j+10*(i-1),] <- incid_SoC[i,]
  }
}

incid_pillboxVSSoC <- incid_pillbox - incidSoC
incid_labelsVSSoC <- incid_label - incidSoC

deadSoC <- data.frame(NA_row = rep(NA, nrow(dead_label)), NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA)#, NA_col = rep(NA, 25))
for (i in 1:nrow(dead_SoC)){ 
  for (j in 1:10){
    deadSoC[j+10*(i-1),] <- dead_SoC[i,]
  }
}

dead_pillboxVSSoC <- dead_pillbox - deadSoC
dead_labelsVSSoC <- dead_label - deadSoC

#CALCULATE MEAN AND CONFIDENCE INTERVALS OF INCIDENCE DIFFERENCE AND DEATHS DIFFERENCE
incide_pillboxVSSoC <- data.frame(NA_row = rep(NA, 5))
incide_labelsVSSoC <- data.frame(NA_row = rep(NA, 5))
dea_pillboxVSSoC <- data.frame(NA_row = rep(NA, 5))
dea_labelsVSSoC <- data.frame(NA_row = rep(NA, 5))

for (i in 1:25){ 
  incide_pillboxVSSoC[i,1]<-mean(incid_pillboxVSSoC[,i])   #mean
  standard_error <- sd(incid_pillboxVSSoC[,i])/sqrt(n)    #standard error
  t_score = qt(p=0.05/2, df=n-1,lower.tail=F)  
  margin_error <- t_score * standard_error
  incide_pillboxVSSoC[i,2] <- incide_pillboxVSSoC[i,1] - margin_error  #lower boundary
  incide_pillboxVSSoC[i,3] <- incide_pillboxVSSoC[i,1] + margin_error  #upper boundary
  incide_pillboxVSSoC[i,4]<-quantile(incid_pillboxVSSoC[,i], probs=c(.025))
  incide_pillboxVSSoC[i,5]<-quantile(incid_pillboxVSSoC[,i], probs=c(.975))
  
  incide_labelsVSSoC[i,1]<-mean(incid_labelsVSSoC[,i])   #mean
  standard_error <- sd(incid_labelsVSSoC[,i])/sqrt(n)    #standard error
  t_score = qt(p=0.05/2, df=n-1,lower.tail=F)  
  margin_error <- t_score * standard_error
  incide_labelsVSSoC[i,2] <- incide_labelsVSSoC[i,1] - margin_error  #lower boundary
  incide_labelsVSSoC[i,3] <- incide_labelsVSSoC[i,1] + margin_error  #upper boundary
  incide_labelsVSSoC[i,4]<-quantile(incid_labelsVSSoC[,i], probs=c(.025))
  incide_labelsVSSoC[i,5]<-quantile(incid_labelsVSSoC[,i], probs=c(.975))
  
  dea_pillboxVSSoC[i,1]<-mean(dead_pillboxVSSoC[,i])   #mean
  standard_error_dead <- sd(dead_pillboxVSSoC[,i])/sqrt(n)    #standard error
  t_score = qt(p=0.05/2, df=n-1,lower.tail=F)
  margin_error_dead <- t_score * standard_error_dead
  dea_pillboxVSSoC[i,2] <- dea_pillboxVSSoC[i,1] - margin_error_dead  #lower boundary
  dea_pillboxVSSoC[i,3] <- dea_pillboxVSSoC[i,1] + margin_error_dead
  dea_pillboxVSSoC[i,4]<-quantile(dead_pillboxVSSoC[,i], probs=c(.025))
  dea_pillboxVSSoC[i,5]<-quantile(dead_pillboxVSSoC[,i], probs=c(.975))
  
  dea_labelsVSSoC[i,1]<-mean(dead_labelsVSSoC[,i])   #mean
  standard_error_dead <- sd(dead_labelsVSSoC[,i])/sqrt(n)    #standard error
  t_score = qt(p=0.05/2, df=n-1,lower.tail=F)
  margin_error_dead <- t_score * standard_error_dead
  dea_labelsVSSoC[i,2] <- dea_labelsVSSoC[i,1] - margin_error_dead  #lower boundary
  dea_labelsVSSoC[i,3] <- dea_labelsVSSoC[i,1] + margin_error_dead
  dea_labelsVSSoC[i,4]<-quantile(dead_labelsVSSoC[,i], probs=c(.025))
  dea_labelsVSSoC[i,5]<-quantile(dead_labelsVSSoC[,i], probs=c(.975))
  
}

#PLOT MEANs AND CIs
par(mar = c(2, 2, 2, 2))
plot(seq(2011,2035), incide_pillboxVSSoC[,1],type = "l",
     col = c('black'),
     xlab = "Year",
     ylab = "Incidence", lty = 1,
     xlim = c(2013, 2035),
     ylim = c(-6, 2.7))
title(main = "TB Difference in Incidence Pillbox vs SoC", cex.main = 0.8,   font.main= 2, )
lines(seq(2011,2035), incide_pillboxVSSoC[,4],type = "l",
      col = c('blue'),     xlab = "Year",
      ylab = "Incidence",  lty = 1)
lines(seq(2011,2035), incide_pillboxVSSoC[,5],type = "l",
      col = c('blue'),     xlab = "Year",
      ylab = "Incidence",  lty = 1)


par(mar = c(2, 2, 2, 2))
plot(seq(2011,2035), incide_labelsVSSoC[,1],type = "l",
     col = c('black'),
     xlab = "Year",
     ylab = "Difference in Incidence", lty = 1,
     xlim = c(2013, 2035),
     ylim = c(-8, 2.7))
title(main = "TB Incidence Label vs SoC", cex.main = 0.8,   font.main= 2, )
lines(seq(2011,2035), incide_labelsVSSoC[,4],type = "l",
      col = c('blue'),     xlab = "Year",
      ylab = "Incidence",  lty = 1)
lines(seq(2011,2035), incide_labelsVSSoC[,5],type = "l",
      col = c('blue'),     xlab = "Year",
      ylab = "Incidence",  lty = 1)



par(mar = c(2, 2, 2, 2))
plot(seq(2011,2035), dea_pillboxVSSoC[,1],type = "l",
     col = c('black'), lty = 1,
     xlab = "Year",
     ylab = "Difference in Deaths",
    xlim = c(2013, 2035),
    ylim = c(-1.5, 1))
title(main = "TB Mortality Pillbox vs SoC", cex.main = 0.8,   font.main= 2, )
lines(seq(2011,2035), dea_pillboxVSSoC[,4],type = "l",
      col = c('blue'), lty = 1)
lines(seq(2011,2035), dea_pillboxVSSoC[,5],type = "l",
      col = c('blue'), lty = 1)

#Estimate the probability of the intervention being better than the SoC i.e. how many runs return lower incidence and lower deaths
par(mar = c(2, 2, 2, 2))
plot(seq(2011,2035), dea_labelsVSSoC[,1],type = "l",
     col = c('black'), lty = 1,
     xlab = "Year",
     ylab = "Difference in Deaths",
     xlim = c(2013, 2035),
     ylim = c(-1.5, 1))
title(main = "TB Mortality Label vs SoC", cex.main = 0.8,   font.main= 2, )
lines(seq(2011,2035), dea_labelsVSSoC[,4],type = "l",
      col = c('blue'), lty = 1)
lines(seq(2011,2035), dea_labelsVSSoC[,5],type = "l",
      col = c('blue'), lty = 1)

for (i in 1:nrow(incid_pillboxVSSoC)){
  incid_pillboxVSSoC[i,26]<-sum(incid_pillboxVSSoC[i,13:25])
  if (incid_pillboxVSSoC[i,26]<=0){incid_pillboxVSSoC[i,27]=1}
  else
  {incid_pillboxVSSoC[i,27]=0}
  incid_labelsVSSoC[i,26]<-sum(incid_labelsVSSoC[i,13:25])
  if (incid_labelsVSSoC[i,26]<=0){incid_labelsVSSoC[i,27]=1}
  else
  {incid_labelsVSSoC[i,27]=0}
  dead_pillboxVSSoC[i,26]<-sum(dead_pillboxVSSoC[i,13:25])
  if (dead_pillboxVSSoC[i,26]<=0){dead_pillboxVSSoC[i,27]=1}
  else
  {dead_pillboxVSSoC[i,27]=0}
  dead_labelsVSSoC[i,26]<-sum(dead_labelsVSSoC[i,13:25])
  if (dead_labelsVSSoC[i,26]<=0){dead_labelsVSSoC[i,27]=1}
  else
  {dead_labelsVSSoC[i,27]=0}
}

p1<-sum(incid_pillboxVSSoC[,27])/nrow(incid_pillboxVSSoC)
p2<-sum(incid_labelsVSSoC[,27])/nrow(incid_labelsVSSoC)
p3<-sum(dead_pillboxVSSoC[,27])/nrow(dead_pillboxVSSoC)
p4<-sum(dead_labelsVSSoC[,27])/nrow(dead_labelsVSSoC)



