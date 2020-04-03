rm(list = ls())

url <- "https://github.com/ryanmaomaomao/wuhan_covid19_project/blob/master/china_cases.rda?raw=true"

download.file(url, destfile= "./china_cases.RData", mode = "wb")
load("./china_cases.RData")



####objective 0#####
# estimate parameters
library(deSolve) 
sir_equations <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    dS <- -beta * I * S
    dI <-  beta * I * S - gamma * I
    dR <-  gamma * I
    return(list(c(dS, dI, dR)))
  })
}

sir_1 <- function(beta, gamma, S0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  as.data.frame(out)
}

flu2<-as.data.frame(cbind(cumsum(china_cases["Wuhan",1:22]),c(1:22)))
flu2<-as.data.frame(cbind(flu2$V2,flu2$V1))
colnames(flu2)<-c("day","cases")
flu2$cases<-flu2$cases+43
model_fit <- function(beta, gamma, data, N = 763, ...) {
  I0 <- data$cases[1] # initial number of infected (from data)
  times <- data$day   # time points (from data)
  # model's predictions:
  predictions <- sir_1(beta = beta, gamma = gamma,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0, # variables' intial values
                       times = times)                # time points
  # plotting the observed prevalences:
  with(data, plot(day, cases, ...))
  # adding the model-predicted prevalence:
  with(predictions, lines(time, I, col = "red"))
}

model_fit(beta =  0.2586/11/10^6, gamma = 0.0821, flu2, N=11*10^6, pch = 19, col = "red", ylim = c(0, 600))

### MSE ####
ss <- function(beta, gamma, data = flu, N = 763) {
  I0 <- data$cases[1]
  times <- data$day
  predictions <- sir_1(beta = beta, gamma = gamma,   # parameters
                       S0 = N - I0, I0 = I0, R0 = 0, # variables' intial values
                       times = times)                # time points
  sum((predictions$I[-1] - data$cases[-1])^2)
}

ss2_2 <- function(x) {
  ss(beta = x[1], gamma = x[2],data = flu2, N=11*10^6)
}

ss(beta = 0.2586/11/10^6, gamma = 0.0821,data = flu2, N=11*10^6)
ss(beta =  0.2/11/10^6, gamma = 0.0961,data = flu2, N=11*10^6)
model_fit(beta =  0.2586/11/10^6, gamma = 0.0821, flu2, N=11*10^6, pch = 19, col = "red", ylim = c(0, 600))
model_fit(beta =  0.2/11/10^6, gamma = 0.0961, flu2, N=11*10^6, pch = 19, col = "red", ylim = c(0, 600))

starting_param_val <- c(0.2/11/10^66, 0.0961)
ss_optim <- optim(starting_param_val, ss2_2)
ss_optim
ss_optim$par
ss(ss_optim$par[1],ss_optim$par[2],data = flu2, N=11*10^6);ss(beta = 0.2586/11/10^6, gamma = 0.0821,data = flu2, N=11*10^6)

predictions <- sir_1(beta = ss_optim$par[1], gamma = ss_optim$par[2],   # parameters
                     S0 = 11*10^6 - 43, I0 = 43, R0 = 0, # variables' intial values
                     times = flu2$day)      
model_fit(beta = ss_optim$par[1], gamma = ss_optim$par[2], N=11*10^6,flu2, pch = 19, col = "red", ylim = c(0, 600))
segments(flu2$day, flu2$cases, predictions$time, predictions$I)


####MLE####
time_points <- seq(min(flu2$day), max(flu2$day), le = 100) # vector of time points
mLL_pois <- function(beta, gamma, day, cases, N = 11*10^6) {
  beta <- exp(beta) # to make sure that the parameters are positive
  gamma <- exp(gamma)
  #  sigma <- exp(sigma)
  I0 <- cases[1] # initial number of infectious
  observations <- cases[-1] # the fit is done on the other data points
  predictions <- sir_1(beta = beta, gamma = gamma,
                       S0 = N - I0, I0 = I0, R0 = 0, times = day)
  predictions <- predictions$I[-1] # removing the first point too
  if (any(predictions < 0)) return(NA) # safety
  # returning minus log-likelihood:
  #  -sum(dnorm(x = observations, mean = predictions, sd = sigma, log = TRUE))
  -sum(dpois(x = observations, lambda = predictions, log = TRUE))
}

mLL_pois(beta = 1.707599e-08, gamma = 9.610071e-02,day = 0:21,cases =flu2$cases,N=11*10^6)
mLL_pois(beta = 2.507599e-08, gamma = 8.410071e-02,day = 0:21,cases =flu2$cases,N=11*10^6)

ss_optim$par
starting_param_val <- list(beta = ss_optim$par[1], gamma = ss_optim$par[2])
estimates_pois <- mle2(minuslogl = mLL_pois,method = "Nelder-Mead",
                       start = lapply(starting_param_val, log),
                       data = c(flu2, N = 11*10^6))
exp(coef(estimates_pois))
vcov(estimates_pois)

# points estimates:
param_hat <- exp(coef(estimates_pois))
# R_0 by our MLE
param_hat[1]/param_hat[2]*11*10^6
param_hat[1]*11*10^6

# model's best predictions:
best_predictions <- sir_1(beta = param_hat["beta"], gamma = param_hat["gamma"],
                          S0 = 11*10^6 - 43, I0 = 43, R0 = 0, time_points)$I
# confidence interval of the best predictions:
cl <- 0.99 # confidence level
cl <- (1 - cl) / 2
lwr <- qpois(p = cl, lambda = best_predictions)
upr <- qpois(p = 1 - cl, lambda = best_predictions)
# layout of the plot:
# figure 1
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure1.png",
    units="px", width=1600, height=1000, res=300)

plot(time_points, time_points, ylim = c(0, max(upr)+5), type = "n",
     xlab = "time (days)", ylab = "Prevalence")
# adding the predictions' confidence interval:
sel <- time_points >= 1 # predictions start from the second data point
polygon(c(time_points[sel], rev(time_points[sel])), c(upr[sel], rev(lwr[sel])),
        border = NA, col = adjustcolor("red", alpha.f = 0.1))
# adding the model's best predictions:
lines(time_points, best_predictions, col = "red")
# adding the observed data:
with(flu2, points(day, cases, pch = 19, col = "red"))
legend("topleft",inset = 0.02,legend=c("Best Prediction (99% CI)","Actual cases"),
       lty=c(1,NA),pch=c(NA,19), cex=0.6, col = "red")
dev.off()

-logLik(estimates_pois)

#compare the R0####
ss_optim$par[1]*11*10^6/ss_optim$par[2]
param_hat["beta"]/param_hat["gamma"]*11*10^6


#### objective 1####
# qurantine effect sir_2 model
# t 0 for taking action
# pct for betta reduction
sir_2 <- function(beta, gamma, S0, I0, R0, times,t0, pct) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values_0 <- c(beta  = beta, gamma = gamma)
  parameters_values_1 <- c(beta  = beta * (1-pct), gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out_0 <- ode(initial_values, seq(0,t0), sir_equations, parameters_values_0)
  out_1 <- ode(c(out_0[t0+1,][2:4]), seq(t0+1,length(times)), sir_equations, parameters_values_1)
  # returning the output:
  as.data.frame(rbind(out_0,out_1[-1,]))
}

# baseline model
par(mfrow = c(1,1))
# 0.167 = param_hat["beta"]*11*10^6

# Figure 2
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure2.png",
  units="px", width=1600, height=1000, res=300)
flu1<-sir_1(beta = 0.167/11/10^6, gamma = param_hat["gamma"],S0 = 11*10^6, I0 = 43, R0 = 0, times = seq(0, 400))
plot(flu1$I,type = "l",ylim = c(0,max(flu1$S)),col = "red", xlim = c(0,400),xlab = "days",ylab= "Number of cases")
with(flu1,lines(S,col="blue"))
with(flu1,lines(R,col="green"))
legend("topright",inset = 0.02, legend=c("Suseptible","Infected","Cumulative Cases"),
       col=c( "blue","red","green"), lty=c(1,1,1), cex=0.6)
dev.off()
plot((flu1$I+flu1$R),type = "l",ylim = c(0,max(flu1$I+flu1$R)),col = "red", xlim = c(0,400))

#1.1
# action at t_0 = 21, beta reduce by pct 20,40,60
# 1.23 quarantine started. so t = 23
flu2_20<-sir_2(beta = 0.167/11/10^6, gamma = param_hat["gamma"],S0 = 11*10^6, I0 = 43,
               R0 = 0, times = seq(0, 400), t0 = 23, pct = 0.2)
flu2_40<-sir_2(beta = 0.167/11/10^6, gamma = param_hat["gamma"],S0 = 11*10^6, I0 = 43,
               R0 = 0, times = seq(0, 400), t0 = 23, pct = 0.4)
flu2_60<-sir_2(beta = 0.167/11/10^6, gamma = param_hat["gamma"],S0 = 11*10^6, I0 = 43,
               R0 = 0, times = seq(0, 400), t0 = 23, pct = 0.6)

flu_eg<-sir_2(beta = 0.167/11/10^6, gamma = param_hat["gamma"],S0 = 11*10^6, I0 = 43,
               R0 = 0, times = seq(0, 400), t0 = 130, pct = 0.6)
# figure 3
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure3.png",
    units="px", width=1600, height=1000, res=300)
plot(flu1$I,type = "l",ylim = c(0,max(flu1$R[1:400])),col = "red", xlim = c(0,400),xlab = "days",ylab= "Number of cases")
with(flu1,lines(R,col="green"))

with(flu_eg,lines(R,col="green",lty=2))
with(flu_eg,lines(I,col="red",lty=2))
legend("topleft",inset = 0.02, legend=c("Inefcted without quarantine","Cumulative cases without qurantine","Infected with quarantine'","Cumulative cases with qurantine"),
       col=c("red","green"), lty=c(1,1,2,2), cex=0.6)
dev.off()

# Figure 4
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure4.png",
    units="px", width=2000, height=1400, res=300)
plot(flu1$I,type = "l",ylim = c(0,max(flu1$I)),col = "black", xlim = c(0,400), ylab = "Epedemic size", xlab = "days")
with(flu2_20,lines(I,col="blue",lty = 2))
with(flu2_40,lines(I,col="green",lty = 2))
with(flu2_60,lines(I,col="red",lty = 2))
legend("topright",inset = 0.02,title="Reduction in beta", legend=c("0", "20%","40%","60%"),
       col=c("black", "blue","green","red"), lty=c(1,2,2,2), cex=0.6)
dev.off()


# Figure 5
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure5.png",
    units="px", width=2000, height=1400, res=300)
plot(flu1$I/flu1$I,type = "l",ylim = c(0,1.1),col = "black",xlim = c(1,100), ylab = "% reduce in epedemic size", xlab = "days")
with(flu2_20,lines(I/flu1$I,col="blue",lty = 2))
with(flu2_40,lines(I/flu1$I,col="green",lty = 2))
with(flu2_60,lines(I/flu1$I,col="red",lty = 2))
legend("topright",inset = 0.02,title="Reduction in beta", legend=c("Unchecked case", "20% reduction in beta","40% reduction in beta","60% reduction in beta"),
       col=c("black", "blue","green","red"), lty=c(1,2,2,2), cex=0.6)
dev.off()

max(flu1$I);max(flu2_20$I);max(flu2_40$I);max(flu2_60$I)

flu1$time[flu1$I == max(flu1$I)];flu2_20$time[flu2_20$I == max(flu2_20$I)];flu2_40$time[flu2_40$I == max(flu2_40$I)];flu2_60$time[flu2_60$I == max(flu2_60$I)]

# 1.2
max(flu1$R);max(flu2_20$R);max(flu2_40$R);max(flu2_60$R)


#### objective 2####
# specify parameters
C_H <- 58.75
C_K <- 97
C_ICU <- 1212
L_H <-12.8
L_ICU <- 28
P_H <- 0.19
P_ICU <- 0.005
P_M <- 0.023
w <- 47808/7/365 # wage per day in $
YPLL <- 81.29-69.9
INCUBATION_PERIODS <- 14
C_FIXED <-1000000000 # for sensitivity test 10 times it
C_FIXED<-100000000 # Referring to Toronto paper.
K<-200 # contacts
K<-100 # contacts

# 2.1 Direct and indirect Cost of baseline model
# direct cost
C_D_p <- P_H * L_H * C_H + P_ICU * L_ICU * C_ICU # per person
C_D_0 <- (C_D_p)*max(flu1$R)
C_D_0
# indirect cost
C_I_p <- (P_H * L_H + P_ICU * L_ICU)*w + P_M*YPLL*w*365  # per person
C_I_0 <- C_I_p*max(flu1$R)
C_I_0

# cost of unchecked case
C_0<-C_D_0+C_I_0



# 2.2 Direct and indirect Cost of action 1 (quaratine) model beta ->20
# direct cost
C_D_20 <- (C_D_p)*max(flu2_20$R)+ C_FIXED
C_D_20
# indirect cost
C_I_20 <- C_I_p*max(flu2_20$R) + (C_K+w*INCUBATION_PERIODS)*(max(flu2_20$R)-flu2_20$R[23])*(K)*(1-0.2)
C_I_20

# 2.3 Direct and indirect Cost of action 2 (quaratine) model beat ->40
# direct cost
C_D_40 <- (C_D_p)*max(flu2_40$R)+ C_FIXED
C_D_40
# indirect cost
C_I_40 <- C_I_p*max(flu2_40$R) + (C_K+w*INCUBATION_PERIODS)*(max(flu2_40$R)-flu2_40$R[23])*(K)*(1-0.4)
C_I_40

# 2.4 Direct and indirect Cost of action 2 (quaratine) model beat ->60
# direct cost
C_D_60 <- (C_D_p)*max(flu2_60$R) + C_FIXED
C_D_60

# indirect cost
C_I_60 <- C_I_p*max(flu2_60$R) + (C_K+w*INCUBATION_PERIODS)*(max(flu2_60$R)-flu2_60$R[23])*(K)*(1-0.6)
C_I_60

C_D_0;C_D_20;C_D_40;C_D_60
C_I_0;C_I_20;C_I_40;C_I_60

C_D_0+C_I_0;C_D_20+C_I_20;C_D_40+C_I_40;C_D_60+C_I_60

# 2.4 find the optimal quarantine level
C_T_pct <- function(pct,option="T") {
  flu2_pct<-sir_2(beta = 0.167/11/10^6, gamma = param_hat["gamma"], S0 = 11*10^6, I0 = 43,
                  R0 = 0, times = seq(0, 400), t0 = 23, pct = pct)
  C_D <- (C_D_p)*max(flu2_pct$R) + C_FIXED
  q<-(max(flu2_pct$R)-flu2_pct$R[23])*(K)*(1-pct)
  i<-max(flu2_pct$R)
  C_I <- C_I_p*max(flu2_pct$R) + ifelse(q>11000000,11000000,q)*(C_K+w*INCUBATION_PERIODS)
  if (option=="T"){
    C_D+C_I
  } else {
    if (option=="I"){C_I} else {C_D}
  }
}
C_D_0+C_I_0;C_D_20+C_I_20;C_D_40+C_I_40;C_D_60+C_I_60

C_T_pct(0,"I");C_T_pct(0.2,"I");C_T_pct(0.4,"I");C_T_pct(0.6,"I")
C_T_pct(0,"D");C_T_pct(0.2,"D");C_T_pct(0.4,"D");C_T_pct(0.6,"D")
C_0;C_T_pct(0.01);C_T_pct(0.05);C_T_pct(0.2);C_T_pct(0.4);C_T_pct(0.6)

# savings
-c(C_T_pct(0),C_T_pct(0.2),C_T_pct(0.4),C_T_pct(0.6))+rep(C_0,4)


# so the optimal policy at t0=23 is to quarantine as much as possible after a certain range
# use optim
optim(0.5,C_T_pct)
# optim > 1 so sould apply the quarantine policy to 100%

#figure 6
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure6.png",
    units="px", width=2000, height=1400, res=300)
a<-c(C_T_pct(0),C_T_pct(0.1),C_T_pct(0.2),C_T_pct(0.3),C_T_pct(0.4),C_T_pct(0.5),C_T_pct(0.6),C_T_pct(0.7),C_T_pct(0.8))
b<-c(C_T_pct(0,"D"),C_T_pct(0.1,"D"),C_T_pct(0.2,"D"),C_T_pct(0.3,"D"),C_T_pct(0.4,"D"),C_T_pct(0.5,"D"),C_T_pct(0.6,"D"),C_T_pct(0.7,"D"),C_T_pct(0.8,"D"))
c<-c(C_T_pct(0,"I"),C_T_pct(0.1,"I"),C_T_pct(0.2,"I"),C_T_pct(0.3,"I"),C_T_pct(0.4,"I"),C_T_pct(0.5,"I"),C_T_pct(0.6,"I"),C_T_pct(0.7,"I"),C_T_pct(0.8,"I"))
plot(x = seq(0,0.8,0.1), y = a, ylim = c(min(c),max(a)),type = "l",xlab ="Reduction in beta",ylab = "Cost in USD",col = "green",lty = 2)
lines(x = seq(0,0.8,0.1), y = b, type = "l",col = "red",lty = 2)
lines(x = seq(0,0.8,0.1), y = c, type = "l",col = "blue",lty = 2)
lines(x = seq(0,0.8,0.1), y = rep(C_0,9),lty = 1)
legend("topright",inset = 0.02, legend=c("Quarantine-total cost","Quarantine-direct cost","Quarantine-indirect cost","Unchecked-total cost"),
       col=c("green","red","blue","black"), lty=c(2,2,2,1), cex=0.6)
# text(x = 0.7, y = 150000000000, label = "")
dev.off()

# Figure 7

png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure7.png",
    units="px", width=1600, height=1400, res=300)
e<- -rep(C_0,10)
a<- -e-c(C_T_pct(0),C_T_pct(0.05),C_T_pct(0.1),C_T_pct(0.2),C_T_pct(0.3),C_T_pct(0.4),C_T_pct(0.5),C_T_pct(0.6),C_T_pct(0.7),C_T_pct(0.8))
plot(x = seq(0,0.8,0.1), y = rep(0,9), ylim = c(min(a),max(a)), type = "l",xlab ="Reduction in beta",ylab = "Net savings in USD",col = "black",lty = 1)
lines(x = c(0,0.05,seq(0.1,0.8,0.1)), y = a, type = "l",col = "black",lty=2)
legend("topleft",inset = 0.02, legend=c("Net savings by quarantine of different levels","Unchecked case"),
       col="black", lty=c(2,1), cex=0.6,box.lwd =0,bg="transparent")
text(x = 0.7, y = 1000000000, label = "Baseline (y=0)",cex=0.6)
dev.off()

#### objective 3####
# simulation game
# now consider different quarantine action time t


C_T_pct_t0 <- function(pct,t0) {
  flu2_pct<-sir_2(beta = 0.167/11/10^6, gamma = param_hat["gamma"], S0 = 11*10^6, I0 = 43,
                  R0 = 0, times = seq(0, 400), t0 = t0, pct = pct)
  C_D <- (C_D_p)*max(flu2_pct$R) + C_FIXED
  q <- (max(flu2_pct$R)-flu2_pct$R[t0])*(K)*(1-pct)
  i<-max(flu2_pct$R)
  C_I <- C_I_p*max(flu2_pct$R) + ifelse(q>11000000,11000000,q)*(C_K+w*INCUBATION_PERIODS)
  C_D+C_I
}


# play around
C_T_pct(0);C_T_pct(0.4);C_T_pct(0.6)
C_T_pct_t0(0,23);C_T_pct_t0(0,43);C_T_pct_t0(0,63);C_T_pct_t0(0,83);C_T_pct_t0(0,103)
C_T_pct_t0(0.1,23);C_T_pct_t0(0.1,43);C_T_pct_t0(0.1,63);C_T_pct_t0(0.1,83);C_T_pct_t0(0.1,83)
C_T_pct_t0(0.2,23);C_T_pct_t0(0.2,43);C_T_pct_t0(0.2,63);C_T_pct_t0(0.2,83);C_T_pct_t0(0.2,83)
C_T_pct_t0(0.4,23);C_T_pct_t0(0.4,43);C_T_pct_t0(0.4,63);C_T_pct_t0(0.4,83);C_T_pct_t0(0.4,83)
C_T_pct_t0(0.6,23);C_T_pct_t0(0.6,43);C_T_pct_t0(0.6,63);C_T_pct_t0(0.6,83);C_T_pct_t0(0.6,83)
C_T_pct_t0(0.8,23);C_T_pct_t0(0.8,43);C_T_pct_t0(0.8,63);C_T_pct_t0(0.8,83);C_T_pct_t0(0.8,83)
# so after some perrod the qurantine is not useful in cost savings
# use optim
# envolope curve




C_T_pct_t_V <- function(pct,t) {
  p<-0
  for (i in 1:length(t)){
    t0=t[i]
    p<-c(p,C_T_pct_t0(pct,t0))
  }
  return(p[2:length(p)])
}
# uppper limit of days = 1000

# Figure 8
png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure8.png",
    units="px", width=1600, height=1400, res=300)
costs_plot <- function(period = 100) {
  x<- seq(23,23+period-1,1)
  h<- rep(C_0,period)
  plot(x = x,y=h, ylim = c(0,50000000000),xlim = c(23,23+period),type = "l",ylab = "Total costs in USD",xlab = "Action time (days since Jan 1)")
  a<- C_T_pct_t_V(0.2,x)
  points(x = x,y=a,col = "orange",type = "l",lty=c(2))
  b<- C_T_pct_t_V(0.4,x)
  points(x = x,y=b,col = "red",type = "l",lty=c(2))
  c<- C_T_pct_t_V(0.6,x)
  points(x = x,y=c,col = "green",type = "l",lty=c(2))
  d<- C_T_pct_t_V(0.8,x)
  points(x = x,y=d,col = "blue",type = "l",lty=c(2))
  legend("topright",inset = 0.01,title="Reduction in beta", legend=c("0 (Unchecked)","20%", "40%","60%","80%"),
         col=c("black","orange", "red","green","blue"), lty=c(1,2,2,2,2,2), cex=0.6,box.lwd =0,bg="transparent")
}
costs_plot(165)
abline(v=153,lty=2)
text(x = 148, y = 30000000000, label = "t=152", srt = 90)
dev.off()


# Figure 9
savings_plot <- function(period = 100) {
  x<- seq(23,23+period-1,1)
  h<- -rep(C_0,period)
  plot(x = x,y=rep(0,period), ylim = c(-10000000000,20000000000),xlim = c(23,23+period),type = "l",ylab = "Net savings in USD",xlab = "Action time (days since Jan 1)")
  a<- -h-C_T_pct_t_V(0.2,x)
  points(x = x,y=a,col = "orange",type = "l",lty=c(2))
  b<- -h-C_T_pct_t_V(0.4,x)
  points(x = x,y=b,col = "red",type = "l",lty=c(2))
  c<- -h-C_T_pct_t_V(0.6,x)
  points(x = x,y=c,col = "green",type = "l",lty=c(2))
  d<- -h-C_T_pct_t_V(0.8,x)
  points(x = x,y=d,col = "blue",type = "l",lty=c(2))
  legend("topright",inset = 0.01,title="Reduction in beta", legend=c("0 (Unchecked)","20%", "40%","60%","80%"),
         col=c("black","orange", "red","green","blue"), lty=c(1,2,2,2,2,2), cex=0.6,box.lwd =0,bg="transparent")
}

png(file="/Users/ryanchen/OneDrive/2020 Spring/comm 415 quantitive methods/project/figures/figure9.png",
    units="px", width=1600, height=1400, res=300)
savings_plot(170)
abline(v=153,lty=2)
text(x = 148, y = 10000000000, label = "t=152", srt = 90)
dev.off()

