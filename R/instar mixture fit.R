library(gulf.graphics)

x <- read.scsbio(2020:2022)

loglike <- function(theta, xi, xp, xm, fixed){
   if (!missing(fixed)) theta <- c(theta, fixed)
   
   # Parse parameter vector:
   logit.p.immature <- theta[grep("logit.p.immature", names(theta))]
   p.immature <- exp(logit.p.immature) / (1 + sum(exp(logit.p.immature)))
   p.immature[length(p.immature)+1] <- 1 - sum(p.immature)
   mu.immature <- c(theta["mu0"], rep(NA, length(p.immature)-1))
   sigma <- exp(theta["log.sigma"])
   for (i in 2:length(p.immature)){
      mu.immature[i]    <- theta["a"] * mu.immature[i-1] + theta["b"]
   }
   
   # Immature mixture density:
   ni <- as.numeric(xi)
   xi <- as.numeric(names(xi))
   vi <- rep(0, length(xi))
   for (i in 1:length(p)) vi <- vi + p.immature[i] * dnorm(xi, mu.immature[i], sigma)
   ll <- - sum(ni * log(vi)) 
   
   # Pubescent:
   if (!missing(xp)){
      logit.p.pubescent <- theta[grep("logit.p.pubescent", names(theta))]
      p.pubescent <- exp(logit.p.pubescent) / (1 + sum(exp(logit.p.pubescent)))
      p.pubescent[length(p.pubescent)+1] <- 1 - sum(p.pubescent)
      
      mu.pubescent <- theta["a.pubescent"] * mu.immature + theta["b.pubescent"]
      
      np <- as.numeric(xp)
      xp <- as.numeric(names(xp)) 
      vp <- rep(0, length(xp))
      for (i in 1:length(p.pubescent)) vp <- vp + p.pubescent[i] * dnorm(xp, mu.pubescent[i], sigma)
      ll <- ll - sum(np * log(vp)) 
   }

   # Mature:
   if (!missing(xm)){
     logit.p.mature <- theta[grep("logit.p.mature", names(theta))]
     p.mature <- exp(logit.p.mature) / (1 + sum(exp(logit.p.mature)))
     p.mature[length(p.mature)+1] <- 1 - sum(p.mature)
     
     mu.mature <- theta["a.mature"] * mu.pubescent + theta["b.mature"]

     nm <- as.numeric(xm)
     xm <- as.numeric(names(xm)) 
     vm <- rep(0, length(xm))
     
     for (i in 1:length(p.mature)) vm <- vm + p.mature[i] * dnorm(xm, mu.mature[i], sigma)
     ll <- ll - sum(nm * log(vm)) 
   }
   
   return(ll)
}

# Prepare data:
xi <- x$carapace.width[which(!is.na(x$carapace.width) & x$sex == 2 & x$gonad.colour != 3 & !is.mature(x))]
xp <- x$carapace.width[which(x$carapace.width >= 28 & !is.na(x$carapace.width) & x$sex == 2 & x$gonad.colour == 3 & !is.mature(x))]
xm <- x$carapace.width[which(x$carapace.width >= 36 & !is.na(x$carapace.width) & x$sex == 2 & is.primiparous(x))]  

ti <- table(round(log(xi), 2))
tp <- table(round(log(xp), 2))
tm <- table(round(log(xm), 2))

theta <- c(mu0 = 2.3020407, 
           log.sigma = -2.44048, 
           logit.p.immature  = rep(1, 7), 
           logit.p.pubescent = rep(1, 7),
           logit.p.mature    = rep(1, 7),
           a = 0.87669, b = 0.6857,
           a.pubescent = 0.777, b.pubescent = 0.92,
           a.mature = 1.0919, b.mature = -0.374)

fixed <- c("a", "b", "mu0", "log.sigma", "a.pubescent", "b.pubescent", "a.mature", "b.mature")
fixed <- theta[fixed]
theta <- theta[setdiff(names(theta), names(fixed))]

loglike(theta, xi = ti, fixed = fixed)
loglike(theta, xi = ti, xp = tp, fixed = fixed)
loglike(theta, xi = ti, xp = tp, xm = tm, fixed = fixed)

#theta <- optim(theta, loglike, xi = ti, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
#theta <- optim(theta, loglike, xi = ti, fixed = fixed, control = list(maxit = 5000, trace = 3))$par

theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 1000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, fixed = fixed, control = list(maxit = 5000, trace = 3))$par

theta <- c(theta, fixed)
#fixed <- c("a.pubescent", "b.pubescent")
#fixed <- theta[fixed]
#theta <- theta[setdiff(names(theta), names(fixed))]

theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, control = list(maxit = 1000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, control = list(maxit = 5000, trace = 3))$par
theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, control = list(maxit = 10000, trace = 3))$par

logit.p.immature <- theta[grep("logit.p.immature", names(theta))]
p.immature <- exp(logit.p.immature) / (1 + sum(exp(logit.p.immature)))
p.immature[length(p.immature)+1] <- 1 - sum(p.immature)
mu.immature <- c(theta["mu0"], rep(NA, length(p.immature)-1))
sigma <- exp(theta["log.sigma"])
for (i in 2:length(p.immature)){
   mu.immature[i]    <- theta["a"] * mu.immature[i-1] + theta["b"]
}

logit.p.pubescent <- theta[grep("logit.p.pubescent", names(theta))]
p.pubescent <- exp(logit.p.pubescent) / (1 + sum(exp(logit.p.pubescent)))
p.pubescent[length(p.pubescent)+1] <- 1 - sum(p.pubescent)

logit.p.mature <- theta[grep("logit.p.mature", names(theta))]
p.mature <- exp(logit.p.mature) / (1 + sum(exp(logit.p.mature)))
p.mature[length(p.mature)+1] <- 1 - sum(p.mature)

mu.pubescent <- theta["a.pubescent"] * mu.immature + theta["b.pubescent"]
#mu.pubescent <- 0.78 * mu.immature + 0.88
mu.mature <- theta["a.mature"] * mu.pubescent + theta["b.mature"]

clg()
t <- seq(0, 5, len = 1000)
m <- kronecker(matrix(c(1,2,3)), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, m, 0), 0)
layout(m)
par(mar = c(0,0,0,0))
res <- 2
gbarplot(ti, xlim = c(2.2, 4.4), xaxs = "i")
grid()
vline(mu.immature, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu)){
   v <- v + 0.01 * sum(ti) * p.immature[i] * dnorm(t, mu.immature[i], sigma)
   lines(t, 0.01 * sum(ti) * p.immature[i] * dnorm(t, mu.immature[i], sigma), lwd = 2, col = "blue")
} 
lines(t, v, lwd = 2, col = "red")
box(col = "grey50")

gbarplot(tp, xlim = c(2.2, 4.4), xaxs = "i")
grid()
vline(mu.pubescent, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu)){
  v <- v + 0.01 * sum(tp) * p.pubescent[i] * dnorm(t, mu.pubescent[i], sigma)
  lines(t, 0.01 * sum(tp) * p.pubescent[i] * dnorm(t, mu.pubescent[i], sigma), lwd = 2, col = "blue")
} 
lines(t, v, lwd = 2, col = "red")
box(col = "grey50")

gbarplot(tm, xlim = c(2.2, 4.4), xaxs = "i")
grid()
vline(mu.mature, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu.mature)){
  v <- v + 0.01 * sum(tm) * p.mature[i] * dnorm(t, mu.mature[i], sigma)
  lines(t, 0.01 * sum(tm) * p.mature[i] * dnorm(t, mu.mature[i], sigma), lwd = 2, col = "blue")
} 
lines(t, v, lwd = 2, col = "red")
box(col = "grey50")

