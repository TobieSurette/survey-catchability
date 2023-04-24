library(gulf.data)
library(gulf.graphics)

x <- read.scsbio(1998:2022)

loglike <- function(theta, xi, xp, xm, fixed){
   if (!missing(fixed)) theta <- c(theta, fixed)
   
   # Parse parameter vector:
   logit.p.immature <- theta[grep("logit.p.immature", names(theta))]
   p.immature <- exp(logit.p.immature) / (1 + sum(exp(logit.p.immature)))
   p.immature[length(p.immature)+1] <- 1 - sum(p.immature)
   mu.immature <- c(theta["mu0"], rep(NA, length(p.immature)-1))
   sigma <- exp(theta["log.sigma"])
   for (i in 2:length(p.immature)) mu.immature[i]    <- theta["a"] * mu.immature[i-1] + theta["b"]
   names(mu.immature) <- as.roman(3 + (1:length(mu.immature)))
   names(p.immature) <- names(mu.immature)
   
   # Immature mixture density:
   ni <- as.numeric(xi)
   xi <- as.numeric(names(xi))
   vi <- rep(0, length(xi))
   for (i in 1:length(p.immature)) vi <- vi + p.immature[i] * dnorm(xi, mu.immature[i], sigma)
   ll <- - sum(ni * log(vi)) 
   
   # Pubescent:
   if (!missing(xp)){
      logit.p.pubescent <- theta[grep("logit.p.pubescent", names(theta))]
      p.pubescent <- exp(logit.p.pubescent) / (1 + sum(exp(logit.p.pubescent)))
      p.pubescent[length(p.pubescent)+1] <- 1 - sum(p.pubescent)
      
      mu.pubescent <- theta["a.pubescent"] * mu.immature[as.character(as.roman(7:9))] + theta["b.pubescent"]
      names(mu.pubescent) <- as.roman(8:10)
      names(p.pubescent) <- as.roman(8:10)
      np <- as.numeric(xp)
      xp <- as.numeric(names(xp)) 
      vp <- rep(0, length(xp))
      sigma.pubescent <- rep(sigma, length(mu.pubescent))
      if ("log.sigma.pubescent" %in% names(theta)) sigma.pubescent <- sigma + exp(theta["log.sigma.pubescent"] * mu.pubescent)
      for (i in 1:length(p.pubescent)) vp <- vp + p.pubescent[i] * dnorm(xp, mu.pubescent[i], sigma.pubescent[i])
      ll <- ll - sum(np * log(vp)) 
   }

   # Mature:
   if (!missing(xm)){
      logit.p.mature <- theta[grep("logit.p.mature", names(theta))]
      p.mature <- exp(logit.p.mature) / (1 + sum(exp(logit.p.mature)))
      p.mature[length(p.mature)+1] <- 1 - sum(p.mature)
     
      mu.mature <- theta["a.mature"] * mu.pubescent[as.character(as.roman(8:10))] + theta["b.mature"]
      names(mu.mature) <- as.roman(9:11)
      names(p.mature) <- as.roman(9:11)
      
      nm <- as.numeric(xm)
      xm <- as.numeric(names(xm)) 
      vm <- rep(0, length(xm))

      sigma.mature <- rep(sigma, length(mu.mature))
      if ("log.sigma.mature" %in% names(theta)) sigma.mature <- sigma + exp(theta["log.sigma.mature"] * mu.mature)
     
      for (i in 1:length(p.mature)) vm <- vm + p.mature[i] * dnorm(xm, mu.mature[i], sigma.mature[i])
      ll <- ll - sum(nm * log(vm)) 
   }
   
   return(ll)
}

# Prepare data:
xi <- x$carapace.width[which(!is.na(x$carapace.width) & x$sex == 2 & x$gonad.colour != 3 & !is.mature(x))]
xp <- x$carapace.width[which(x$carapace.width >= 28 & !is.na(x$carapace.width) & x$sex == 2 & x$gonad.colour == 3 & !is.mature(x))]
xm <- x$carapace.width[which(x$carapace.width >= 36 & !is.na(x$carapace.width) & x$sex == 2 & is.primiparous.scsbio(x))]  

ti <- table(round(log(xi), 2))
tp <- table(round(log(xp), 2))
tm <- table(round(log(xm), 2))

theta <- c(mu0 = 2.2957724, 
           log.sigma = -2.5701350, 
           log.sigma.pubescent = -1,
           log.sigma.mature = -1,
           logit.p.immature  = c(-1.9025520, 0.3302821, 2.1653939, 3.3263875, 2.9821234), 
           logit.p.pubescent = c(3.3735471, 2.9042755),
           logit.p.mature    = c(8.2518782, 8.2026561),
           a = 0.8667117, b = 0.7225728,            
           a.pubescent = 0.7069271, b.pubescent = 1.3312277,          
           a.mature = 0.934210, b.mature = 0.4103975)
         
fixed <- c("a", "b", "mu0", "log.sigma", "log.sigma.pubescent", "log.sigma.mature", "a.pubescent", "b.pubescent", "a.mature", "b.mature")
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
for (i in 1:20) theta <- optim(theta, loglike, xi = ti, xp = tp, xm = tm, control = list(maxit = 5000, trace = 3))$par
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
names(mu.immature) <- as.roman((1:length(mu.immature))+3)
names(mu.immature) <- as.roman(3 + (1:length(mu.immature)))
names(p.immature) <- names(mu.immature)
   
logit.p.pubescent <- theta[grep("logit.p.pubescent", names(theta))]
p.pubescent <- exp(logit.p.pubescent) / (1 + sum(exp(logit.p.pubescent)))
p.pubescent[length(p.pubescent)+1] <- 1 - sum(p.pubescent)
      
mu.pubescent <- theta["a.pubescent"] * mu.immature[as.character(as.roman(7:9))] + theta["b.pubescent"]
names(mu.pubescent) <- as.roman(8:10)
names(p.pubescent) <- as.roman(8:10)
     
sigma.pubescent <- rep(sigma, length(mu.pubescent))
if ("log.sigma.pubescent" %in% names(theta)) sigma.pubescent <- sigma + exp(theta["log.sigma.pubescent"] * mu.pubescent)
      
logit.p.mature <- theta[grep("logit.p.mature", names(theta))]
p.mature <- exp(logit.p.mature) / (1 + sum(exp(logit.p.mature)))
p.mature[length(p.mature)+1] <- 1 - sum(p.mature)
     
mu.mature <- theta["a.mature"] * mu.pubescent[as.character(as.roman(8:10))] + theta["b.mature"]
names(mu.mature) <- as.roman(9:11)
names(p.mature) <- as.roman(9:11)

sigma.mature <- rep(sigma, length(mu.mature))
if ("log.sigma.mature" %in% names(theta)) sigma.mature <- sigma + exp(theta["log.sigma.mature"] * mu.mature)
 
clg()
cols <- c("palegreen3", "skyblue4")
xlim <- c(log(10), log(80))
t <- seq(0, 5, len = 1000)
m <- kronecker(matrix(c(1,2,3)), matrix(1, ncol = 5, nrow = 5))
m <- rbind(0, cbind(0, m, 0), 0, 0)
layout(m)
par(mar = c(0,0,0,0))
res <- 2
gbarplot(ti, xlim = xlim, xaxs = "i", xaxt = "n", col = "grey90", border = "grey80")
vline(log(seq(0, 100, by = 5)), lty = "dashed", col = "grey65")
#vline(mu.immature, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu.immature)){
   v <- v + 0.01 * sum(ti) * p.immature[i] * dnorm(t, mu.immature[i], sigma)
   d <- 0.01 * sum(ti) * p.immature[i] * dnorm(t, mu.immature[i], sigma)
   lines(t, d, lwd = 2, col = cols[1])
   if (p.immature[i] > 0.01) text(t[which.max(d)], max(d), names(mu.immature)[i], pos = 3, font = 2, cex = 1.0)
} 
lines(t, v, lwd = 2, col = cols[2])
mtext("Immature", 4, 1.0)
box(col = "grey50")

gbarplot(tp, xlim = xlim, xaxs = "i", xaxt = "n", col = "grey90", border = "grey80")
vline(log(seq(0, 100, by = 5)), lty = "dashed", col = "grey65")
#vline(mu.pubescent, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu.pubescent)){
   d <- 0.01 * sum(tp) * p.pubescent[i] * dnorm(t, mu.pubescent[i], sigma.pubescent[i])  
   v <- v + d
   lines(t, d, lwd = 2, col = cols[1])
   if (p.pubescent[i] > 0.01) text(t[which.max(d)], max(d), names(mu.pubescent)[i], pos = 3, font = 2, cex = 1.0)
} 
lines(t, v, lwd = 2, col = cols[2])
mtext("Pubescent", 4, 1.0)
mtext("Density (#/km2)", 2, 2.5, cex = 1.25)
box(col = "grey50")

gbarplot(tm, xlim = xlim, xaxs = "i", xaxt = "n", col = "grey90", border = "grey80")
vline(log(seq(0, 100, by = 5)), lty = "dashed", col = "grey65")
#vline(mu.mature, lty = "dashed", lwd = 2, col = "blue")
v <- rep(0, length(t))
for (i in 1:length(mu.mature)){
   d <- 0.01 * sum(tm) * p.mature[i] * dnorm(t, mu.mature[i], sigma.mature[i])
   v <- v + d
   lines(t, d, lwd = 2, col = cols[1])
   if (p.mature[i] > 0.01) text(t[which.max(d)], max(d), names(mu.mature)[i], pos = 3, font = 2, cex = 1.0)
} 
lines(t, v, lwd = 2, col = cols[2])
box(col = "grey50")
axis(1, at = log(seq(0, 100, by = 5)), labels = seq(0, 100, by = 5), lty = "dashed", col = "grey65")
mtext("Carapace width (mm)", 1, 2.5, cex = 1.25)
mtext("Primiparous", 4, 1.0)

