
instar.moments <- function(x, species = 2526, sex = 2, maturity){
      
   if (species == 2526 & sex == 2){
       # Female snow crab, survey data 1998 to 2022:
       theta <- c(a = 0.89594187, 
                  b = 0.62825172,
                  mu0 = 2.30276852,
                  log.sigma = -2.38146892,
                  log.sigma.pubescent = -3,
                  log.sigma.mature = -1.22662952,
                  a.pubescent = 0.65044384,
                  b.pubescent = 1.56508817,
                  a.mature =  1.00392668,  
                  b.mature = 0.15810622)
 
       # Calculate instar means:
       mu.immature <- theta["mu0"]
       for (i in 2:6) mu.immature[i] <- theta["a"] * mu.immature[i-1] + theta["b"]
       names(mu.immature) <- as.roman(3 + (1:length(mu.immature)))
   
       # Pubescent:
       mu.pubescent <- theta["a.pubescent"] * mu.immature[as.character(as.roman(7:9))] + theta["b.pubescent"]
       names(mu.pubescent) <- as.roman(8:10)
   
       # Mature:
       mu.mature <- theta["a.mature"] * mu.pubescent[as.character(as.roman(8:10))] + theta["b.mature"]
       names(mu.mature) <- as.roman(9:11)

       # Instar errors:
       sigma <- exp(theta["log.sigma"])
       sigma.pubescent <- sigma
       sigma.pubescent <- sigma + exp(theta["log.sigma.pubescent"] * mu.pubescent)
       sigma.mature <- sigma
       sigma.mature <- sigma + exp(theta["log.sigma.mature"] * mu.mature)
       
       # Build result variable:
       r <- list(immature = list(mu = mu.immature, sigma = rep(sigma, length(mu.immature))),
                 pubescent = list(mu = mu.pubescent, sigma = sigma.pubescent),
                 mature = list(mu = mu.mature, sigma = sigma.mature))
       
       if (!missing(maturity)) return(r[maturity]) else return(r)
    }
}

instar.loglike <- function(theta, xi, xp, xm, fixed){
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

y <- x
x <- xi

tuesday morning 

instar.default <- function(x, species = 2526, sex = 2, maturity = "immature", log = TRUE){
   if (log) x <- log(x)
   m <- instar.moments(species = 2526, sex = 2)
   
   # Immature:
   r <- NULL
   if (maturity == "immature") for (i in 1:length(m$immature$mu)) r <- cbind(r, dnorm(x, m$immature$mu[i], m$immature$sigma[i]))
   colnames(r) <- names(mu.immature)
   p.immature <- apply(r, 2, sum) / sum(r)
   
   # Pubescent:
   r <- NULL
   if (maturity == "pubescent") for (i in 1:length(m$pubescent$mu)) r <- cbind(r, dnorm(x, m$pubescent$mu[i], m$pubescent$sigma[i]))
   colnames(r) <- names(mu.pubescent)
   p.pubescent <- apply(r, 2, sum) / sum(r)  
   
   # Mature:
   r <- NULL
   if (maturity == "mature") for (i in 1:length(m$mature$mu)) r <- cbind(r, dnorm(x, m$mature$mu[i], m$mature$sigma[i]))
   colnames(r) <- names(mu.mature)
   p.mature <- apply(r, 2, sum) / sum(r)  
   
   # Calculate the instar classification for each observation:
   
   
   
}
