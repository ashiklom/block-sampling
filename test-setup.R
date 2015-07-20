# Test block sampling on a simple cubic model
library(minpack.lm)
library(PEcAnRTM)
library(MASS)
true.vals <- c("N"=1.7, "Cab"=32, "Cw"=0.01, "Cm"=0.008)
obs <- prospect(true.vals, 4)

ngibbs <- 10000
adapt <- 100
accept <- 0
target <- 0.44
adjmin <- 0.1
pmin <- c(1,0,0,0)
#inits <- true.vals
inits <- c(2, 40, 0.001, 0.001)
initsd <- inits*0.01
Jump <- diag(initsd)

result <- matrix(NA, nrow=ngibbs, ncol=4)
for(n in 1:ngibbs){
    if(n > adapt && n %% adapt == 1){
        if(accept == 0){
            rescale <- diag(rep(adjmin, 4))
            Jump <- rescale %*% Jump %*% rescale
        }
        else {
            arate <- accept/adapt
            region <- (n-adapt):(n-1)
            stdev <- apply(result[region,], 2, sd)
            stdev[is.na(stdev)] <- initsd[is.na(stdev)]
            adjust <- max(arate/target, adjmin)
            rescale <- diag(stdev * adjust)
            corr <- cor(result[region,])
            if(any(is.na(corr))) corr <- diag(rep(1,4))
            Jump <- rescale %*% corr %*% rescale
        }
        print(diag(Jump))
        accept <- 0
    }
    guess <- mvrnorm(1, inits, Jump)
    if(any(guess < pmin)){
        result[n,] <- inits
        next
    }
    guess.obs <- prospect(guess, 4)
    guess.error <- guess.obs - obs
    guess.like <- sum(dnorm(guess.error, 0, 0.1, log=TRUE))
    prev.obs <- prospect(inits, 4)
    prev.error <- prev.obs - obs
    prev.like <- sum(dnorm(prev.error, 0, 0.1, log=TRUE))
    a <- exp(guess.like - prev.like)
    if(a > runif(1)){
        accept <- accept + 1
        inits <- guess
    }
    result[n,] <- inits
}

par(mfrow=c(4,2))
for(i in 1:4){
    plot(result[,i], type='l')
    abline(h = true.vals[i], col="red")
    plot(density(result[-5000:0,i]))
    abline(v = true.vals[i], col="red")
}

