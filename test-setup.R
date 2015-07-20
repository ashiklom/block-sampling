# Test block sampling on a simple cubic model
model <- function(b,x) b[1]*x^3 + b[2]*x^2 + b[3]*x + b[4]
coefs <- c(4,3,2,1)
x <- seq(1, 100, by=0.1)
obs <- model(coefs, x) + rnorm(length(x), 0, 200*x)
plot(x,obs, pch=".")

loglike <- function(b){
    error <- model(b,x) - obs
    logl <- sum(dnorm(error, 0, 0.5))
    return(logl)
}
