## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
## Reproduce results of Eubank (1997)

## Last Updated: 03/05/2015

## This is the initial code I wrote, it has some bugs.
## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
eubank <- function(K, n, p, T){
    rej <- 0                           #initial value of number of rejections
    rej_classical <- 0                 #initial value of rejects of classical test
    Q <- rep(NA, T)                    #vector storage of all qhat's
    p0 <- rep(1/K, K)                  #value of H0, using capital K
    J <- K-1                           #J is always K-1
    q <- seq(0, K-1, by=1)             #q is from 0 to 9
    M <- rep(0, length(q))             #initial M with M(q=0)=0
    coef1 <- (n+1)/(n-1)               #coeficient 1 of M
    coef2 <- 2/(n-1)                   #coeficient 2 of M
    X <- matrix(NA, J, K)              #matrix of Fourier coefficients
    for (j in 1:J){
        for (k in 1:K){
            X[j, k] <- sqrt(2/K)*cos(j*pi*(k-0.5)/K)
        }
    }
    for (t in 1:T){
        data <- rmultinom(1, n, prob=p)         #generate data using real prob
        phat <- as.vector(data/n)               #\hat{p} is calculated by simulated data
        ## classical chi-sqaure test
        Xsquare_classical <- n*sum(((phat-p0)^2)/p0)
        if (Xsquare_classical > qchisq(0.95, K-1)){
            rej_classical <- rej_classical + 1
        }else{
            rej_classical <- rej_classical
        }
        ## paper's chi-square test
        fhat <- as.matrix((phat - p0)/sqrt(p0)) #matrix of \hat{f(k)}
        b <- X %*% fhat                         #b1-b9 are generated
        v <- (X^2) %*% (phat/p0)                #v_jj for M(q)
        for (q in 1:J){
            M[q+1] <- coef1*sum(b[1:q, 1]) - coef2*sum(v[1:q, 1])
        }
        qhat <- which.max(M) - 1                #q is from 0 to 9
        Q[t] <- qhat
        if (qhat != 0){
            Xsquare <- n*sum((b^2)[1:qhat, 1])      #test statistic
            critical <- qchisq(0.95, qhat)       #critial chi-sq value for alpha=0.05
            ## critical <- 4.18     #paper's critical value, which doesn't work well
            if (Xsquare > critical){
                rej <- rej + 1
            }else{
                rej <- rej
            }
        }else{
            rej <- rej
        }
        ## show progress
        if (t==1 | t==T){
            cat ("Iteration:", 100*(t/T), "%, ", t, "of", T, "\n")
            print(Sys.time())
            flush.console()
        }
        if (t %% 5000==0){
            cat ("Iteration:", 100*(t/T), "%, ", t, "of", T, "\n")
            flush.console()
        }
    }
    reject <- rej/T; reject_classical <- rej_classical/T
    list(reject = reject, Q = Q, reject_classical = reject_classical)
}



## The following code simulates alpha=0.05
## individual case with T=600000
K <- 10                                 #capital K
n <- 150                                #total number draws
p <- rep(1/K, K)                        #set real prob for simulation alpha
T <- 600000
sim <- eubank(K, n, p, T)
sink("individual.Rout",  append=FALSE)
## print simlulated type I error
print(sim$reject); print(sim$reject_classical)
## print tables for all qhat's
print(table(sim$Q)); print(round(table(sim$Q)/T, 3))
sink()

## plot K vs. Alpha and n vs. Alpha
## K vs. Alpha, when n=150
n <- 150                                #total number draws
T <- 600000
testK <- seq(3, 50, by=1)
myreject1 <- rep(NA, length(testK))
for (i in 1:length(testK)){
    K <- testK[i]; p <- rep(1/K, K)
    sim <- eubank(K, n, p, T)
    myreject1[i] <- sim$reject    
}
## n vs. Alpha when K=10
K <- 10                                 #capital K
p <- rep(1/K, K)                        #set real prob for simulation alpha
T <- 600000
testn <- seq(50, 500, by=5)
myreject2 <- rep(NA, length(testn))
for (i in 1:length(testn)){
    n <- testn[i]
    sim <- eubank(K, n, p, T)
    myreject2[i] <- sim$reject
}
## Plot
pdf("alpha.pdf", height=5, width=10)
par(mfrow=c(1, 2))
plot(testK, myreject1, type="o", main="K vs. Alpha, when n=150", xlab="K", ylab="Alpha")
plot(testn, myreject2, type="o", main="n vs. Alpha, when K=10", xlab="n", ylab="Alpha")
dev.off()

## The following code simulates power when alpha=0.05 with n=75, 150 and K=10
K <- 10                                 #capital K
T <- 600000
beta <- seq(0, 0.07, by=0.01)
testp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        testp[i, k] <- 0.1+beta[i]*(k-5.5)/10  #set real prob for simulation alpha
    }
}
## beta vs. power.
n <- 75                                #total number draws
myreject3 <- rep(NA, length(beta))
for (i in 1:length(beta)){
    p <- as.vector(testp[i,])
    sim <- eubank(K, n, p, T)
    myreject3[i] <- sim$reject
}
n <- 150                                #total number draws
myreject4 <- rep(NA, length(beta))
for (i in 1:length(beta)){
    p <- as.vector(testp[i,])
    sim <- eubank(K, n, p, T)
    myreject4[i] <- sim$reject
}
## plot
pdf("power.pdf", height=5, width=10)
par(mfrow=c(1, 2))
plot(beta, myreject3, type="o", main="Beta vs. Power, n=75", xlab="Beta", ylab="Power")
plot(beta, myreject4, type="o", main="Beta vs. Power, n=150", xlab="Beta", ylab="Power")
dev.off()
############################# END ###########################################
