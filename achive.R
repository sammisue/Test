## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
## Reproduce results of Eubank (1997)

## Last Updated: 03/19/2015

## This code is based code written on 03/05/2015, some bugs are fixed.
## Simulation of a_{0.05}=4.18 is done for verification purpose
## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
eubank <- function(K, n, p, T){
    rej_new<- 0                        #initial of rejects of the new data-driven method
    rej_newq <- 0                      #initial of rejects of new method using q
    rej_classical <- 0                 #initial of rejects of classical test
    rej_transform <- 0                 #initial of rejects of transformed classcial test
    Q <- rep(NA, T)                    #vector storage of all qhat's
    p0 <- rep(1/K, K)                  #value of H0, using capital K
    J <- K-1                           #J is always K-1
    q <- seq(0, K-1, by=1)             #q is from 0 to 9
    M <- rep(0, length(q))             #initial M with M(q=0)=0
    coef1 <- (n+1)/(n-1)               #coeficient 1 of M
    coef2 <- 4.18/(n-1)                #coeficient 2 of M, 4.18 is for a_{0.05}
    X <- matrix(NA, J, K)              #matrix of Fourier coefficients
    for (j in 1:J){
        for (k in 1:K){
            X[j, k] <- sqrt(2/K)*cos(j*pi*(k-0.5)/K)
        }
    }
    for (t in 1:T){
        data <- rmultinom(1, n, prob=p)         #generate data using real prob
        phat <- as.vector(data/n)               #\hat{p} is calculated by simulated data
        fhat <- as.matrix((phat - p0)/sqrt(p0)) #matrix of \hat{f(k)}
        b <- X %*% fhat                         #b1-b9 are generated                
        ## classical chi-sqaure test, alpha=0.05
        Xsquare_classical <- n*sum(fhat^2)
        if (Xsquare_classical > qchisq(0.95, K-1)){
            rej_classical <- rej_classical + 1
        }else{
            rej_classical <- rej_classical
        }
        ## transformed chi-square test, alpha=0.05
        Xsquare_transform <- n*sum(b^2)
        if (Xsquare_transform > qchisq(0.95, K-1)){
            rej_transform <- rej_transform + 1
        }else{
            rej_transform <- rej_transform
        }
        ## paper's new chi-square test, alpha=0.05
        v <- (X^2) %*% (phat/p0)                #v_jj for M(q)
        for (qq in 1:J){
            M[qq+1] <- coef1*sum((b[1:qq, 1])^2) - coef2*sum(v[1:qq, 1])
        }
        qhat <- which.max(M) - 1                #q is from 0 to 9
        Q[t] <- qhat
        if (qhat != 0){
            Xsquare_new <- n*sum((b[1:qhat, 1])^2)      #test statistic
            rej_newq <- rej_newq + 1
            if (Xsquare_new > qchisq(0.95, qhat)){
                rej_new <- rej_new+ 1
            }else{
                rej_new <- rej_new
            }
        }else{
            rej_new<- rej_new
            rej_newq <- rej_newq
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
    reject_new <- rej_new/T
    reject_newq <- rej_newq/T
    reject_classical <- rej_classical/T
    reject_transform <- rej_transform/T
    list(reject_new = reject_new, reject_newq = reject_newq, Q = Q,
         reject_classical = reject_classical, reject_transform = reject_transform)
}



## simulation of a_{0.05}=4.18 for verification
K <- 10
a <- 4.18
T <- 5000
count <- 0
for (t in 1:T){
    test <- rep(NA, K-1)
    for (k in 1:(K-1)){
        test[k] <- (1/k)*rchisq(1, df=k)
    }
    if (max(test) >= a){
        count <- count + 1
    }else{
        count <- count
    }
}
alpha <- count/T
print(alpha)




## The following code simulates alpha=0.05
## individual case with T=10,000
K <- 10                                 #capital K
T <- 10000                              #Paper used 10,000 iterations
beta <- seq(0, 0.07, by=0.01)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        allp[i, k] <- 0.1+beta[i]*(k-5.5)/10  #set real prob for simulation alpha
    }
}

## Simulate the power of original and transformed classical chi-square
## test using both my code and R build-in function
power_myclass <- rep(NA, length(beta))
power_R <- rep(NA, length(beta))
power_mytransform <- rep(NA, length(beta))
for (kk in 1: length(beta)){
    p <- allp[kk,]
    n <- 75
    ## n <- 150
    ## my code
    sim <- eubank(K, n, p, T)
    power_myclass[kk] <- sim$reject_classical
    power_mytransform[kk] <- sim$reject_transform
    ## R build-in function "chisq.test()"
    rej_R <- 0
    for (i in 1:T){
        data <- as.vector(rmultinom(1, n, prob=p))   #generate data using real prob
        Rtest <- chisq.test(data, p=rep(1/K, K))
        if (Rtest$p.value < 0.05){
            rej_R <- rej_R + 1
        }else{
            rej_R <- rej_R
        }
    }
    power_R[kk] <- rej_R/T
}

pdf("1_classical.pdf", height=5, width=7.5)
par(mfrow=c(1, 1))
plot(beta, power_myclass, ylim=c(0,1), ylab="Power", col=1, lty=1, type="b",
     main="Comparison of Classical Chi-Square, n=")
points(beta, power_R, col=2, lty=2, type="b")
points(beta, power_mytransform, col=3, lty=3, type="b")
legend("topright", c("My Classical", "R Function", "My Transform"), col=c(1,2,3),
       lty=1:3)
dev.off()


## Simulate the power of the new data-driven method
power_new <- rep(NA, length(beta))
power_newq <- rep(NA, length(beta))
power_class <- rep(NA, length(beta))
for (kk in 1:length(beta)){
    p <- allp[kk,]
    n <- 75
    ## n <- 150
    sim <- eubank(K, n, p, T)
    power_new[kk] <- sim$reject_new
    power_newq[kk] <- sim$reject_newq
    power_class[kk] <- sim$reject_classical
}

pdf("3_new.pdf", height=5, width=7.5)
par(mfrow=c(1, 1))
plot(beta, power_new, ylim=c(0,1), ylab="Power", col=1, lty=1, type="b",
     main="Comparison of New and Classical Chi-Square, n=, a_(0.05)=4.18")
points(beta, power_newq, col=2, lty=2, type="b")
points(beta, power_class, col=3, lty=3, type="b")
legend("topright", c("New Method", "New Method using q", "Classical"), col=c(1,2,3),
       lty=1:3)
dev.off()
############################# END ###########################################
