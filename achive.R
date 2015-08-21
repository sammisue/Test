## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
## Reproduce results of Eubank (1997)

## Last Updated: 03/26/2015

## This code is based on modification of the code of 03/19/2015
## 1. I don't use R build-in function "chisq.test" for the classical
##    Chi-square test, because it is slower than my manual code.
## 2. Add W statistic in order to reproduce figures on Eubank's paper.
## 3. The simulation of W under H0 is done which is consistent with
##    the results on Eubank's paper
## 4. It makes sense now that W statistic depends on qhat, and the
##    other new method depends on qhat_{alpha}. Notice that "qhat" and
##    "qhat_{alpha}" are different.
## 5. We can input value of alpha and a_alpha in the function now.
## 6. qhat and qhat_{alpha} are both recorded in the function.
## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
eubank <- function(K, n, p, T, alpha, a_alpha){
    rej_new<- 0                        #initial of rejects of the new data-driven method
    rej_newq <- 0                      #initial of rejects of new method using q
    rej_W <- 0                         #initial of rejects of method using W
    rej_classical <- 0                 #initial of rejects of classical test
    rej_transform <- 0                 #initial of rejects of transformed classcial test
    Q <- rep(NA, T)                    #vector storage of all qhat's
    Q_alpha <- rep(NA, T)              #vector storage of all qhat_alpha's
    p0 <- rep(1/K, K)                  #value of H0, using capital K
    J <- K-1                           #J is always K-1
    q <- seq(0, K-1, by=1)             #q is from 0 to 9
    M <- rep(0, length(q))             #initial M with M(q=0)=0
    X <- matrix(NA, J, K)              #matrix of Fourier coefficients
    for (j in 1:J){
        for (k in 1:K){
            X[j, k] <- sqrt(2/K)*cos(j*pi*(k-0.5)/K)
        }
    }
    ## find critical value for W
    TT <- 100000
    data_W0 <- rmultinom(TT, n, p=p0)
    phat_W0 <- data_W0/n
    fhat_W0 <- apply(phat_W0, 2, function(x) (x-p0)/sqrt(p0))
    X_W0 <- matrix(NA, J, K)              #matrix of Fourier coefficients
    for (j in 1:J){
        for (k in 1:K){
            X_W0[j, k] <- sqrt(2/K)*cos(j*pi*(k-0.5)/K)
        }
    }
    b_W0 <- X_W0 %*% fhat_W0                         #b1-b9 are generated
    v_W0 <- apply(phat_W0/p0, 2, function(x) (X_W0^2) %*% x) #v_jj for M(q) or M_{alpha}(q)
    M_W0 <- matrix(0, K, TT)
    for (tt in 1:TT){
        for (qq in 1:J){
            M_W0[qq+1, tt] <- ((n+1)/(n-1))*sum((b_W0[1:qq, tt])^2) -
                (2/(n-1))*sum(v_W0[1:qq, tt])
        }
    }
    qhat_W0 <- apply(M_W0, 2, function(x) (which.max(x)-1))
    W0 <- rep(NA, TT)
    for (i in 1:TT){
        if (qhat_W0[i] != 0){
            W0[i] <- (n*sum((b_W0[1:qhat_W0[i], i])^2) - qhat_W0[i])/sqrt(2*qhat_W0[i])
        }else{
            W0[i] <- 0
        }
    }
    critical_W <- quantile(W0, probs=(1-alpha)) #critical_W is found in here
    for (t in 1:T){
        data <- as.vector(rmultinom(1, n, prob=p)) #generate data using real prob
        phat <- data/n                          #\hat{p} is calculated by simulated data
        fhat <- as.matrix((phat - p0)/sqrt(p0)) #matrix of \hat{f(k)}
        b <- X %*% fhat                         #b1-b9 are generated                
        ## classical chi-sqaure test, alpha=0.05
        Xsquare_classical <- n*sum(fhat^2)
        if (Xsquare_classical > qchisq((1-alpha), K-1)){
            rej_classical <- rej_classical + 1
        }else{
            rej_classical <- rej_classical
        }
        ## transformed chi-square test, alpha=0.05
        Xsquare_transform <- n*sum(b^2)
        if (Xsquare_transform > qchisq((1-alpha), K-1)){
            rej_transform <- rej_transform + 1
        }else{
            rej_transform <- rej_transform
        }
        ## paper's new chi-square test for W using qhat
        v <- (X^2) %*% (phat/p0)                #v_jj for M(q) or M_{alpha}(q)
        for (qq in 1:J){
            M[qq+1] <- ((n+1)/(n-1))*sum((b[1:qq, 1])^2) - (2/(n-1))*sum(v[1:qq, 1])
        }
        qhat<- which.max(M) - 1                #q is from 0 to 9
        Q[t] <- qhat
        if (qhat != 0){
            Xsquareq <- n*sum((b[1:qhat, 1])^2)
            W <- (Xsquareq - qhat)/sqrt(2*qhat)      #test statistic W
        }else{
            W <- 0
        }
        if (W > critical_W){        #2.99 is for 0.05 level, 2.3 is for 0.1 level
            rej_W <- rej_W + 1      #reject using test statistic "W"
        }else{
            rej_W <- rej_W
        }
        ## paper's new chi-square test for qhat_alpha
        for (qq in 1:J){
            M[qq+1] <- ((n+1)/(n-1))*sum((b[1:qq, 1])^2) - (a_alpha/(n-1))*sum(v[1:qq, 1])
        }
        qhat_alpha<- which.max(M) - 1                #q is from 0 to 9
        Q_alpha[t] <- qhat_alpha
        if (qhat_alpha != 0){
            rej_newq <- rej_newq + 1    #reject using \hat{q}            
            Xsquare_new <- n*sum((b[1:qhat_alpha, 1])^2)      #test statistic
            if (Xsquare_new > qchisq((1-alpha), qhat_alpha)){
                rej_new <- rej_new + 1   #reject using test statistic "Xsquare_new"
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
    reject_W <- rej_W/T
    reject_classical <- rej_classical/T
    reject_transform <- rej_transform/T
    list(reject_new = reject_new, reject_newq = reject_newq, reject_W = reject_W,
         Q = Q, Q_alpha = Q_alpha, 
         reject_classical = reject_classical, reject_transform = reject_transform)
}



## The following code simulates alpha=0.05
## individual case with T=10,000

## Alternative (21)
K <- 10                                 #capital K
beta <- seq(0, 0.14, by=0.01)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        #set real prob for simulation alpha
        allp[i, k] <- 0.1+beta[i]*(k-5.5)/10  
    }
}

## Alternative (22)
K <- 10                                 #capital K
beta <- seq(0, 0.1, by=0.01)
j <- 2                                  #set j=2 for Figure 3
## j <- 4                                  #set j=4 for Figure 4
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        #set real prob for simulation alpha
        allp[i, k] <- 0.1+beta[i]*cos(j*pi*(k-0.5)/10)  
    }
}
print(allp)

## Alternative (23)
K <- 10                                 #capital K
beta <- seq(0.6, 1.4, by=0.1)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        #set real prob for simulation alpha
        allp[i, k] <- pnorm(beta[i]*qnorm(k/10))-pnorm(beta[i]*qnorm((k-1)/10))  
    }
}
print(allp)




## Simulate the power of the new data-driven method, alpha=0.05
n <- 75                                 #75 draws for multinomial distribution
## n <- 150                                #150 draws for multinomial distribution
T <- 10000                              #Paper used 10,000 iterations
alpha <- 0.05; a_alpha <- 4.18          #for level 0.05
## alpha <- 0.1; a_alpha <- 3.22           #for level 0.1
power_new <- rep(NA, length(beta))
power_newq <- rep(NA, length(beta))
power_W <- rep(NA, length(beta))
power_class <- rep(NA, length(beta))
qhat_rec <- rep(NA, length(beta))       #mean value of qhat for every 10,000 iterations
for (kk in 1:length(beta)){
    cat ("beta=", beta[kk])        
    p <- allp[kk,]
    sim <- eubank(K, n, p, T, alpha, a_alpha)
    power_new[kk] <- sim$reject_new
    power_newq[kk] <- sim$reject_newq
    power_W[kk] <- sim$reject_W
    power_class[kk] <- sim$reject_classical
    qhat_rec[kk] <- mean(sim$Q)
}

## For alternative (21)
## pdf("f2b.pdf", height=5, width=7.5)
## par(mfrow=c(1, 1))
plot(beta, power_new, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0,0.14), col=1, lty=1, type="l", lwd=1,
     main=paste("Comparison of New and Classical Chi-Square, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0, 0.14, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_newq, col=2, lty=2, type="l", lwd=3)
points(beta, power_class, col=3, lty=3, type="l", lwd=4)
points(beta, power_W, col=4, lty=4, type="l", lwd=2)
legend("topleft", c("New Method", "New Method using q", "Classical", "W"),
       col=1:4, lty=1:4, lwd=c(1,3,4,2))
box()
## dev.off()


## For alternative (22)
## pdf("f3b.pdf", height=7.5, width=7.5)
## par(mfrow=c(1, 1))
plot(beta, power_new, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0,0.07), col=1, lty=1, type="l", lwd=1,
     main=paste("Power, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0, 0.07, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_newq, col=2, lty=2, type="l", lwd=3)
points(beta, power_class, col=3, lty=3, type="l", lwd=4)
points(beta, power_W, col=4, lty=4, type="l", lwd=2)
legend("bottomright", c("New Method", "New Method using q", "Classical", "W"),
       col=1:4, lty=1:4, lwd=c(1,3,4,2))
box()
## dev.off()

## pdf("f4c.pdf", height=7.5, width=7.5)
## par(mfrow=c(1, 1))
plot(beta, qhat_rec, ylim=c(0,5), ylab="Average qhat",
     xlim=c(0,0.07), col=1, lty=1, type="l", lwd=1,
     main=paste("qhat vs. beta, a_(0.05)=4.18, n=",n))
## dev.off()


## For alternative (23)
## pdf("f5b.pdf", height=7.5, width=11.25)
## par(mfrow=c(1, 1))
plot(beta, power_new, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0.4, 1.4), col=1, lty=1, type="l", lwd=1,
     main=paste("Power, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0.4, 1.4, by=0.1))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_newq, col=2, lty=2, type="l", lwd=3)
points(beta, power_class, col=3, lty=3, type="l", lwd=4)
points(beta, power_W, col=4, lty=4, type="l", lwd=2)
legend("bottomleft", c("New Method", "New Method using q", "Classical", "W"),
       col=1:4, lty=1:4, lwd=c(1,3,4,2))
box()
## dev.off()

#################################################################################


#################################################################################
## simulate W under H0 with 100,000 iterations
## method 1, this method is NOT correct!!!!!! Thus it does NOT work good!!!!!!!!
K <- 10
J <- K-1
T <- 100000
critical_W <- matrix(NA, 4, J)
W0 <- matrix(0, T, K)
for (qhat in 1:9){
    chisq_qhat <- rchisq(T, df=qhat)       #Chisq_qhat under H0
    W0[,qhat+1] <- (chisq_qhat-qhat)/sqrt(2*qhat)        #distr of W under H0
    ## critical_W[,qhat] <- quantile(W0, probs=0.95)
}

quantile(as.vector(W0), prob=0.95)
print(round(critical_W, 3))
## write.table(round(critical_W, 3), file = "cao.csv", sep=",", append = FALSE)


## method 2, the results of this method are consistent with Eubank's paper!!!!!!
critical_W <- c(10, 10)
while (critical_W[1] > 3){                        #just want to get 2.99 and 2.3
    K <- 10
    J <- K-1
    TT <- 100000
    p0 <- rep(1/K, K)
    n <- 75
    ## n <- 150
    data <- rmultinom(TT, n, p=p0)
    phat <- data/n
    fhat <- apply(phat, 2, function(x) (x-p0)/sqrt(p0))
    X <- matrix(NA, J, K)              #matrix of Fourier coefficients
    for (j in 1:J){
        for (k in 1:K){
            X[j, k] <- sqrt(2/K)*cos(j*pi*(k-0.5)/K)
        }
    }
    b <- X %*% fhat                         #b1-b9 are generated
    v <- apply(phat/p0, 2, function(x) (X^2) %*% x)      #v_jj for M(q) or M_{alpha}(q)
    M <- matrix(0, K, TT)
    for (tt in 1:TT){
        for (qq in 1:J){
            M[qq+1, tt] <- ((n+1)/(n-1))*sum((b[1:qq, tt])^2) - (2/(n-1))*sum(v[1:qq, tt])
        }
    }
    qhat <- apply(M, 2, function(x) (which.max(x)-1))
    W0 <- rep(NA, TT)
    for (i in 1:TT){
        if (qhat[i] != 0){
            W0[i] <- (n*sum((b[1:qhat[i], i])^2) - qhat[i])/sqrt(2*qhat[i])
        }else{
            W0[i] <- 0
        }
    }
    critical_W <- quantile(W0, probs=c(0.95, 0.90))
    ## write.table(round(critical_W, 3), file = "cao.csv", sep=",", append = FALSE)
}
print(critical_W)
#################################################################################


#################################################################################
## Simulation of a_{\alpha} using method in Eubank's 1997 paper, which
## doesn't work good, since a_{0.05} is approx. 4.5, not 4.18
M <- 100                                #create 100 emprical distributions
a <- matrix(NA, M, 5)
for (m in 1:M){
    ## One empirical distribution
    K <- 10                                 #K=10 categories in multinomial distribution
    T <- 5000                               #T=5000 values for one empirical distribution
    stat <- rep(NA, T)
    for (t in 1:T){
        ## find maximum for an individual run
        test <- rep(NA, K-1)
        for (k in 1:(K-1)){
            test[k] <- (1/k)*rchisq(1, df=k)
        }
        stat[t] <- max(test)
    }
    ## quantile(stat, probs = seq(0.9, 1, 0.01)) #check quantiles for one empirical distr.
    a[m,] <- quantile(stat, probs=c(0.99, 0.95, 0.9, 0.8, 0.71))
}
apply(a, 2, mean)                          #find average of a
## hist(a_0.05)                            #plot all a_{0.05}
############################# END ###########################################
