## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
## Applying Eubank's method (1997) to Weighted multinomial data.

## Last Updated: 05/05/2015

## This code is based on modification of the code of 04/02/2015
## 1. Equivalence of Sum version and Matrix version of Chi-squared
##    test statistic is verified.
## 2. "rej_newq" is now renamed as "rej_q", which is tested by "Q_alpha".
## 3. "rej_new" is removed because it is equivalent to "rej_q".
## 4. "Q" is now renamed as "Q_W" because it is the value of q for W.
## 5. The generation of weighted multinomial data is now integrated in
##    the "eubank" function
## 6. Generating a batch of "pdf" plots using loop is now available in
##    "eubank_run.R". A simple sample code is included
## 7. Rao's 1st and 2nd order chi-squared tests are added, but with concerns.
## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+

################## Generating Weighted Multinomial Data #####################
## Generates a multinomial single frame population with clusters of size psusize
## input:
## numpsu = number of psus
## psusize = size of clusters
## rho = intraclass correlation coefficient, generates dependence within psus
gensample <- function(p, numpsu, psusize, rho){
    ## set.seed(1346)
    ## if (sum(p) != 1) stop("sum of domain probabilities must equal 1")
    n <- numpsu*psusize  #population size, total number of obs in the population
    pmat <- p
    qcatpop <- rep(1,n)
    for (i in 1:numpsu){
        permcat <- sample(1:length(p), length(p))
        ## Randomly permute the categories so the ordering
        ## does not determine the clustering
        pi <- pmat[permcat]
        ## Find probit quantities
        probcut <- cumsum(pi[1:(length(p)-1)]) #(length(p)-1) must be in parenthesis
        qcut <- qnorm(probcut)
        ## We generate z_{ij} = \alpha_i + \eps_{ij},
        ## where \alpha_i ~ N(0,rho) and \eps_{ij} ~ N(0,1-rho).
        ## Then z_{ij} are clustered N(0,1) random variables.
        zval <- rep(rnorm(1, 0, rho), psusize) + rnorm(psusize, 0, 1-rho)
        qcat <- rep(permcat[length(p)], psusize) 
        ## Now use probit to assort obs into categories
        for (k in (length(p)-1):1){     #(length(p)-1) must be in parenthesis
            qcat[zval < qcut[k]] <- permcat[k]
        }
        qcatpop[( (i-1)*psusize+1 ):(i*psusize)] <- qcat
    }
    weightedData <- rep(0, length(p))
    for (uu in 1:length(p)){
        weightedData[uu] <- length(qcatpop[qcatpop==uu])
    }
    return(weightedData)    #return to number of trials in each category
}
############################################################################


#################### Function of Eubank's paper ############################
## Input:
## K = number of categories in multinomial distribution
## numpsu = number of psus
## psusize = size of clusters
## n = numpsu * psusize, total number of trials in a multinomial data
## rho = ICC in each cluster
## p = vector of probabilities for each category
## T = number of simulations
## alpha = level of significance of the test
## a_alpha = a value which is paired with alpha for q_hat
## Please switch between Weighed and Unweighted data for particular goal.
eubank <- function(K, numpsu, psusize, rho, p, T, alpha, a_alpha){
    n <- numpsu*psusize                #total number of trials
    rej_q <- 0                         #initial of rejects of new method using q
    rej_W <- 0                         #initial of rejects of method using W
    rej_classical <- 0                 #initial of rejects of classical test
    rej_transform <- 0                 #initial of rejects of transformed classcial test
    rej_1st <- 0                       #initial of rejects of 1st order test
    rej_2nd <- 0                       #initial of rejects of 2nd order test
    Q_W<- rep(NA, T)                    #vector storage of all qhat's for W
    Q_alpha <- rep(NA, T)              #vector storage of all qhat_alpha's for rej_q
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
    TT <- 10000                        #Enbank's paper used 100,000 for this


    ## data_W0 <- rmultinom(TT, n, p=p0)   #generate data under H0
    data_W0 <- matrix(NA, length(p0), TT)
    for (tt in 1:TT){
        data_W0[,tt] <- gensample(p0, numpsu, psusize, rho) #generate data under H0
    }
    


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
    ## Loop starts here
    for (t in 1:T){

        
        ## data <- as.vector(rmultinom(1,n,prob=p))#generate unweighted data by real prob
        data <- gensample(p, numpsu, psusize, rho) #generate weighted data by real prob
        
        
        phat <- data/n                          #\hat{p} is calculated by simulated data
        fhat <- as.matrix((phat - p0)/sqrt(p0)) #matrix of \hat{f(k)}
        b <- X %*% fhat                         #b1-b9 are generated                
        ## classical chi-sqaure test, at level of alpha
        Xsquare_classical <- n*sum(fhat^2)
        if (Xsquare_classical > qchisq((1-alpha), (K-1))){
            rej_classical <- rej_classical + 1
        }else{
            rej_classical <- rej_classical
        }
        ## transformed chi-square test, at level of alpha
        Xsquare_transform <- n*sum(b^2)
        if (Xsquare_transform > qchisq((1-alpha), (K-1))){
            rej_transform <- rej_transform + 1
        }else{
            rej_transform <- rej_transform
        }
        ## Rao chi-square test with the 1st order correction, at level of alpha
        ## d <- ((phat*(1-phat))/n)/((phat*(1-phat))/n)
        d <- 1
        delta_dot <- sum((phat*(1-phat)*d/p0))/(K-1)
        Xsquare_1st <- Xsquare_classical/delta_dot
        if (Xsquare_1st > qchisq((1-alpha), (K-1))){
            rej_1st <- rej_1st + 1
        }else{
            rej_1st <- rej_1st
        }
        ## Rao chi-square test with the 2nd order correction, at level of alpha
        nV <- matrix(NA, (K-1), (K-1))          #Covariance matrix is (K-1) by (K-1)
        for (i in 1:(K-1)){
            for (j in 1:(K-1)){
                if (i==j){
                    nV[i, j] <- phat[i]*(1-phat[i])
                }else{
                    nV[i, j] <- -phat[i]*phat[j]
                }
            }
        }
        for (i in 1:(K-1)){
            for (j in 1:(K-1)){
                temp <- nV/sqrt(p0[i]*p0[j])
            }
        }
        asquare <- sum(temp^2)/((K-1)*delta_dot^2) - 1
        Xsquare_2nd <- Xsquare_1st/(1+asquare)
        df_2nd <- (K-1)/(1+asquare)
        if (Xsquare_2nd > qchisq((1-alpha), df_2nd)){
            rej_2nd <- rej_2nd + 1
        }else{
            rej_2nd <- rej_2nd
        }
        ## paper's new chi-square test for W using qhat
        v <- (X^2) %*% (phat/p0)                #v_jj for M(q) or M_{alpha}(q)
        for (qq in 1:J){
            M[qq+1] <- ((n+1)/(n-1))*sum((b[1:qq, 1])^2) - (2/(n-1))*sum(v[1:qq, 1])
        }
        qhat<- which.max(M) - 1                #q is from 0 to 9
        Q_W[t] <- qhat
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
            rej_q <- rej_q + 1    #reject using q_hat            
        }else{
            rej_q <- rej_q
        }
        ## show progress
        if (t==1 | t==T){
            cat("Iteration:", 100*(t/T), "%, ", t, "of", T, "\n")
            print(Sys.time())
            flush.console()
        }
        if (t %% 5000==0){
            cat("Iteration:", 100*(t/T), "%, ", t, "of", T, "\n")
            flush.console()
        }
    }
    reject_q <- rej_q/T
    reject_W <- rej_W/T
    reject_classical <- rej_classical/T
    reject_transform <- rej_transform/T
    reject_1st <- rej_1st/T
    reject_2nd <- rej_2nd/T
    list(reject_q = reject_q, reject_W = reject_W, Q_W= Q_W, Q_alpha = Q_alpha, 
         reject_classical = reject_classical, reject_transform = reject_transform,
         reject_1st = reject_1st, reject_2nd = reject_2nd)
}

#################################################################################
##
## Test for Eubank's methods, as well as 1st and 2nd order of the chi-square test
##
rho <- 0.8
K <- 10
beta <- seq(0, 0.14, by=0.01)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        ## set real prob for simulation alpha
        allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K
    }
}
print(allp)
print(apply(allp, 1, sum))
numpsu <- 15                            #number of clusters
psusize <- 50                            #number of ssus in each cluster
n <- numpsu*psusize                     #equivalent to n=15*5=75
T <- 10000
alpha <- 0.05; a_alpha <- 4.18
power_q <- rep(NA, length(beta))
power_W <- rep(NA, length(beta))
power_class <- rep(NA, length(beta))
power_1st <- rep(NA, length(beta))
power_2nd <- rep(NA, length(beta))
for (kk in 1:length(beta)){
    cat ("beta =", beta[kk], "\n")        
    p <- allp[kk,]
    sim <- eubank(K, numpsu, psusize, rho, p, T, alpha, a_alpha)
    power_q[kk] <- sim$reject_q
    power_W[kk] <- sim$reject_W
    power_class[kk] <- sim$reject_classical
    power_1st[kk] <- sim$reject_1st
    power_2nd[kk] <- sim$reject_2nd
}

plot(beta, power_q, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0,0.14), col=1, lty=1, type="l", lwd=2,
     main=paste("Comparison of New and Classical Chi-Square, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0, 0.14, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_class, col=2, lty=2, type="l", lwd=1)
points(beta, power_W, col=3, lty=3, type="l", lwd=3)
points(beta, power_1st, col=4, lty=4, type="l", lwd=1.5)
points(beta, power_2nd, col=5, lty=5, type="l", lwd=2.5)
legend("bottomright", c("New Method using q", "Classical", "W", "1st", "2nd"),
       col=1:5, lty=1:5, lwd=c(2, 1, 3, 1.5, 2.5))
box()


##
## Test Unweighted data for K=3:10, with 1st, 2nd and Eubank's method
## Need to modify "eubank_current.R" to generate Unweighted data.
##
## Run the function and generating corresponding pdf plots using loop
pdf(paste("unweighted", ".pdf", sep=""), height=20, width=7.5)
par(mfrow=c(4, 2))                 #10 plots for 10 values of rho
for (K in 3:10){                        #number of categories in multinomial
    beta <- seq(0, 0.14, by=0.01)
    allp <- matrix(NA, length(beta), K)
    for (i in 1:length(beta)){
        for (k in 1:K){
            ## set real prob for simulation alpha
            allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K
        }
    }
    print(allp); print(apply(allp, 1, sum))    #check if the sum of probabilities is 1
    rho <- 0.9999       #we don't need rho for unweighted, just an input for my function
    ## Simulate the power of the new data-driven method, alpha=0.05
    numpsu <- 15                            #number of clusters
    psusize <- 50                            #number of ssus in each cluster
    n <- numpsu*psusize                     #equivalent to n=15*5=75
    ## rho <- 0.3                              #ICC in cluster
    ## n <- 75                                 #75 draws for multinomial distribution
    ## n <- 150                                #150 draws for multinomial distribution
    T <- 10000                              #Paper used 10,000 iterations
    alpha <- 0.05; a_alpha <- 4.18          #for level 0.05
    ## alpha <- 0.1; a_alpha <- 3.22           #for level 0.1
    power_q <- rep(NA, length(beta))
    power_W <- rep(NA, length(beta))
    power_class <- rep(NA, length(beta))
    power_1st <- rep(NA, length(beta))
    power_2nd <- rep(NA, length(beta))
    qhat_rec <- rep(NA, length(beta))  #mean value of qhat for every 10,000 iterations
    for (kk in 1:length(beta)){
        cat ("K =", K, "rho =", rho, "beta =", beta[kk], "\n")        
        p <- allp[kk,]
        sim <- eubank(K, numpsu, psusize, rho, p, T, alpha, a_alpha)
        power_q[kk] <- sim$reject_q
        power_W[kk] <- sim$reject_W
        power_class[kk] <- sim$reject_classical
        power_1st[kk] <- sim$reject_1st
        power_2nd[kk] <- sim$reject_2nd            
        qhat_rec[kk] <- mean(sim$Q_W)
    }
    plot(beta, power_q, axes=FALSE, ylim=c(0,1), ylab="Power",
         xlim=c(0,0.14), col=1, lty=1, type="l", lwd=2,
         main=paste("a_(0.05)=4.18, n=",n, "," , "K=", K))
    axis(side=1, at=seq(0, 0.14, by=0.01))
    axis(side=2, at=seq(0, 1, by=0.05))
    points(beta, power_class, col=2, lty=2, type="l", lwd=1)
    points(beta, power_W, col=3, lty=3, type="l", lwd=3)
    points(beta, power_1st, col=4, lty=4, type="l", lwd=1.5)
    points(beta, power_2nd, col=5, lty=5, type="l", lwd=2.5)
    legend("bottomright", c("New Method using q", "Classical", "W", "1st", "2nd"),
           col=1:5, lty=1:5, lwd=c(2, 1, 3, 1.5, 2.5))
    box()
}
dev.off()


##
## Test Weighted data for K=3:10, with 1st, 2nd and Eubank's method
## Need to modify "eubank_current.R" to generate Weighted data.
##
## Run the function and generating corresponding pdf plots using loop
for (K in 3:10){                        #number of categories in multinomial

    pdf(paste("K=",K, ".pdf", sep=""), height=50, width=7.5)
    par(mfrow=c(10, 1))                 #10 plots for 10 values of rho

    beta <- seq(0, 0.14, by=0.01)
    allp <- matrix(NA, length(beta), K)
    for (i in 1:length(beta)){
        for (k in 1:K){
            ## set real prob for simulation alpha
            allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K
        }
    }
    print(allp); print(apply(allp, 1, sum))    #check if the sum of probabilities is 1
    for (rho in c(0.01, seq(0.1, 0.9, by=0.1))){
        ## Simulate the power of the new data-driven method, alpha=0.05
        numpsu <- 15                            #number of clusters
        psusize <- 50                            #number of ssus in each cluster
        n <- numpsu*psusize                     #equivalent to n=15*5=75
        ## rho <- 0.3                              #ICC in cluster
        ## n <- 75                                 #75 draws for multinomial distribution
        ## n <- 150                                #150 draws for multinomial distribution
        T <- 10000                              #Paper used 10,000 iterations
        alpha <- 0.05; a_alpha <- 4.18          #for level 0.05
        ## alpha <- 0.1; a_alpha <- 3.22           #for level 0.1
        power_q <- rep(NA, length(beta))
        power_W <- rep(NA, length(beta))
        power_class <- rep(NA, length(beta))
        power_1st <- rep(NA, length(beta))
        power_2nd <- rep(NA, length(beta))
        qhat_rec <- rep(NA, length(beta))  #mean value of qhat for every 10,000 iterations
        for (kk in 1:length(beta)){
            cat ("K =", K, "rho =", rho, "beta =", beta[kk], "\n")        
            p <- allp[kk,]
            sim <- eubank(K, numpsu, psusize, rho, p, T, alpha, a_alpha)
            power_q[kk] <- sim$reject_q
            power_W[kk] <- sim$reject_W
            power_class[kk] <- sim$reject_classical
            power_1st[kk] <- sim$reject_1st
            power_2nd[kk] <- sim$reject_2nd            
            qhat_rec[kk] <- mean(sim$Q_W)
        }
        plot(beta, power_q, axes=FALSE, ylim=c(0,1), ylab="Power",
             xlim=c(0,0.14), col=1, lty=1, type="l", lwd=2,
             main=paste("a_(0.05)=4.18, n=",n, "," , "K=", K, ",", "rho=", rho))
        axis(side=1, at=seq(0, 0.14, by=0.01))
        axis(side=2, at=seq(0, 1, by=0.05))
        points(beta, power_class, col=2, lty=2, type="l", lwd=1)
        points(beta, power_W, col=3, lty=3, type="l", lwd=3)
        points(beta, power_1st, col=4, lty=4, type="l", lwd=1.5)
        points(beta, power_2nd, col=5, lty=5, type="l", lwd=2.5)
        legend("bottomright", c("New Method using q", "Classical", "W", "1st", "2nd"),
               col=1:5, lty=1:5, lwd=c(2, 1, 3, 1.5, 2.5))
        box()
    }
    dev.off()
}
#################################################################################


#################################################################################
#################### Above code are based on Code below #########################
## The following code simulates when alpha=0.05
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
print(allp)

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
numpsu <- 15                            #15 clusters
psusize <- 50                            #5 ssus in each cluster
n <- numpsu*psusize                     #equivalent to n=15*5=75
rho <- 0.1                              #ICC in cluster
## n <- 75                                 #75 draws for multinomial distribution
## n <- 150                                #150 draws for multinomial distribution
T <- 10000                              #Paper used 10,000 iterations
alpha <- 0.05; a_alpha <- 4.18          #for level 0.05
## alpha <- 0.1; a_alpha <- 3.22           #for level 0.1
power_q <- rep(NA, length(beta))
power_W <- rep(NA, length(beta))
power_class <- rep(NA, length(beta))
qhat_rec <- rep(NA, length(beta))       #mean value of qhat for every 10,000 iterations
for (kk in 1:length(beta)){
    cat ("beta =", beta[kk], "\n")        
    p <- allp[kk,]
    sim <- eubank(K, numpsu, psusize, rho, p, T, alpha, a_alpha)
    power_q[kk] <- sim$reject_q
    power_W[kk] <- sim$reject_W
    power_class[kk] <- sim$reject_classical
    qhat_rec[kk] <- mean(sim$Q_W)
}


## For alternative (21)
## pdf("f2a.pdf", height=5, width=7.5)
## par(mfrow=c(1, 1))
plot(beta, power_q, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0,0.14), col=1, lty=1, type="l", lwd=2,
     main=paste("Comparison of New and Classical Chi-Square, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0, 0.14, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_class, col=2, lty=2, type="l", lwd=1)
points(beta, power_W, col=3, lty=3, type="l", lwd=3)
legend("topleft", c("New Method using q", "Classical", "W"),
       col=1:3, lty=1:3, lwd=c(2, 1, 3))
box()
## dev.off()


## For alternative (22)
## pdf("f3b.pdf", height=7.5, width=7.5)
## par(mfrow=c(1, 1))
plot(beta, power_q, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0,0.07), col=1, lty=1, type="l", lwd=2,
     main=paste("Power, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0, 0.07, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_class, col=2, lty=2, type="l", lwd=1)
points(beta, power_W, col=3, lty=3, type="l", lwd=3)
legend("bottomright", c("New Method using q", "Classical", "W"),
       col=1:3, lty=1:3, lwd=c(2, 1, 3))
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
plot(beta, power_q, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0.4, 1.4), col=1, lty=1, type="l", lwd=2,
     main=paste("Power, a_(0.05)=4.18, n=",n))
axis(side=1, at=seq(0.4, 1.4, by=0.1))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_class, col=2, lty=2, type="l", lwd=1)
points(beta, power_W, col=3, lty=3, type="l", lwd=3)
legend("bottomleft", c("New Method using q", "Classical", "W"),
       col=1:3, lty=1:3, lwd=c(2, 1, 3))
box()
## dev.off()
#################################################################################


#################################################################################
##
## Sample code for generating a lot of "pdf" plots using a loop
## 
for (k in 3:5){
    pdf(paste("K=", k, ".pdf", sep=""), height=10, width=7.5)
    par(mfrow=c(2, 1))
    for (mu in 1:2){
        plot(rnorm(100, mu, k), main=paste("mu=", mu))
    }
    dev.off()
}
#################################################################################


#################################################################################
## Simulation of a_{\alpha} using Eubank's method in his 1992 paper,
## which works much better than the method in his 1997 paper. This
## simulation can give exact a_{0.05}=4.18.
K <- 15   #number of terms for summation, K should be infinity, but even 10 works good.
alpha <- c(0.01, 0.05, 0.10, 0.20, 0.29) #alpha level
a_alpha <- seq(0, 10, by=0.01)            #create possible a_{alpha} for search
test <- matrix(NA, length(a_alpha), K)   #matrix that store each term for each a_{alpha}
for (i in 1:length(a_alpha)){
    for (k in 1:K){
        test[i, k] <- ((1-pchisq(k*a_alpha[i], df=k))/k) #probability of each term
    }
}
stat <- matrix(NA, length(a_alpha), length(alpha))
for (j in 1:length(alpha)){
    stat[,j] <- abs(exp(-apply(test, 1, sum))-(1-alpha[j])) #fomula in Eubank's 1992 paper
}
index <- apply(stat, 2, function(x) which.min(x))#find index which has minimum error
cat("alpha=", alpha, "\n", "a_alpha=", a_alpha[index], "\n")
#################################################################################


#################################################################################
## Verify the Sum version and Matrix version of X^2 in chi-squared
## test are equivalent.
n <- 75                                 #total number of trials in multinomial
K <- 3                                  #total number of categories in multinomial
p0 <- rep(1/K, K)                       #H0, all probabilities are the same
set.seed(1)
data <- rmultinom(1, n, prob=p0)
phat <- data/n
## Sum version
Xsquare_sum <- n*sum(((phat-p0)^2)/p0)
## Matrix version
phat <- phat[1:(K-1)]                   #Matrix version has only (K-1) dimension
p0 <- p0[1:(K-1)]                       #Matrix version has only (K-1) dimension
P0 <- matrix(NA, (K-1), (K-1))          #Covariance matrix is (K-1) by (K-1)
for (i in 1:(K-1)){
    for (j in 1:(K-1)){
        if (i==j){
            P0[i, j] <- p0[i]*(1-p0[i])
        }else{
            P0[i, j] <- -p0[i]*p0[j]
        }
    }
}
Xsquare_matrix <- n*(t((phat-p0))%*%solve(P0)%*%(phat-p0))
## print both of versions
print(Xsquare_sum); print(Xsquare_matrix)
#################################################################################



#################################################################################
## Test: Is the generated weighted multinomial data equivalent
## between my revised version and Dr. Lu's original version?
gensampleorg<-function(p, ntilde,psusize, rho){
    ## Generates a multinomial single frame population with clusters of size psusize
    ## Input:
    ## ntilde<-number of psu
    ## psusize = s-ize of clusters 
    ## rho = intraclass correlation coefficient, generates dependence within psus
    ## set.seed(1346)
    ## if (sum(p) != 1) stop("sum of domain probabilities must equal 1")
    N<-ntilde*psusize  #population size
    ## numpsu <- ntilde
    pmat <- p
    qcatpop <-rep(1,N)
    for (i in 1:ntilde) {
        permcat <- sample(1:length(p),length(p))
        ## Randomly permute the categories so the ordering
        ## does not determine the clustering
        pi <- pmat[permcat]
        ## Find probit quantities
        probcut <- cumsum(pi[1:(length(p)-1)])
        qcut <- qnorm(probcut)
        ## We generate z_{ij} = \alpha_i + \eps_{ij},
        ## where \alpha_i ~ N(0,rho) and \eps_{ij} ~ N(0,1-rho).
        ## Then z_{ij} are clustered N(0,1) random variables.
        zval <- rep(rnorm(1,0,rho),psusize)+rnorm(psusize,0,1-rho)
        qcat <- rep(permcat[length(p)],psusize) 
        ## Now use probit to assort into categories
        for (k in (length(p)-1):1) {
            qcat[zval < qcut[k]] <- permcat[k]
        }
        qcatpop[( (i-1)*psusize+1 ): (i*psusize)] <- qcat
    }  
    list(qcatpop=qcatpop,N=N)
}

K <- 10
p <- rep(1/K, K)
numpsu <- 15; psusize <- 50; rho <- 0.9
for (i in 1:100){
    set.seed(i)
    org <- table(gensampleorg(p, numpsu, psusize, rho)$qcatpop)
    set.seed(i)
    my <- gensample(p, numpsu, psusize, rho)
    print(org); print(my)
}
################################# END ###########################################
