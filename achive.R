## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+
## Applying Eubank's method (1997) to Weighted multinomial data.

## Last Updated: 9/1/2015

## This code is based on modification of the code of 08/24/2015

## 1. Updated the simulation of equation to get value of a_alpha to
##    the version which takes advantage of eigenvalues. This updated
##    simulation is designed to choose a_alpha for the weighted
##    multinomial data. The simulation is based on the equation in
##    Eubank's 1992 paper.
## 2. Now, the code for selection of a_alpha is in a seperate file.
## +--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+--+



######## Function to get \hat{p} for each repeated generated sample ###########
get_count_sample <- function(K, sample){
    count_sample <- rep(NA, K)
    for (uu in 1:K){
        count_sample[uu] <- length(sample[sample==uu])
    }
    return(count_sample)
}
############################################################################



################ Test for Rao's 1st and 2nd coorection #####################
rao <- function(K, numpsu, psusize, rho, p, T, alpha){
    n <- numpsu*psusize                #total number of trials
    rej_classical <- 0                 #initial of rejects of classical test    
    rej_1st <- 0                       #initial of rejects of 1st order test
    rej_2nd <- 0                       #initial of rejects of 2nd order test
    p0 <- rep(1/K, K)                  #value of H0, using capital K
    ## ##################################################################
    ## find covariance matrix P0 under H0 and SRS
    P0 <- matrix(NA, (K-1), (K-1))          #Covariance matrix is (K-1) by (K-1)
    for (i in 1:(K-1)){
        for (j in 1:(K-1)){
            if (i==j){
                P0[i, j] <- p0[i]*(1-p0[i])
            }else{
                P0[i, j] <- -p0[i]*p0[j] #negative
            }
        }
    }
    P0 <- P0/n
    ## ##################################################################
    ## using repeated samples to get \hat{V} and eigenvalues
    samplenum <- T[3]                       #number of samples
    phat_sample <- matrix(NA, samplenum, length(p))    
    for (ii in 1:samplenum){

        ## ------------------------------------------------ ##
        ## Unweighted option
        ## temp_data <- as.vector(rmultinom(1, n, prob=p))
        ## Weighted option 1
        temp_data <- gensample(p, numpsu, psusize, rho)$weightedData
        ## Weighted option 2
        ## popn <- gensample(p, numpsu, psusize, rho)$popn #weighted population
        ## temp_data <- get_count_sample(K, popn)
        ## ------------------------------------------------ ##
        
        temp_data[temp_data==0] <- 0.01 #make empty cell a small count, 0.01        
        phat_sample[ii,] <- temp_data/n
    }
    pbar_sample <- t(replicate(samplenum, apply(phat_sample, 2, mean)))
    ## Covariance matrix \hat{V} is obtained, which is (K-1)x(K-1)
    V <- (t((phat_sample-pbar_sample)[, -K])%*%(phat_sample-pbar_sample)[, -K])/(samplenum-1)
    ## Get 1st order correction using estimation of eigenvalues
    v_kk <- sum(V)                      #calculate variance of \hat{p_k} based on \hat{V}
    V_diag <- c(as.vector(diag(V)), v_kk) #vector of all variances of \hat{p_1}, ..., \hat{p_K}
    delta_dot<- n*sum(V_diag/p0)/(K-1)  #delta.dot for 1st order correction

    ## ----------------------------------------------------------------------------- ##
    ## Get full estimated covariance matrix for 2nd order correction
    ## full estimated covariance matrix Vfull is (K x K), the following 2 are the same
    Vfull <- rbind(cbind(V, -apply(V, 1, sum)), c(-apply(V, 1, sum), v_kk))
    ## Vfull <- (t((phat_sample-pbar_sample))%*%(phat_sample-pbar_sample))/(samplenum-1)
    ## ----------------------------------------------------------------------------- ##

    temp_matrix<- matrix(NA, K, K) #the temporary matrix is to calculate numerator of asquare
    for (i in 1:K){
        for (j in 1:k){
            temp_matrix[i, j] <- (Vfull[i, j]^2)/(p0[i]*p0[j])
        }
    }
    ## 2nd order correction will be (1+asquare)
    asquare <- (n^2)*sum(temp_matrix)/((K-1)*delta_dot^2) - 1
    ## check value of "asquare" and eigenvalues
    eigenvalues <- eigen(solve(P0) %*% V)$values
    cat("asquare=", asquare, "\n", "eigenvalues=", eigenvalues, "\n")
    ## #################################################################
    for (t in 1:T[1]){

        ## ------------------------------------------------ ##        
        ## Unweighted option
        ## data <- as.vector(rmultinom(1, n, prob=p))#real unweighted data
        ## Weighted option
        data <- gensample(p, numpsu, psusize, rho)$weightedData  #real weighted data
        ## ------------------------------------------------ ##
        
        data[data==0] <- 0.01 #make empty cell a small count, 0.01
        phat <- data/n                          #\hat{p} is calculated by simulated data
        fhat <- as.matrix((phat - p0)/sqrt(p0)) #matrix of \hat{f(k)}        
        ## ##########################################################
        ## Classical chi-sqaure test, at level of alpha
        Xsquare_classical <- n*sum(fhat^2)
        if (Xsquare_classical > qchisq((1-alpha), (K-1))){
            rej_classical <- rej_classical + 1
        }else{
            rej_classical <- rej_classical
        }
        ## ##########################################################
        ## Rao chi-square test with the 1st order correction, at level of alpha
        Xsquare_1st <- Xsquare_classical/delta_dot
        if (Xsquare_1st > qchisq((1-alpha), (K-1))){
            rej_1st <- rej_1st + 1
        }else{
            rej_1st <- rej_1st
        }
        ## ##########################################################        
        ## Rao chi-square test with the 2nd order correction, at level of alpha
        Xsquare_2nd <- Xsquare_1st/(1+asquare)
        df_2nd <- (K-1)/(1+asquare)
        if (Xsquare_2nd > qchisq((1-alpha), df_2nd)){
            rej_2nd <- rej_2nd + 1
        }else{
            rej_2nd <- rej_2nd
        }
        ## ##########################################################
    }
    reject_classical <- rej_classical/T[1]
    reject_1st <- rej_1st/T[1]
    reject_2nd <- rej_2nd/T[1]
    list(reject_classical = reject_classical,
         reject_1st = reject_1st, reject_2nd = reject_2nd)
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
## T = (number of simulations, iteration for W, number of repeated samples for \hat{V})
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
    Q_W<- rep(NA, T[1])                    #vector storage of all qhat's for W
    Q_alpha <- rep(NA, T[1])              #vector storage of all qhat_alpha's for rej_q
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
    ## #############################################################
    ## find covariance matrix P0 under H0 and SRS
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
    P0 <- P0/n
    ## #############################################################
    ## using repeated samples to get \hat{V} and eigenvalues
    samplenum <- T[3]                       #number of samples
    phat_sample <- matrix(NA, samplenum, length(p))    
    for (ii in 1:samplenum){
        
        ## ------------------------------------------------ ##        
        ## Unweighted option
        ## temp_data <- as.vector(rmultinom(1, n, prob=p))
        ## Weighted option 1
        temp_data <- gensample(p, numpsu, psusize, rho)$weightedData
        ## Weighted option 2
        ## popn <- gensample(p, numpsu, psusize, rho)$popn #weighted population
        ## temp_data <- get_count_sample(K, popn)
        ## ------------------------------------------------ ##
        
        temp_data[temp_data==0] <- 0.01 #make empty cell a small count, 0.01        
        phat_sample[ii,] <- temp_data/n
    }
    pbar_sample <- t(replicate(samplenum, apply(phat_sample, 2, mean)))
    ## Covariance matrix \hat{V} is obtained, which is (K-1)x(K-1)
    V <- (t((phat_sample-pbar_sample)[, -K])%*%(phat_sample-pbar_sample)[, -K])/(samplenum-1)
    ## Get 1st order correction using estimation of eigenvalues
    v_kk <- sum(V)                      #calculate variance of \hat{p_k} based on \hat{V}
    V_diag <- c(as.vector(diag(V)), v_kk) #vector of all variances of \hat{p_1}, ..., \hat{p_K}
    delta_dot<- n*sum(V_diag/p0)/(K-1)  #delta.dot for 1st order correction

    ## ----------------------------------------------------------------------------- ##    
    ## Get full estimated covariance matrix for 2nd order correction
    ## full estimated covariance matrix Vfull is (K x K), the following 2 are the same
    Vfull <- rbind(cbind(V, -apply(V, 1, sum)), c(-apply(V, 1, sum), v_kk))
    ## Vfull <- (t((phat_sample-pbar_sample))%*%(phat_sample-pbar_sample))/(samplenum-1)
    ## ----------------------------------------------------------------------------- ##    
    
    temp_matrix<- matrix(NA, K, K) #the temporary matrix is to calculate numerator of asquare
    for (i in 1:K){
        for (j in 1:k){
            temp_matrix[i, j] <- (Vfull[i, j]^2)/(p0[i]*p0[j])
        }
    }
    ## 2nd order correction will be (1+asquare)
    asquare <- (n^2)*sum(temp_matrix)/((K-1)*delta_dot^2) - 1
    ## #################################################################
    ## find critical value for W
    TT <- T[2]                        #Enbank's paper used 100,000 for this
    ## data_W0 <- rmultinom(TT, n, p=p0)   #generate data under H0
    data_W0 <- matrix(NA, length(p0), TT)
    for (tt in 1:TT){
        ## generate data under H0        
        data_W0[,tt] <- gensample(p0, numpsu, psusize, rho)$weightedData 
    }
    data_W0[data_W0==0] <- 0.01         #make empty cell a small count, 0.01
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
    ################################################################    
    ## Loop starts here
    for (t in 1:T[1]){
        
        ## ------------------------------------------------ ##
        ## Unweighted option
        ## data <- as.vector(rmultinom(1, n, prob=p))#real unweighted data
        ## Weighted option
        data <- gensample(p, numpsu, psusize, rho)$weightedData  #real weighted data
        ## ------------------------------------------------ ##
        
        data[data==0] <- 0.01 #make empty cell a small count, 0.01
        phat <- data/n                          #\hat{p} is calculated by simulated data
        fhat <- as.matrix((phat - p0)/sqrt(p0)) #matrix of \hat{f(k)}
        b <- X %*% fhat                         #b1-b9 are generated
        ## ##########################################################
        ## classical chi-sqaure test, at level of alpha
        Xsquare_classical <- n*sum(fhat^2)
        if (Xsquare_classical > qchisq((1-alpha), (K-1))){
            rej_classical <- rej_classical + 1
        }else{
            rej_classical <- rej_classical
        }
        ## ##########################################################        
        ## transformed chi-square test, at level of alpha
        Xsquare_transform <- n*sum(b^2)
        if (Xsquare_transform > qchisq((1-alpha), (K-1))){
            rej_transform <- rej_transform + 1
        }else{
            rej_transform <- rej_transform
        }
        ## ##########################################################
        ## Rao chi-square test with the 1st order correction, at level of alpha
        Xsquare_1st <- Xsquare_classical/delta_dot
        if (Xsquare_1st > qchisq((1-alpha), (K-1))){
            rej_1st <- rej_1st + 1
        }else{
            rej_1st <- rej_1st
        }
        ## ##########################################################        
        ## Rao chi-square test with the 2nd order correction, at level of alpha
        Xsquare_2nd <- Xsquare_1st/(1+asquare)
        df_2nd <- (K-1)/(1+asquare)
        if (Xsquare_2nd > qchisq((1-alpha), df_2nd)){
            rej_2nd <- rej_2nd + 1
        }else{
            rej_2nd <- rej_2nd
        }
        ## ##########################################################
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
        ## ##########################################################        
        ## show progress
        if (t==1 | t==T[1]){
            cat("Iteration:", 100*(t/T[1]), "%, ", t, "of", T[1], "\n")
            print(Sys.time())
            flush.console()
        }
        if (t %% 5000==0){
            cat("Iteration:", 100*(t/T[1]), "%, ", t, "of", T[1], "\n")
            flush.console()
        }
    }
    reject_q <- rej_q/T[1]
    reject_W <- rej_W/T[1]
    reject_classical <- rej_classical/T[1]
    reject_transform <- rej_transform/T[1]
    reject_1st <- rej_1st/T[1]
    reject_2nd <- rej_2nd/T[1]
    list(reject_q = reject_q, reject_W = reject_W, Q_W = Q_W, Q_alpha = Q_alpha, 
         reject_classical = reject_classical, reject_transform = reject_transform,
         reject_1st = reject_1st, reject_2nd = reject_2nd)
}
############################################################################




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
    list(popn = qcatpop, weightedData = weightedData)
}
############################################################################    






## a_table.R
#################################################################################
numpsu <- 50                            #number of psus
psusize <- 10                           #number of ssus in each psu
rho <- 0.3                              #correlation ("ICC") among ssus in each psu
K <- 5                                  #number of categories in multinomial
p0 <- rep(1/K, K)                  #value of prob under H0, using capital K
p <- p0                            #real p is p0, for probability of type I error
n <- numpsu*psusize                #total number of trials

## find covariance matrix P0 under H0 and SRS
P0 <- matrix(NA, (K-1), (K-1))          #Covariance matrix is (K-1) by (K-1)
for (i in 1:(K-1)){
    for (j in 1:(K-1)){
        if (i==j){
            P0[i, j] <- p0[i]*(1-p0[i])
        }else{
            P0[i, j] <- -p0[i]*p0[j] #negative
        }
    }
}
P0 <- P0/n

## using repeated samples to get \hat{V} and eigenvalues
samplenum <- 10000                       #number of samples
phat_sample <- matrix(NA, samplenum, length(p))    
for (ii in 1:samplenum){
    
    ## ------------------------------------------------ ##
    ## Unweighted option
    ## temp_data <- as.vector(rmultinom(1, n, prob=p))
    ## Weighted option 1
    temp_data <- gensample(p, numpsu, psusize, rho)$weightedData
    ## Weighted option 2
    ## popn <- gensample(p, numpsu, psusize, rho)$popn #weighted population
    ## temp_data <- get_count_sample(K, popn)
    ## ------------------------------------------------ ##
    
    temp_data[temp_data==0] <- 0.01 #make empty cell a small count, 0.01        
    phat_sample[ii,] <- temp_data/n
}
pbar_sample <- t(replicate(samplenum, apply(phat_sample, 2, mean)))
## Covariance matrix \hat{V} is obtained, which is (K-1)x(K-1)
V <- (t((phat_sample-pbar_sample)[, -K])%*%(phat_sample-pbar_sample)[, -K])/(samplenum-1)
## Get 1st order correction using estimation of eigenvalues
v_kk <- sum(V)                      #calculate variance of \hat{p_k} based on \hat{V}
V_diag <- c(as.vector(diag(V)), v_kk) #vector of all variances of \hat{p_1}, ..., \hat{p_K}



## It seems these 3 are the same
delta_dot<- n*sum(V_diag/p0)/(K-1)  #delta.dot for 1st order correction
delta_dot1 <- sum(diag(solve(P0) %*% V))/(K-1)
delta_dot2 <- sum(eigen(solve(P0) %*% V)$values)/(K-1)
cat(delta_dot, delta_dot1, delta_dot2, "\n")

eigenvalues <- eigen(solve(P0) %*% V)$values
cat(eigenvalues, "\n")


## Simulation of a_{\alpha} using Eubank's method in his 1992 paper,
## which works much better than the method in his 1997 paper. This
## simulation is trying to get a_alpha for the Weighted data.
alpha <- c(0.01, 0.05, 0.10, 0.20, 0.29) #alpha level
a_alpha <- seq(0, 10, by=0.01)            #create possible a_{alpha} for search
test <- matrix(NA, length(a_alpha), K-1)   #matrix that store each term for each a_{alpha}
for (i in 1:length(a_alpha)){
    for (k in 2:K){

        ## temp is the vector to store values of sum of eigenvalues weighted chi-square 1
        temp <- rep(NA, length=samplenum)
        for (uu in 1:length(temp)){
            chi1 <- rchisq(k-1, df=1)
            temp[uu] <- sum(eigenvalues[1:k-1]*chi1)
        }

        ## numerator_prob is the prob that in the numerator
        numerator_prob <- sum(1*(temp > (k-1)*a_alpha[i]))/length(temp)
        test[i, k-1] <- numerator_prob/(k-1) #probability of each term
    }
    cat("Current a_alpha =", a_alpha[i], "\n")
}

## stat find all differences of possible a_alpha's
stat <- matrix(NA, length(a_alpha), length(alpha))
for (j in 1:length(alpha)){
    stat[,j] <- abs(exp(-apply(test, 1, sum))-(1-alpha[j])) #fomula in Eubank's 1992 paper
}
index <- apply(stat, 2, function(x) which.min(x))#find index which has minimum error
cat("alpha=", alpha, "\n", "a_alpha=", a_alpha[index], "\n")
#################################################################################





#################################################################################
## 
## Simulation of a_{\alpha} using method in Eubank's 1997 paper, which
## doesn't work good, since a_{0.05} is approx. 4.5, not 4.18
## 
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
## hist(a <- 0.05)                            #plot all a_{0.05}
#################################################################################





#################################################################################
##
## Search a_alpha with different value of rho to control alpha=0.05 by simulation
##
numpsu <- 100                            #number of clusters
psusize <- 30                           #number of ssus in each cluster
K <- 5
beta <- seq(0, 0.14, by=0.01)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        ## set real prob for simulation alpha
        allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K #Alternative (21)
    }
}
p <- allp[1,]                       #test for alpha, only need prob under H0
n <- numpsu*psusize                     #total number of units in the popn
T <- c(10000, 100000, 10000)            #for "total loop", "W", "\hat{V}" respectively
alpha <- 0.05
rho <- seq(0.1, 0.9, by=0.1)
## possible values of a_alpha, and these values are manually set.
a_alpha <- rbind(seq(4.35, 4.46, length=10), seq(4.93, 5.16, length=10),
             seq(5.97, 6.13, length=10), seq(7.63, 7.93, length=10),
             seq(10.08, 10.39, length=10), seq(13.63, 14.16, length=10),
             seq(19, 20, length=10), seq(25.1, 26.43, length=10),
             seq(33.15, 34.74, length=10))
## Pr(type I error) for corresponding value of a_alpha 
power_q <- matrix(NA, nrow(a_alpha), ncol(a_alpha))
for (kk in 1:nrow(a_alpha)){
    cat ("rho =", rho[kk], "\n")
    for (jj in 1:ncol(a_alpha)){
        cat("a_alpha=", a_alpha[kk, jj], "\n")
        sim <- eubank(K, numpsu, psusize, rho[kk], p, T, alpha, a_alpha[kk, jj])
        power_q[kk, jj] <- sim$reject_q
    }
}
print(power_q)
#################################################################################






## eubank_run.R
#################################################################################
##
## Test for 1st and 2nd order correction and classicial only.
##
numpsu <- 50                            #number of clusters
psusize <- 10                           #number of ssus in each cluster
rho <- 0.3
K <- 5
beta <- seq(0, 0.14, by=0.01)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        ## set real prob for simulation alpha
        allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K #Alternative (21)
    }
}
print(allp)
print(apply(allp, 1, sum))
n <- numpsu*psusize                 #total number of units in the popn
T <- c(10000, 100000, 10000)        #for "total loop", "W", "\hat{V}" respectively
alpha <- 0.05; a_alpha <- 4.18      #Actually we don't need "a_alpha" in here
power_class <- rep(NA, length(beta))
power_1st <- rep(NA, length(beta))
power_2nd <- rep(NA, length(beta))
for (kk in 1:length(beta)){
    cat ("beta =", beta[kk], "\n")        
    p <- allp[kk,]
    sim <- rao(K, numpsu, psusize, rho, p, T, alpha)
    power_class[kk] <- sim$reject_classical
    power_1st[kk] <- sim$reject_1st
    power_2nd[kk] <- sim$reject_2nd
}

plot(beta, power_class, axes=FALSE, ylim=c(0,1), ylab="Power",
     xlim=c(0,0.14), col=1, lty=1, type="l", lwd=2,
     main=paste("Psu=", numpsu, "," , "Ssu=", psusize, ",", "K=", K, ",", "rho=", rho))
axis(side=1, at=seq(0, 0.14, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_1st, col=2, lty=2, type="l", lwd=1)
points(beta, power_2nd, col=3, lty=3, type="l", lwd=3)
legend("bottomright", c("Classical", "1st", "2nd"),
       col=1:3, lty=1:3, lwd=c(2, 1, 3))
box()
#################################################################################







#################################################################################
##
## Test for Eubank's methods, as well as 1st and 2nd order of the chi-square test
##
alpha <- 0.05
a_alpha <- 6.19
numpsu <- 50                            #number of clusters
psusize <- 10                           #number of ssus in each cluster
rho <- 0.3
K <- 5
beta <- seq(0, 0.14, by=0.01)
allp <- matrix(NA, length(beta), K)
for (i in 1:length(beta)){
    for (k in 1:K){
        ## set real prob for simulation alpha
        allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K #Alternative (21)
    }
}
print(allp)
print(apply(allp, 1, sum))
n <- numpsu*psusize                     #total number of units in the popn
T <- c(10000, 100000, 10000)            #for total loop, W, \hat{V} respectively
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
     main=paste("a_(0.05)=", a_alpha, ",", "n=", n, "," , "K=", K, ",", "rho=", rho))
axis(side=1, at=seq(0, 0.14, by=0.01))
axis(side=2, at=seq(0, 1, by=0.05))
points(beta, power_class, col=2, lty=2, type="l", lwd=1)
points(beta, power_W, col=3, lty=3, type="l", lwd=3)
points(beta, power_1st, col=4, lty=4, type="l", lwd=1.5)
points(beta, power_2nd, col=5, lty=5, type="l", lwd=2.5)
legend("bottomright", c("New Method using q", "Classical", "W", "1st", "2nd"),
       col=1:5, lty=1:5, lwd=c(2, 1, 3, 1.5, 2.5))
box()
#################################################################################



#################################################################################
##
## Test Weighted data for K=3:10, with 1st, 2nd and Eubank's method
## Need to modify "eubank_current.R" to generate Weighted data.
##
## Run the function and generating corresponding pdf plots using loop
for (K in 3:10){                        #number of categories in multinomial

    pdf(paste("K=",K, ".pdf", sep=""), height=20, width=15)
    par(mfrow=c(4, 2))                 #10 plots for 10 values of rho

    beta <- seq(0, 0.14, by=0.01)
    allp <- matrix(NA, length(beta), K)
    for (i in 1:length(beta)){
        for (k in 1:K){
            ## set real prob for simulation alpha
            allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K #Alternative (21)
        }
    }
    print(allp); print(apply(allp, 1, sum))    #check if the sum of probabilities is 1
    for (rho in c(0.01, seq(0.1, 0.9, by=0.1))){
        ## Simulate the power of the new data-driven method, alpha=0.05
        numpsu <- 15                            #number of clusters
        psusize <- 5                            #number of ssus in each cluster
        n <- numpsu*psusize                     #equivalent to n=15*5=75
        T <- c(10000, 100000, 10000)            #for total loop, W, \hat{V} respectively
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
##
## Test Unweighted data for K=3:10, with 1st, 2nd and Eubank's method
## Need to modify "eubank_current.R" to generate Unweighted data.
##
## Run the function and generating corresponding pdf plots using loop
pdf(paste("unweighted", ".pdf", sep=""), height=20, width=15)
par(mfrow=c(4, 2))                 #10 plots for 10 values of rho
for (K in 3:10){                        #number of categories in multinomial
    beta <- seq(0, 0.14, by=0.01)
    allp <- matrix(NA, length(beta), K)
    for (i in 1:length(beta)){
        for (k in 1:K){
            ## set real prob for simulation alpha
            allp[i, k] <- 1/K + beta[i]*(k-median(1:K))/K #Alternative (21)
        }
    }
    print(allp); print(apply(allp, 1, sum))    #check if the sum of probabilities is 1
    rho <- 0.001       #we don't need rho for unweighted, just an input for my function
    ## Simulate the power of the new data-driven method, alpha=0.05
    numpsu <- 15                            #number of clusters
    psusize <- 5                            #number of ssus in each cluster
    n <- numpsu*psusize                     #equivalent to n=15*5=75
    ## rho <- 0.3                              #ICC in cluster
    ## n <- 75                                 #75 draws for multinomial distribution
    ## n <- 150                                #150 draws for multinomial distribution
    T <- c(10000, 100000, 10000)        #for total loop, W, \hat{V} respectively
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
#################################################################################



#################################################################################
## The following code simulates when alpha=0.05 for Unweighted data
## Need to modify "eubank_current.R" to generate Unweighted data.
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
T <- c(10000, 100000, 10000)            #for total loop, W, \hat{V} respectively
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
