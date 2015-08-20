## Deleted 2 lines and added this new line

#####################################################
#
#		Title: Homework 1: Introduction to R
#		Created: 2011/08/25
#		Due Date: 2011/09/06
#  	Description: The purpose of this file is the introduce 
#			     you to the R programing language and to 
#			     show you how to run simple Simulations in R. 
#
#####################################################

#Start by declaring a working directory.  This is the directory where you 
# are storing your files (ie, hw1.R). This is also where R will look for files.
setwd('/homes/sonksen/courses/470/')

#note, in Windows, the will probably be something like setwd('c://courses/470/')


# Variables in R do not need to be declared. ``<-'' is used to assign 
# a value to a variable. ``='' may also be used. 

n <- 5  #n will denote our sample size here

# next, we will randomly simulate from Model (1) in the homework.
# Here, we will assume that the Cov(x, epsilon) is zero.  
# The rnorm function simulates normal variables. The first argument 
# is the number of random variables to be simulated. The second argument
# gives the mean and the third argument the standard deviation. 

x <- rnorm(n,0,1)
epsilon <- rnorm(n,3,4)

# running ?functionname will give the help file for functionname. 
?rnorm

y <- x+epsilon

# note that y is a vector. Like Matlab, R is ``vectorized''
print(y)

# most standard statitics (mean, median, variance) have functions in R. 
# The var function calculates the variance of a vector. 

var(y)

# Based on y and epsilon, an estimate of sigma_x is
# note that we used the sqrt, var, and max functions. 
  
hat_sigma_x <- sqrt(max(0,var(y)-var(epsilon)))
print(hat_sigma_x)

#recall, the true value of sigma_x is 1



# A concern of hat_sigma_x is that it can be zero. 
# To see how often this happens, let us repeat this 
# simulation 1000 times. 


#T will denote the number of iterations
T<-1000

# the rep function creates a vector. The first argument is the value 
# of each item in the vector and the second argument is the length of the 
# vector. 

hat_sigma_x<-rep(0,T)

# for loops are run using the for function. The commands to be 
# repeated at each iteration of the loop are contained in curly bracks
 
for(t in 1:T){

# this reads, for t in (1, 2, ..., T) do the following. 

  x<-rnorm(n,0,1)
  epsilon<-rnorm(n,3,4)
  y<-x+epsilon

#square braces allow us to assign a value to an element of a vector.

   hat_sigma_x[t] <-sqrt(max(0,var(y)-var(epsilon)))
}# end loop

#finally, lets see what percent of our estimators are equal to zero. 

# Logical statements in R return a TRUE or FALSE. Note ``=='' denotes equality

hat_sigma_x == 0

# To see how many of these are equal to zero, try the sum command. 
# The sum command converts ``TRUE'' to 1 and ``FALSE'' to 0. 

sum(hat_sigma_x ==0)

# Finally, divide this by T to get the percent which are zero. 

sum(hat_sigma_x ==0)/T

# Note, this simulation is based on random draws from distributions. Thus, 
# everyone's answers can be different.  Also, this study is heavily 
# dependent on the assumed values of the sample size and the standard deviations 
# of x and  epsilon. If we change these values, the estimator may work better or 
# worse.  
