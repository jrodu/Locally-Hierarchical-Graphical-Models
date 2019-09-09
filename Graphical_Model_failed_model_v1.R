


#  Graphical Model
library(glmnet)
#install.packages("rmutil")
library(rmutil)
library(stats4)
#  read data: TCGA endometrial cancer gene mutation data with 59 genes
df <- read.csv("C:\\Users\\mxy00\\Desktop\\endometrial_cancer.csv", header = T, 
               sep = "\t")
df <- df[, -1]

#  for loop for regressing each gene with the rest genes by GLMNET
#  with L1 regularization
lambda <- 10
epsilon <- 10000
beta_mat <- matrix(0, nrow = ncol(df), ncol = ncol(df))
#lambda_L <- matrix(0, nrow = ncol(df), ncol = 1)
for (i in 1:ncol(df)){
  aa <- 1
  while (epsilon > 0.01) {
    #print(length(epsilon))
    fit <- glmnet(x = data.matrix(df[, -i]), y = data.matrix(df[, i]), 
                  family = "poisson", alpha = 0, 
                  #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                  lambda = lambda)
    #lambda_L[i, 1] <- lambda
    beta_mat[-i, i] <- as.data.frame(as.matrix(fit$beta))[, 1]
    sigma_squared <- sum(beta_mat[-i, i]^2)/(nrow(beta_mat) - 1)
    penalty <- 1/(2*sigma_squared)
    epsilon <- abs(lambda - penalty)/lambda
    print(c(aa, epsilon, lambda, penalty, sigma_squared))
    lambda <- penalty
    aa <- aa + 1
  }
}    
#  penalty goes to infinite

    
    ###################################
    # #  assume beta, i.e. coefficients follow Laplacian Distribution
    # #  write MLE function for the prior
    # #  formulate the log-likelihood function
    # y <- fit$beta
    # mu <- mean(y)
    # sigma <- sd(y)
    # LL <- function(mu, sigma){
    #   R = dlaplace(y, m = mu, s = sigma) - sum(log(R))
    # }
    # #  Apply MLE to estimate the parameters
    # MLE <- mle(LL(mu, sigma), start = list(mu = 1, sigma = 1)) #  problem: sd cannot be zero
    # #  method1: compare two loglikelihood function
    # #logL <- summary(MLE)$-2 log L
    # #  method2: use estimated parameters to implement the penalty on glm
    # estimated_mu <- MLE$coefficient[1]
    # estimated_sig <- MLE$coefficient[2]
    #####################################
    #  assume beta, i.e. coefficients follow Gaussian Distribution
    # y <- fit$beta
    # sigma <- sd(y)
    # LL <- function(mu, sigma){
    #   mu <- 0
    #   loglikelihood <- - sum(log(dnorm(as.vector(as.matrix(y)), mu, sigma)))
    #   print(loglikelihood)
    #   }
    # #  Apply MLE to estimate the parameters
    # MLE <- mle(LL(mu, sigma), start = list(lambda = 1)) 
    # #  method1: compare two loglikelihood function
    # #logL <- summary(MLE)$-2 log L
    # #  method2: use estimated parameters to implement the penalty on glm
    # estimated_mu <- 0  #  find sigma
    # estimated_sig <- MLE$coefficient[2]

