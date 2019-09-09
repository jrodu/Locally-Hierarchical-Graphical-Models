


#  Graphical Model
library(glmnet)
#install.packages("rmutil")
library(rmutil)
library(stats4)
#  read data: TCGA endometrial cancer gene mutation data with 59 genes
df <- read.csv("endometrial_cancer.csv", header = T,sep = "\t")
df <- df[, -1]

for(ii in 1:ncol(df)){
  df[, ii] <- as.integer(df[, ii] != 0)
}

glm_family <- "binomial"



best_model_given_prior <- function(X, y, lambda.vec, prior.a, prior.b) {
  X <- data.matrix(X)
  log_lik.vec <- c()
  fit.lasso <- glmnet(x = X, y = y, 
                  family = glm_family, alpha = 1, 
                  #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                  lambda = lambda.vec)
  
  for(ii in c(1:ncol(fit.lasso$beta))){
    is_in_model <- as.vector(fit.lasso$beta[, ii] != 0)
    p <- mean(is_in_model)
    if(sum(is_in_model) == 0){
      fit <- glm(y ~ 1, family = glm_family)
    }
    else{
      tmp_X <- X[, is_in_model]
      fit <- glm(y ~ tmp_X, family = glm_family)
    }
    log_prior <- dbeta(p, prior.a, prior.b, log = T)
    log_lik <- logLik(fit) + log_prior
    log_lik.vec <- c(log_lik.vec, log_lik)
    print(c(lambda.vec[ii], log_prior, logLik(fit), p))
  }
  if(ncol(fit.lasso$beta) != length(lambda.vec)){
    print("not converged")
  }
    
  lambda_final <- lambda.vec[order(log_lik.vec, decreasing = T)[1]]
  fit <- glmnet(x = X, y = y, 
                family = glm_family, alpha = 1, 
                #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                lambda = lambda_final)
  if(sum(as.vector(fit$beta != 0)) == 0){
    fit <- glm(y ~ 1, family = glm_family)
  }
  else{
    tmp_X <- X[, as.vector(fit$beta != 0)] 
    fit <- glm(y ~ tmp_X, family = glm_family)
  }
  p <- mean(fit.lasso$beta[, order(log_lik.vec, decreasing = T)[1]] != 0)
  return(list(fit = fit, lambda_final = lambda_final, p.connection = p))
}




lambda_vec <- c(10, 1, 0.1, 0.001, 0.0001)
beta_mat <- matrix(0, nrow = ncol(df), ncol = ncol(df))
for (i in 1:ncol(df)){

    #print(length(epsilon))
    res <- best_model_given_prior(
           X = data.matrix(df[, -i]), y = data.matrix(df[, i]),
            lambda.vec = lambda_vec, prior.a = 2, prior.b = 2)
    fit <- res$fit
    lambda_final <- res$lambda_final
    p.connection <- res$p.connection
    #beta_mat[-i, i] <- fit$coefficients
    print(c(lambda_final, p.connection))
}  



