

instability_measure <- function(L_adj_mat.array) {
  mean_adj_mat <- apply(L_adj_mat.array, c(1, 2), mean)
  p <- nrow(mean_adj_mat)
  tot_instability <- sum(2 * mean_adj_mat * (1 - mean_adj_mat)) / (p * (p - 1))
  return(tot_instability)
}



allen_method <- function(df, seed = 434, n_lambda = 100, B = 30){
  
  lambda_min <- 10 ^ (-4)
  
  tmp_mat <- data.matrix(df)
  tmp_mat <- t(tmp_mat) %*% tmp_mat
  diag(tmp_mat) <- 0
  lambda_max <- max(tmp_mat)
  lambda.vec <- exp(seq(from = log(lambda_min), to = log(lambda_max), length.out = n_lambda))
  
  glm_family <- "binomial"
  
  B <- B
  
  n <- nrow(df)
  n_subsample <- floor(10 * sqrt(n))
  p <- ncol(df)
  
  set.seed(seed)
  
  instability.vec <- c()
  n_edge.vec <- c()
  
  for(lambda in lambda.vec){
    L_adj_mat.array <- array(0, dim = c(p, p, B))
    for(b in c(1:B)){
      df_sample <- df[sample.int(n, n_subsample),]
      beta_mat <- matrix(0, nrow = p, ncol = p)
      for (i in 1:ncol(df)){
        fit <- glmnet(x = data.matrix(df_sample[, -i]), y = data.matrix(df_sample[, i]), 
                      family = glm_family, alpha = 1, 
                      #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                      lambda = lambda)
        beta_mat[-i, i] <- as.data.frame(as.matrix(fit$beta))[, 1]
      }
      L_adj_mat.array[, , b] <- as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)
    }
    instability <- instability_measure(L_adj_mat.array = L_adj_mat.array)
    instability.vec <- c(instability.vec, instability)
    
    beta_mat <- matrix(0, nrow = p, ncol = p)
    for (i in 1:ncol(df)){
      fit <- glmnet(x = data.matrix(df[, -i]), y = data.matrix(df[, i]), 
                    family = glm_family, alpha = 1, 
                    #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                    lambda = lambda)
      beta_mat[-i, i] <- as.data.frame(as.matrix(fit$beta))[, 1]
    }
    n_edge.vec <- c(n_edge.vec, sum(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)) / 2)
    
    print(c(lambda, instability))
    
  }
  plot(log(lambda.vec), instability.vec)
  lambda_candidate.vec <- c()
  for(i in c(1:length(instability.vec))){
    if(instability.vec[i] <= 0.05){
      lambda_candidate.vec <- c(lambda_candidate.vec, lambda.vec[i])
    }
  }
  lambda_opt <- min(lambda_candidate.vec)
  print(lambda_opt)
  
  beta_mat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:ncol(df)){
    fit <- glmnet(x = data.matrix(df[, -i]), y = data.matrix(df[, i]), 
                  family = glm_family, alpha = 1, 
                  #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                  lambda = lambda_opt)
    beta_mat[-i, i] <- as.data.frame(as.matrix(fit$beta))[, 1]
  }
  adj_mat <- as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)
  df_summary <- data.frame(lambda = lambda.vec, instability = instability.vec,
                           n_edge = n_edge.vec)
  
  res <- list(adj_mat = adj_mat, lambda_opt = lambda_opt, summary = df_summary)
  return(res)
}


best_model_given_prior <- function(X, y, lambda.vec, prior.a, prior.b, glm_family, prior_weight = 10) {
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
    log_lik <- logLik(fit) + log_prior * prior_weight
    log_lik.vec <- c(log_lik.vec, log_lik)
    #print(c(lambda.vec[ii], log_prior, logLik(fit), p))
  }
  if(ncol(fit.lasso$beta) != length(lambda.vec)){
    #print("not converged")
  }
  
  lambda_final <- lambda.vec[order(log_lik.vec, decreasing = T)[1]]
  fit <- glmnet(x = X, y = y, 
                family = glm_family, alpha = 1, 
                #lambda = c(10^(-3), 10^(-2), 10^(-1), 1, 10^(1), 10^(2), 10^(3)))
                lambda = lambda_final)
  is_in_model <- as.integer(as.vector(fit$beta != 0))
  if(sum(is_in_model) == 0){
    fit <- glm(y ~ 1, family = glm_family)
  }
  else{
    fit <- glm(y ~ X[, is_in_model] , family = glm_family)
  }
  p <- mean(fit.lasso$beta[, order(log_lik.vec, decreasing = T)[1]] != 0)
  return(list(fit = fit, lambda_final = lambda_final, p.connection = p, edges = is_in_model))
}


my_prior_stability_search <- function(df, prior_a.vec, prior_b.vec, seed = 434, n_lambda = 100, B = 30){
  
  lambda_min <- 10 ^ (-4)
  
  tmp_mat <- data.matrix(df)
  tmp_mat <- t(tmp_mat) %*% tmp_mat
  diag(tmp_mat) <- 0
  lambda_max <- max(tmp_mat)
  lambda.vec <- exp(seq(from = log(lambda_min), to = log(lambda_max), length.out = n_lambda))
  lambda.vec <- sort(lambda.vec, decreasing = T)
  
  glm_family <- "binomial"
  
  B <- B
  
  n <- nrow(df)
  n_subsample <- floor(10 * sqrt(n))
  p <- ncol(df)
  
  set.seed(seed)
  
  instability.vec <- c()
  n_edge.vec <- c()
  
  for(ii in c(1:length(prior_a.vec))){
    prior.a <- prior_a.vec[ii]
    prior.b <- prior_b.vec[ii]
    L_adj_mat.array <- array(0, dim = c(p, p, B))
    for(b in c(1:B)){
      df_sample <- df[sample.int(n, n_subsample),]
      beta_mat <- matrix(0, nrow = p, ncol = p)
      for (i in 1:ncol(df)){
        fit <- best_model_given_prior(X = data.matrix(df_sample[, -i]), y = data.matrix(df_sample[, i]),
                               lambda.vec = lambda.vec, glm_family = glm_family,
                               prior.a = prior.a, prior.b = prior.b)
        beta_mat[-i, i] <- fit$edges
      }
      L_adj_mat.array[, , b] <- as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)
    }
    instability <- instability_measure(L_adj_mat.array = L_adj_mat.array)
    instability.vec <- c(instability.vec, instability)
    
    beta_mat <- matrix(0, nrow = p, ncol = p)
    lambda_final.vec <- c()
    for (i in 1:ncol(df)){
      fit <- best_model_given_prior(X = data.matrix(df[, -i]), y = data.matrix(df[, i]),
                                    lambda.vec = lambda.vec, glm_family = glm_family,
                                    prior.a = prior.a, prior.b = prior.b)
      lambda_final <- fit$lambda_final
      lambda_final.vec <- c(lambda_final.vec, lambda_final)
      beta_mat[-i, i] <- fit$edges
    }
    print(lambda_final.vec)
    n_edge <- sum(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)) / 2
    n_edge.vec <- c(n_edge.vec, n_edge)
    print(c(prior.a, prior.b, n_edge, instability))
  }
  
  df_summary <- data.frame(prior_a = prior_a.vec, prior_b = prior_b.vec, 
                           instability = instability.vec,
                           n_edge = n_edge.vec)
  res <- list(summary = df_summary)
  return(res)
}

# #  Graphical Model
# library(glmnet)
# #install.packages("rmutil")
# library(rmutil)
# library(stats4)
# #  read data: TCGA endometrial cancer gene mutation data with 59 genes
# df <- read.csv("endometrial_cancer.csv", header = T,sep = "\t")
# df <- df[, -1]
# 
# for(ii in 1:ncol(df)){
#   df[, ii] <- as.integer(df[, ii] != 0)
# }
