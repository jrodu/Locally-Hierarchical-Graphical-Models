library(glmnet)
library(glm2)
library(mpath)
#install.packages("rmutil")
library(rmutil)
library(stats4)
library(MASS)
library(ggplot2)
library(pROC) # only for ROC curves and AUC
library(PoisNor) # only for generating multivariate Poisson
library(igraph) # only for generating the adj mat of scale-free network
library(xtable) # only for exporting latex table
library(XMRF) # LPGM and Data generation


#---- Functions ----

instability_measure <- function(L_adj_mat.array) {
  mean_adj_mat <- apply(L_adj_mat.array, c(1, 2), mean)
  p <- nrow(mean_adj_mat)
  tot_instability <- sum(mean_adj_mat * (1 - mean_adj_mat)) / (p * (p - 1))
  return(tot_instability)
}


allen_method <- function(df, seed = 434, B = 30, lambda.vec = NULL, n_lambda = 100,
                         lambda_min = 10^-4, lambda_max = NULL, glm_family = "binomial",
                         is_print = T, is_sparse_graph = T){

  if(is.null(lambda.vec)){
    if(is.null(lambda_max)){
      tmp_mat <- data.matrix(df)
      tmp_mat <- t(tmp_mat) %*% tmp_mat
      diag(tmp_mat) <- 0
      lambda_max <- max(tmp_mat)
    }
    lambda.vec <- exp(seq(from = log(lambda_max), to = log(lambda_min), length.out = n_lambda))
  }

  B <- B

  n <- nrow(df)
  n_subsample <- floor(10 * sqrt(n))
  p <- ncol(df)

  set.seed(seed)

  fit_subsample.list <- list()
  for(b in c(1:B)){
    cat("subsample", b, "starting", "\n")
    flag <- T
    while(flag){
      flag <- F
      df_sample <- df[sample.int(n, n_subsample),]
      fit_subsample.list[[b]] <- list()
      for (i in 1:ncol(df)){
        fit_subsample.list[[b]][[i]] <- glmnet(
          x = data.matrix(df_sample[, -i]), y = data.matrix(df_sample[, i]),
          family = glm_family, alpha = 1,  lambda = lambda.vec)

        if(ncol(fit_subsample.list[[b]][[i]]$beta) != length(lambda.vec)){
          flag <- T
          break
        }
      }
    }
  }

  fit_df.list <- list()
  for(i in 1:ncol(df)){
    fit_df.list[[i]] <- glmnet(x = data.matrix(df[, -i]), y = data.matrix(df[, i]),
                               family = glm_family, alpha = 1,
                               lambda = lambda.vec)
  }

  instability.vec <- c()
  n_edge.vec <- c()
  adj_mat.list <- list()
  for(jj in 1:length(lambda.vec)){
    lambda <- lambda.vec[jj]
    L_adj_mat.array <- array(0, dim = c(p, p, B))
    for(b in c(1:B)){
      df_sample <- df[sample.int(n, n_subsample),]
      beta_mat <- matrix(0, nrow = p, ncol = p)
      for (i in 1:ncol(df)){
        fit <- fit_subsample.list[[b]][[i]]
        beta_mat[-i, i] <- fit$beta[, jj]
      }
      L_adj_mat.array[, , b] <- as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)
    }
    instability <- instability_measure(L_adj_mat.array = L_adj_mat.array)
    instability.vec <- c(instability.vec, instability)

    beta_mat <- matrix(0, nrow = p, ncol = p)
    for (i in 1:ncol(df)){
      beta_mat[-i, i] <- fit_df.list[[i]]$beta[, jj]
    }
    n_edge <- sum(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)) / 2
    adj_mat <- matrix(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0), nrow = nrow(beta_mat))
    n_edge.vec <- c(n_edge.vec, n_edge)
    adj_mat.list[[jj]] <- adj_mat

    if(is_print){
      print(c(lambda, instability, n_edge))
    }
  }
  #plot(log(lambda.vec), instability.vec)
  lambda_candidate.vec <- c()
  for(i in c(1:length(instability.vec))){
    if(instability.vec[i] <= 0.05){
      if(is_sparse_graph){
        if(n_edge.vec[i] < p * (p - 1) / 4){
          lambda_candidate.vec <- c(lambda_candidate.vec, lambda.vec[i])
        }
      }
      else{
        lambda_candidate.vec <- c(lambda_candidate.vec, lambda.vec[i])
      }

    }
  }
  lambda_opt <- min(lambda_candidate.vec)
  index_opt <- c(1:length(lambda.vec))[lambda.vec == lambda_opt]
  if(is_print){
    print(lambda_opt)
  }

  beta_mat <- matrix(0, nrow = p, ncol = p)
  for (i in 1:ncol(df)){
    beta_mat[-i, i] <- fit_df.list[[i]]$beta[, index_opt]
  }
  adj_mat <- matrix(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0), nrow = nrow(beta_mat))
  df_summary <- data.frame(lambda = lambda.vec, instability = instability.vec,
                           n_edge = n_edge.vec)

  res <- list(adj_mat = adj_mat, lambda_opt = lambda_opt, summary = df_summary, adj_mat.list = adj_mat.list)
  return(res)
}


my_method_given_prior <- function(df, prior.a = 1, prior.b = 1, glm_family = "binomial",
                                  lambda.vec = NULL, n_lambda = 100, lambda_min = 10^-4, lambda_max = NULL,
                                  prior_weight = NULL){
  if(is.null(lambda.vec)){
    if(is.null(lambda_max)){
      tmp_mat <- data.matrix(df)
      tmp_mat <- t(tmp_mat) %*% tmp_mat
      diag(tmp_mat) <- 0
      lambda_max <- max(tmp_mat)
    }
    lambda.vec <- exp(seq(from = log(lambda_max), to = log(lambda_min), length.out = n_lambda))
  }

  if(is.null(prior_weight)){
    prior_weight <- log(nrow(df))
  }

  p <- ncol(df)

  beta_mat <- matrix(0, nrow = p, ncol = p)
  for(ii in 1:p){
    fit <- best_model_given_prior(df[, -ii], df[, ii], prior.a = prior.a, prior.b = prior.b,
                                  glm_family = glm_family, theta = theta, lambda.vec = lambda.vec,
                                  prior_weight = prior_weight)
    beta_mat[-ii, ii] <- fit$edges
  }
  adj_mat <- matrix(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0), nrow = nrow(beta_mat))
  return(list(adj_mat = adj_mat))
}

best_model_given_prior <- function(X, y, prior.a = 1, prior.b = 1, glm_family = "binomial", theta = NULL,
                                   lambda.vec = NULL, n_lambda = 100, lambda_min = 10^-4, lambda_max = NULL,
                                   prior_weight = NULL) {

  if(is.null(lambda.vec)){
    if(is.null(lambda_max)){
      tmp_mat <- data.matrix(df)
      tmp_mat <- t(tmp_mat) %*% tmp_mat
      diag(tmp_mat) <- 0
      lambda_max <- max(tmp_mat)
    }
    lambda.vec <- exp(seq(from = log(lambda_max), to = log(lambda_min), length.out = n_lambda))
  }

  X <- data.matrix(X)
  if(glm_family != "negbin"){
    fit.lasso <- glmnet(x = X, y = y,
                        family = glm_family, alpha = 1,
                        lambda = lambda.vec)
  }
  else{
    fit.lasso <- glmregNB(y ~ X, alpha = 1,
                        lambda = lambda.vec, maxit.theta = 50)
  }

  if(is.null(prior_weight)){
    prior_weight <- log(nrow(X))
  }

  log_lik.vec <- c()
  log_prior.vec <- c()
  log_samplelik.vec <- c()
  for(ii in c(1:ncol(fit.lasso$beta))){
    is_in_model <- as.vector(fit.lasso$beta[, ii] != 0)
    p <- mean(is_in_model)
    if(sum(is_in_model) == 0){
      if(glm_family != "negbin"){
        fit <- glm(y ~ 1, family = glm_family)
        log_samplelik <- logLik(fit)
      }
      else{
        fit <- glm.nb(y ~ 1)
        log_samplelik <- logLik(fit)
      }
    }
    else{
      if(glm_family != "negbin"){
        fit <- glm2(y ~ X[, is_in_model], family = glm_family)
        log_samplelik <- logLik(fit)
      }
      else{
        #beta_start <- fit.lasso$beta[, ii]
        #fit <- glm.nb(y ~ X[, is_in_model], method = "glm.fit2",start = c(fit.lasso$b0[ii], beta_start[beta_start != 0]), init.theta = fit.lasso$theta[ii])
        log_samplelik <- fit.lasso$resdev[ii]
      }
    }
    log_prior <- dbeta(p, prior.a, prior.b, log = T)
    log_lik <- log_samplelik + log_prior * prior_weight
    log_lik.vec <- c(log_lik.vec, log_lik)
    log_prior.vec <- c(log_prior.vec, log_prior)
    log_samplelik.vec <- c(log_samplelik.vec, log_samplelik)
  }
  if(ncol(fit.lasso$beta) != length(lambda.vec)){
    #print("not converged")
  }

  index_opt <- order(log_lik.vec, decreasing = T)[1]
  lambda_final <- lambda.vec[index_opt]
  is_in_model <- as.vector(fit.lasso$beta[, index_opt] != 0)
  if(sum(is_in_model) == 0){
    if(glm_family != "negbin"){
      fit <- glm(y ~ 1, family = glm_family)
    }
    else{
      fit <- glm.nb(y ~ 1)
    }
  }
  else{
    if(glm_family != "negbin"){
      fit <- glm2(y ~ X[, is_in_model], family = glm_family)
    }
    else{
      #beta_start <- fit.lasso$beta[, index_opt]
      #fit <- glm.nb(y ~ X[, is_in_model], method = "glm.fit2", start = c(fit.lasso$b0[index_opt], beta_start[beta_start != 0]), init.theta = fit.lasso$theta[index_opt])
      fit <- fit.lasso
    }
  }
  p <- mean(is_in_model)
  return(list(fit = fit, lambda_final = lambda_final, p.connection = p, edges = is_in_model,
              log_lik.vec = log_lik.vec,
              log_prior.vec = log_prior.vec, log_samplelik.vec = log_samplelik.vec,
              edges.vec = fit.lasso$beta != 0))
}

best_model_given_prior_brute <- function(X, y, prior.a = 1, prior.b = 1, glm_family = "binomial", theta = NULL,
                                   prior_weight = NULL) {


  X <- data.matrix(X)

  if(is.null(prior_weight)){
    prior_weight <- log(nrow(X))
  }

  tmp_args <- list()
  for(ii in 1:ncol(X)){
    tmp_args[[ii]] <- c(F, T)
  }
  is_in_model.df <- do.call(expand.grid, tmp_args)

  log_lik.vec <- c()
  log_prior.vec <- c()
  log_samplelik.vec <- c()
  for(ii in 1:nrow(is_in_model.df)){
    is_in_model <- as.logical(is_in_model.df[ii, ])
    p <- mean(is_in_model)
    if(sum(is_in_model) == 0){
      if(glm_family != "negbin"){
        fit <- glm(y ~ 1, family = glm_family)
        log_samplelik <- logLik(fit)
      }
      else{
        fit <- glm.nb(y ~ 1)
        log_samplelik <- logLik(fit)
      }
    }
    else{
      if(glm_family != "negbin"){
        fit <- glm2(y ~ X[, is_in_model], family = glm_family)
        log_samplelik <- logLik(fit)
      }
      else{
        #beta_start <- fit.lasso$beta[, ii]
        #fit <- glm.nb(y ~ X[, is_in_model], method = "glm.fit2",start = c(fit.lasso$b0[ii], beta_start[beta_start != 0]), init.theta = fit.lasso$theta[ii])
        log_samplelik <- fit.lasso$resdev[ii]
      }
    }
    log_prior <- dbeta(p, prior.a, prior.b, log = T)
    log_lik <- log_samplelik + log_prior * prior_weight
    log_lik.vec <- c(log_lik.vec, log_lik)
    log_prior.vec <- c(log_prior.vec, log_prior)
    log_samplelik.vec <- c(log_samplelik.vec, log_samplelik)
  }
  if(ncol(fit.lasso$beta) != length(lambda.vec)){
    #print("not converged")
  }

  index_opt <- order(log_lik.vec, decreasing = T)[1]
  is_in_model <- as.logical(is_in_model.df[index_opt, ])
  if(sum(is_in_model) == 0){
    if(glm_family != "negbin"){
      fit <- glm(y ~ 1, family = glm_family)
    }
    else{
      fit <- glm.nb(y ~ 1)
    }
  }
  else{
    if(glm_family != "negbin"){
      fit <- glm2(y ~ X[, is_in_model], family = glm_family)
    }
    else{
      #beta_start <- fit.lasso$beta[, index_opt]
      #fit <- glm.nb(y ~ X[, is_in_model], method = "glm.fit2", start = c(fit.lasso$b0[index_opt], beta_start[beta_start != 0]), init.theta = fit.lasso$theta[index_opt])
      fit <- fit.lasso
    }
  }
  p <- mean(is_in_model)
  return(list(fit = fit, p.connection = p, edges = is_in_model,
              log_lik.vec = log_lik.vec,
              log_prior.vec = log_prior.vec, log_samplelik.vec = log_samplelik.vec))
}


bruteforce_model_given_prior <- function(X, y, prior.a, prior.b, glm_family, prior_weight = 5) {
  X <- data.matrix(X)
  n <- ncol(X)

  bin.list <- list()
  for(ii in 1:n){
    bin.list[[ii]] <- c(F, T)
  }
  bin.df <- t(expand.grid(bin.list))

  log_lik.vec <- c()
  log_prior.vec <- c()
  log_samplelik.vec <- c()
  for(ii in 1:ncol(bin.df)){
    is_in_model <- bin.df[, ii]
    if(sum(is_in_model) != 0){
      tmp_X <- data.matrix(X[, is_in_model])
      fit <- glm(y ~ tmp_X, family = glm_family)
    }
    else{
      fit <- glm(y ~ 1, family = glm_family)
    }
    p <- mean(is_in_model)
    log_prior <- dbeta(p, prior.a, prior.b, log = T)
    log_lik <- logLik(fit) + log_prior * prior_weight
    log_lik.vec <- c(log_lik.vec, log_lik)
    log_prior.vec <- c(log_prior.vec, log_prior)
    log_samplelik.vec <- c(log_samplelik.vec, logLik(fit))
  }

  is_in_model <- bin.df[, order(log_lik.vec, decreasing = T)[1]]
  if(sum(is_in_model) != 0){
    tmp_X <- data.matrix(X[, is_in_model])
    fit <- glm(y ~ tmp_X, family = glm_family)
  }
  else{
    fit <- glm(y ~ 1, family = glm_family)
  }
  p <- mean(is_in_model)
  return(list(fit = fit, p.connection = p, edges = is_in_model,
              bin.df = bin.df, log_lik.vec = log_lik.vec,
              log_prior.vec = log_prior.vec, log_samplelik.vec = log_samplelik.vec))

}


my_prior_stability_search <- function(df, prior_a.vec, prior_b.vec, seed = 434, n_lambda = 100, B = 30, is_print = T){

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
  adj_mat.list <- list()

  for(ii in c(1:length(prior_a.vec))){
    prior.a <- prior_a.vec[ii]
    prior.b <- prior_b.vec[ii]
    L_adj_mat.array <- array(0, dim = c(p, p, B))
    for(b in c(1:B)){
      df_sample <- df[sample.int(n, n_subsample),]
      beta_mat <- matrix(0, nrow = p, ncol = p)
      for (i in 1:ncol(df)){
        fit <- best_model_given_prior(X = data.matrix(df_sample[, -i]), y = data.matrix(df_sample[, i]),
                               lambda.vec = lambda.vec, glm_family = glm_family, theta = theta,
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
                                    lambda.vec = lambda.vec, glm_family = glm_family, theta = theta,
                                    prior.a = prior.a, prior.b = prior.b)
      lambda_final <- fit$lambda_final
      lambda_final.vec <- c(lambda_final.vec, lambda_final)
      beta_mat[-i, i] <- fit$edges
    }

    adj_mat <- matrix(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0), nrow = nrow(beta_mat))
    n_edge <- sum(adj_mat) / 2
    n_edge.vec <- c(n_edge.vec, n_edge)
    adj_mat.list[[ii]] <- adj_mat
    if(is_print){
      print(lambda_final.vec)
      print(c(prior.a, prior.b, n_edge, instability))
    }
  }

  df_summary <- data.frame(prior_a = prior_a.vec, prior_b = prior_b.vec,
                           instability = instability.vec,
                           n_edge = n_edge.vec)

  tmp <- df_summary
  tmp$index <- c(1: nrow(tmp))
  tmp <- tmp[tmp$instability <= 0.05, ]
  tmp <- tmp[order(tmp$n_edge, decreasing = T), ]
  tmp_index <- tmp[1, "index"]
  adj_mat <- adj_mat.list[[tmp_index]]

  res <- list(summary = df_summary, adj_mat.list = adj_mat.list, adj_mat = adj_mat)
  return(res)
}


my_prior_stability_search_samesubsample <- function(df, prior_a.vec, prior_b.vec, glm_family = "binomial", theta = NULL,
                                                    seed = 434, n_lambda = 100, B = 30, prior_weight = NULL, is_print = T,
                                                    is_sparse_graph = T){

  if(length(prior_a.vec) == 1){
    prior_a.vec <- rep(prior_a.vec, length(prior_b.vec))
  }

  lambda_min <- 10 ^ (-4)

  tmp_mat <- data.matrix(df)
  tmp_mat <- t(tmp_mat) %*% tmp_mat
  diag(tmp_mat) <- 0
  lambda_max <- max(tmp_mat)
  lambda.vec <- exp(seq(from = log(lambda_min), to = log(lambda_max), length.out = n_lambda))
  lambda.vec <- sort(lambda.vec, decreasing = T)

  B <- B

  n <- nrow(df)
  n_subsample <- floor(10 * sqrt(n))
  p <- ncol(df)

  if(is.null(prior_weight)){
    prior_weight <- log(nrow(df))
  }

  set.seed(seed)

  fit_subsample.list <- list()
  for(b in c(1:B)){
    flag <- T
    while(flag){
      flag <- F
      df_sample <- df[sample.int(n, n_subsample),]
      fit_subsample.list[[b]] <- list()
      for (i in 1:ncol(df)){
        fit_subsample.list[[b]][[i]] <- best_model_given_prior(X = data.matrix(df_sample[, -i]), y = data.matrix(df_sample[, i]),
                                                               lambda.vec = lambda.vec, glm_family = glm_family, theta = theta,
                                                               prior.a = prior_a.vec[1], prior.b = prior_b.vec[1])
        if(length(fit_subsample.list[[b]][[i]]$log_lik.vec) != length(lambda.vec)){
          flag <- T
          break
        }
      }
    }
    cat("subsample", b, "finished", "\n")
  }

  instability.vec <- c()
  n_edge.vec <- c()
  adj_mat.list <- list()

  for(ii in c(1:length(prior_a.vec))){
    prior.a <- prior_a.vec[ii]
    prior.b <- prior_b.vec[ii]
    L_adj_mat.array <- array(0, dim = c(p, p, B))
    for(b in c(1:B)){
      beta_mat <- matrix(0, nrow = p, ncol = p)
      for (i in 1:ncol(df)){
        fit <- fit_subsample.list[[b]][[i]]
        log_lik.vec <- fit$log_samplelik.vec + dbeta(apply(fit$edges.vec, 2, mean), prior.a, prior.b, log = T) * prior_weight
        index_opt <- order(log_lik.vec, decreasing = T)[1]
        beta_mat[-i, i] <- fit$edges.vec[, index_opt]
      }
      L_adj_mat.array[, , b] <- as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0)
    }
    instability <- instability_measure(L_adj_mat.array = L_adj_mat.array)
    instability.vec <- c(instability.vec, instability)

    beta_mat <- matrix(0, nrow = p, ncol = p)
    lambda_final.vec <- c()
    for (i in 1:ncol(df)){
      fit <- best_model_given_prior(X = data.matrix(df[, -i]), y = data.matrix(df[, i]),
                                    lambda.vec = lambda.vec, glm_family = glm_family, theta = theta,
                                    prior.a = prior.a, prior.b = prior.b)
      lambda_final <- fit$lambda_final
      lambda_final.vec <- c(lambda_final.vec, lambda_final)
      beta_mat[-i, i] <- fit$edges
    }

    adj_mat <- matrix(as.integer((abs(beta_mat) + t(abs(beta_mat))) > 0), nrow = nrow(beta_mat))
    n_edge <- sum(adj_mat) / 2
    n_edge.vec <- c(n_edge.vec, n_edge)
    adj_mat.list[[ii]] <- adj_mat
    if(is_print){
      print(lambda_final.vec)
      print(c(prior.a, prior.b, n_edge, instability))
    }
  }

  df_summary <- data.frame(prior_a = prior_a.vec, prior_b = prior_b.vec,
                           instability = instability.vec,
                           n_edge = n_edge.vec)

  tmp <- df_summary
  tmp$index <- c(1: nrow(tmp))
  tmp <- tmp[tmp$instability <= 0.05, ]
  if(is_sparse_graph){
    tmp <- tmp[tmp$n_edge < p * (p - 1) / 4, ]
  }
  tmp <- tmp[order(tmp$n_edge, decreasing = T), ]
  tmp_index <- tmp[1, "index"]
  adj_mat <- adj_mat.list[[tmp_index]]

  res <- list(summary = df_summary, adj_mat.list = adj_mat.list, adj_mat = adj_mat)
  return(res)
}


poisson_transform <- function(df, alpha.vec = NULL, n_alpha = 100){

  if(is.null(alpha.vec)){
    alpha.vec <- seq(0, 1, length.out = n_alpha)
  }

  res.df <- data.frame(matrix(nrow = nrow(df), ncol = 0))
  for(ii in 1:ncol(df)){
    tmp_col <- df[, ii]
    mean_var_diff.vec <- c()
    for(jj in 1:length(alpha.vec)){
      tmp_out <- round(tmp_col ^ alpha.vec[jj])
      mean_var_diff.vec <- c(mean_var_diff.vec, abs(mean(tmp_out) - var(tmp_out)))
    }
    alpha_opt <- alpha.vec[order(mean_var_diff.vec)[1]]
    tmp_out <- round(tmp_col ^ alpha_opt)
    res.df <- cbind(res.df, data.frame(X = tmp_out))
  }
  colnames(res.df) <- colnames(df)
  return(res.df)
}


prob_predicting <- function(adj_mat.list){
  adj_mat.array <- array(dim = c(nrow(adj_mat.list[[1]]), ncol(adj_mat.list[[1]]), length(adj_mat.list)))
  for(ii in 1:length(adj_mat.list)){
    adj_mat.array[, , ii] <- adj_mat.list[[ii]]
  }
  n_edge <- apply(adj_mat.array, 3, sum)
  adj_mat.array <- adj_mat.array[, , order(n_edge, decreasing = T)]
  res <- apply(adj_mat.array, c(1, 2), function(x) ifelse(sum(x) > 0, max(c(1:length(x))[x > 0]), 0))
  res <- res / length(adj_mat.list)
  return(res)
}


#---- Real-word data: TCGA UECE ----
# read data: TCGA endometrial cancer gene mutation data with 59 genes
df_tcgauece <- read.csv("data\\endometrial_cancer.csv", header = T,sep = "\t")
df_tcgauece <- df_tcgauece[, -1]
df_tcgauecePois <- df_tcgauece

for(ii in 1:ncol(df_tcgauece)){
  df_tcgauece[, ii] <- as.integer(df_tcgauece[, ii] != 0)
}

ptm <- proc.time()
allen_method(df, B = 100, n_lambda = 100)
proc.time() - ptm
ptm <- proc.time()
for(ii in 1:10){
  my_method_given_prior(df, n_lambda = 100, prior.a = 1, prior.b = 45)
}
proc.time() - ptm

res_allen_tcgauece <- allen_method(df_tcgauece, B = 100, n_lambda = 300)
res2_allen_tcgauece <- allen_method(df_tcgauece, B = 100, n_lambda = 300)
res_my_tcgauece <- my_prior_stability_search_samesubsample(df_tcgauece,  n_lambda = 100, B = 100,
                                                           prior_a.vec = 1,
                                                           prior_b.vec = c(5, 10:35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
plot(res_allen_tcgauece$summary$n_edge, res_allen_tcgauece$summary$instability,
     xlim = c(0, 400), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res2_allen_tcgauece$summary$n_edge, res2_allen_tcgauece$summary$instability, col = "red")
points(res_my_tcgauece$summary$n_edge, res_my_tcgauece$summary$instability, col = "red")
abline(h = 0.05)

plot(NULL, xlim = c(0, 13), ylim = c(0, 20), xlab = "degree", ylab = "Freq",
     main = "TCGA.UECE distribution comparison LHGM vs LLGM")
for(ii in 1:1){
  tmp_df <- data.frame(table(apply(res_allen$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:13){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_my_tcgauece$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:13){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")
}

adj_mat_my_tcgauece <- res_my_tcgauece$adj_mat
rownames(adj_mat_my_tcgauece) <- colnames(df_tcgauece)
colnames(adj_mat_my_tcgauece) <- colnames(df_tcgauece)
write.csv(adj_mat_my_tcgauece, "adj_mat_my_tcgauece.csv")

adj_mat_allen_tcgauece <- res_allen_tcgauece$adj_mat
rownames(adj_mat_allen_tcgauece) <- colnames(df_tcgauece)
colnames(adj_mat_allen_tcgauece) <- colnames(df_tcgauece)
write.csv(adj_mat_allen_tcgauece, "adj_mat_allen_tcgauece.csv")


# Poisson
df_tcgauecePois_transformed <- poisson_transform(df_tcgauecePois)

res_allen_tcgauecePois <- allen_method(df_tcgauecePois_transformed, B = 100, n_lambda = 300, glm_family = "poisson")
res_my_tcgauecePois <- my_prior_stability_search_samesubsample(
  df_tcgauecePois_transformed, glm_family = "poisson", n_lambda = 100, B = 100,
  prior_a.vec = 1,
  prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgauecePois$summary$n_edge, res_allen_tcgauecePois$summary$instability,
     xlim = c(0,400), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgauecePois$summary$n_edge, res_my_tcgauecePois$summary$instability, col = "red")
abline(h=0.05)

adj_mat_my_tcgauecePois <- res_my_tcgauecePois$adj_mat
rownames(adj_mat_my_tcgauecePois) <- colnames(df_tcgauecePois)
colnames(adj_mat_my_tcgauecePois) <- colnames(df_tcgauecePois)
write.csv(adj_mat_my_tcgauecePois, "adj_mat_my_tcgauecePois.csv")

adj_mat_allen_tcgauecePois <- res_allen_tcgauecePois$adj_mat
rownames(adj_mat_allen_tcgauecePois) <- colnames(df_tcgauecePois)
colnames(adj_mat_allen_tcgauecePois) <- colnames(df_tcgauecePois)
write.csv(adj_mat_allen_tcgauecePois, "adj_mat_allen_tcgauecePois.csv")

plot(NULL, xlim = c(0, 20), ylim = c(0, 15), xlab = "degree", ylab = "Freq",
     main = "TCGA.UECE distribution comparison LHGM vs LLGM")
for(ii in 1:1){
  tmp_df <- data.frame(table(apply(adj_mat_allen_tcgauecePois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 0:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")
  tmp_df <- data.frame(table(apply(adj_mat_my_tcgauecePois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
}

dat <- rbind(data.frame(xx = apply(adj_mat_allen_tcgauecePois, 2, sum), yy = "LLGM"),
             data.frame(xx = apply(adj_mat_my_tcgauecePois, 2, sum), yy = "LHGM"))
ggplot(dat,aes(x=xx, fill = facotr(yy))) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="group", values=c("red","darkgray"),labels=c("a","b")) +
  xlab("degree") + ylab("#nodes")
ggplot(dat, aes(x=xx, fill = yy)) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="Method",values=c("red","blue"),labels=c("LHGM","LLGM"))

ggplot(dat,aes(x=xx))+
  geom_histogram(data=subset(dat,yy == "LHGM"),aes(fill=yy),alpha=0.8, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),aes(fill=yy),alpha=0.2, binwidth = 1)+
  scale_fill_manual(name="method", values=c("pink","blue"),labels=c("LHGM","LLGM"))+
  xlab("degree") + ylab("#nodes")

# XMRF
library(XMRF)
res_tcgauece_xmrf <- XMRF(t(df_tcgauecePois_transformed), method = "LPGM", stability = "star", nlams = 300, beta = 0.1)

# NegBin
res_allen_tcgaueceNB <- allen_method(df_tcgauecePois, B = 10, n_lambda = 100, glm_family = "negbin")
res_my_tcgaueceNB <- my_prior_stability_search_samesubsample(df_tcgauecePois, B = 100, n_lambda = 100, glm_family = "negbin",
                                                             prior_a.vec = 1,
                                                             prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgaueceNB$summary$n_edge, res_allen_tcgaueceNB$summary$instability,
     xlim = c(0,400), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgaueceNB$summary$n_edge, res_my_tcgaueceNB$summary$instability, col = "red")
abline(h=0.05)

#---- Real-word data: TCGA UECE MSKCC ----
# read data: TCGA endometrial cancer gene mutation data with 59 genes
df_tcgauecemskcc <- read.csv("data\\UCEC_bmi_gene_mut_Elsie_onlyMSKCCgene.txt", header = T,sep = "\t", stringsAsFactors = F)
df_tcgauecemskcc <- df_tcgauecemskcc[-c(1, 2), ]
gene_names <- df_tcgauecemskcc[, 1]
df_tcgauecemskcc <- df_tcgauecemskcc[, -1]
df_tcgauecemskcc <- data.frame(t(df_tcgauecemskcc))
colnames(df_tcgauecemskcc) <- gene_names
df_tcgauecemskccPois <- df_tcgauecemskcc

var.vec <- c()
mean.vec <- c()
var_bin.vec <- c()
mean_bin.vec <- c()
for(ii in 1:ncol(df_tcgauecemskccPois)){
  var.vec <- c(var.vec, var(df_tcgauecemskccPois[, ii]))
  mean.vec <- c(mean.vec, mean(df_tcgauecemskccPois[, ii]))
  var_bin.vec <- c(var_bin.vec, var(df_tcgauecemskccPois[, ii] != 0))
  mean_bin.vec <- c(mean_bin.vec, mean(df_tcgauecemskccPois[, ii] != 0))
}
tcgauecemskcc_gene_summary.df <- data.frame(gene = colnames(df_tcgauecemskccPois), var = var.vec, mean = mean.vec, var_bin = var_bin.vec, mean_bin = mean_bin.vec)
tcgauecemskcc_gene_summary.df$is_in_59 <- tcgauecemskcc_gene_summary.df$gene %in% colnames(df_tcgauecePois)
tcgauecemskcc_gene_summary.df <- tcgauecemskcc_gene_summary.df[order(tcgauecemskcc_gene_summary.df$mean_bin, decreasing = T), ]
tcgauecemskcc_gene_summary.df$rank <- c(1:nrow(tcgauecemskcc_gene_summary.df))
tcgauecemskcc_gene_summary.df[tcgauecemskcc_gene_summary.df$is_in_59, ]
df_tcgauecemskccPois <- df_tcgauecemskccPois[, order(var.vec, decreasing = T)[1:93]]
# Want to choose the largest 90 variance, but here 89-93 are tied, so we choose 93.
# Also, the 91st is JAK1, which is known very important.

df_tcgauecemskccPois_transformed <- poisson_transform(df_tcgauecemskccPois)

res_allen_tcgauecemskccPois <- allen_method(df_tcgauecemskccPois_transformed, B = 100, n_lambda = 300, glm_family = "poisson")
res_my_tcgauecemskccPois <- my_prior_stability_search_samesubsample(
  df_tcgauecemskccPois_transformed, glm_family = "poisson", n_lambda = 100, B = 100,
  prior_a.vec = 1,
  prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgauecemskccPois$summary$n_edge, res_allen_tcgauecemskccPois$summary$instability,
     xlim = c(0,400), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgauecemskccPois$summary$n_edge, res_my_tcgauecemskccPois$summary$instability, col = "red")
abline(h=0.05)

adj_mat_my_tcgauecemskccPois <- res_my_tcgauecemskccPois$adj_mat
rownames(adj_mat_my_tcgauecemskccPois) <- colnames(df_tcgauecemskccPois)
colnames(adj_mat_my_tcgauecemskccPois) <- colnames(df_tcgauecemskccPois)
write.csv(adj_mat_my_tcgauecemskccPois, "adj_mat_my_tcgauecemskccPois.csv")

adj_mat_allen_tcgauecemskccPois <- res_allen_tcgauecemskccPois$adj_mat
rownames(adj_mat_allen_tcgauecemskccPois) <- colnames(df_tcgauecemskccPois)
colnames(adj_mat_allen_tcgauecemskccPois) <- colnames(df_tcgauecemskccPois)
write.csv(adj_mat_allen_tcgauecemskccPois, "adj_mat_allen_tcgauecemskccPois.csv")

plot(NULL, xlim = c(0, 20), ylim = c(0, 15), xlab = "degree", ylab = "Freq",
     main = "TCGA.UECE distribution comparison LHGM vs LLGM")
for(ii in 1:1){
  tmp_df <- data.frame(table(apply(adj_mat_allen_tcgauecemskccPois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 0:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")
  tmp_df <- data.frame(table(apply(adj_mat_my_tcgauecemskccPois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
}

dat <- rbind(data.frame(xx = apply(adj_mat_allen_tcgauecemskccPois, 2, sum), yy = "LLGM"),
             data.frame(xx = apply(adj_mat_my_tcgauecemskccPois, 2, sum), yy = "LHGM"))
ggplot(dat,aes(x=xx, fill = facotr(yy))) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="group", values=c("red","darkgray"),labels=c("a","b")) +
  xlab("degree") + ylab("#nodes")
ggplot(dat, aes(x=xx, fill = yy)) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="Method",values=c("red","blue"),labels=c("LHGM","LLGM"))

ggplot(dat,aes(x=xx))+
  geom_histogram(data=subset(dat,yy == "LHGM"),aes(fill=yy),alpha=0.8, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),aes(fill=yy),alpha=0.2, binwidth = 1)+
  scale_fill_manual(name="method", values=c("pink","blue"),labels=c("LHGM","LLGM"))+
  xlab("degree") + ylab("#nodes")

#---- Real-word data: TCGA UECE Filter with BMI ----
# read data: TCGA endometrial cancer gene mutation data with all genes and then filtering
df_tcgaueceallgene <- read.csv("data\\UCEC_bmi_gene_mut_Elsie_allgene.txt", header = T,sep = "\t", stringsAsFactors = F)
df_tcgaueceallgene <- df_tcgaueceallgene[-c(1, 2), ]
gene_names <- df_tcgaueceallgene[, 1]
df_tcgaueceallgene <- df_tcgaueceallgene[, -1]
df_tcgaueceallgene <- data.frame(t(df_tcgaueceallgene))
colnames(df_tcgaueceallgene) <- gene_names
df_tcgaueceallgene <- subset(df_tcgaueceallgene, select = -Unknown)

genelist_nature2016 <- c("PTEN", "PIK3CA", "PIK3R1", "CTNNB1", "TP53", "KRAS",
                         "FBXW7", "SPOP", "CTCF", "ARID1A", "PPP2R1A", "CCND1",
                         "FGFR2", "ARHGAP35", "CHD4", "ZFHX3", "SOX17", "ERBB3",
                         "ARID5B", "NRAS", "NFE2L2", "MAX", "SOS1", "SGK1", "ERBB2",
                         "BCOR", "ESR1", "RRAS2", "SIN3A", "MYCN", "AFMID", "MECOM",
                         "FOXA2", "ALPK2", "AKT1", "METTL14", "SERHL2", "WDR45",
                         "U2AF1", "TAB3", "ADNP", "MARK3", "EDC4", "DICER1", "RBM39",
                         "POLE", "RNF43", "JAK1", "NRIP1")

var.vec <- c()
mean.vec <- c()
var_bin.vec <- c()
mean_bin.vec <- c()
for(ii in 1:ncol(df_tcgaueceallgene)){
  var.vec <- c(var.vec, var(df_tcgaueceallgene[, ii]))
  mean.vec <- c(mean.vec, mean(df_tcgaueceallgene[, ii]))
  var_bin.vec <- c(var_bin.vec, var(df_tcgaueceallgene[, ii] != 0))
  mean_bin.vec <- c(mean_bin.vec, mean(df_tcgaueceallgene[, ii] != 0))
}
tcgaueceallgene_gene_summary.df <- data.frame(gene = colnames(df_tcgaueceallgene), var = var.vec, mean = mean.vec, var_bin = var_bin.vec, mean_bin = mean_bin.vec)
tcgaueceallgene_gene_summary.df$is_in_59 <- tcgaueceallgene_gene_summary.df$gene %in% colnames(df_tcgauecePois)
tcgaueceallgene_gene_summary.df$is_in_49 <- tcgaueceallgene_gene_summary.df$gene %in% genelist_nature2016
tcgaueceallgene_gene_summary.df <- tcgaueceallgene_gene_summary.df[order(tcgaueceallgene_gene_summary.df$mean_bin, decreasing = T), ]
tcgaueceallgene_gene_summary.df$rank <- c(1:nrow(tcgaueceallgene_gene_summary.df))
tcgaueceallgene_gene_summary.df[tcgaueceallgene_gene_summary.df$is_in_59 , ]

df_tcgauecefiltered <- df_tcgaueceallgene[, order(mean_bin.vec, decreasing = T)[1:100]]
#df_tcgauecefiltered <- df_tcgaueceallgene[, colnames(df_tcgaueceallgene) %in% genelist_nature2016]

df_tcgauecefilteredPois <- df_tcgauecefiltered
df_tcgauecefilteredPois_transformed <- poisson_transform(df_tcgauecefilteredPois)

res_allen_tcgauecefilteredPois <- allen_method(df_tcgauecefilteredPois_transformed, B = 100, n_lambda = 300, glm_family = "poisson")
res_my_tcgauecefilteredPois <- my_prior_stability_search_samesubsample(
  df_tcgauecefilteredPois_transformed, glm_family = "poisson", n_lambda = 100, B = 100,
  prior_a.vec = 1,
  prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgauecefilteredPois$summary$n_edge, res_allen_tcgauecefilteredPois$summary$instability,
     xlim = c(0,800), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgauecefilteredPois$summary$n_edge, res_my_tcgauecefilteredPois$summary$instability, col = "red")
abline(h=0.05)

adj_mat_my_tcgauecefilteredPois <- res_my_tcgauecefilteredPois$adj_mat
rownames(adj_mat_my_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
colnames(adj_mat_my_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
write.csv(adj_mat_my_tcgauecefilteredPois, "adj_mat_my_tcgauecefilteredPois.csv")

adj_mat_allen_tcgauecefilteredPois <- res_allen_tcgauecefilteredPois$adj_mat
rownames(adj_mat_allen_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
colnames(adj_mat_allen_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
write.csv(adj_mat_allen_tcgauecefilteredPois, "adj_mat_allen_tcgauecefilteredPois.csv")

plot(NULL, xlim = c(0, 20), ylim = c(0, 15), xlab = "degree", ylab = "Freq",
     main = "TCGA.UECE distribution comparison LHGM vs LLGM")
for(ii in 1:1){
  tmp_df <- data.frame(table(apply(adj_mat_allen_tcgauecefilteredPois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 0:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")
  tmp_df <- data.frame(table(apply(adj_mat_my_tcgauecefilteredPois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
}

dat <- rbind(data.frame(xx = apply(adj_mat_allen_tcgauecefilteredPois, 2, sum), yy = "LLGM"),
             data.frame(xx = apply(adj_mat_my_tcgauecefilteredPois, 2, sum), yy = "LHGM"))
ggplot(dat,aes(x=xx, fill = facotr(yy))) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="group", values=c("red","darkgray"),labels=c("a","b")) +
  xlab("degree") + ylab("#nodes")
ggplot(dat, aes(x=xx, fill = yy)) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="Method",values=c("red","blue"),labels=c("LHGM","LLGM"))

ggplot(dat,aes(x=xx))+
  geom_histogram(data=subset(dat,yy == "LHGM"),aes(fill=yy),alpha=0.8, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),aes(fill=yy),alpha=0.2, binwidth = 1)+
  scale_fill_manual(name="method", values=c("pink","blue"),labels=c("LHGM","LLGM"))+
  xlab("degree") + ylab("#nodes")


#---- Real-word data: TCGA UECE Filter without BMI ----
# read data: TCGA endometrial cancer gene mutation data with all genes and then filtering
df_tcgaueceallgene <- read.csv("data\\UCEC_nobmi_gene_mut_Elsie_allgene.txt", header = T,sep = "\t", stringsAsFactors = F)
df_tcgaueceallgene <- df_tcgaueceallgene[-c(1), ]
gene_names <- df_tcgaueceallgene[, 1]
df_tcgaueceallgene <- df_tcgaueceallgene[, -1]
df_tcgaueceallgene <- data.frame(t(df_tcgaueceallgene))
colnames(df_tcgaueceallgene) <- gene_names
df_tcgaueceallgene <- subset(df_tcgaueceallgene, select = -Unknown)

var.vec <- c()
mean.vec <- c()
var_bin.vec <- c()
mean_bin.vec <- c()
for(ii in 1:ncol(df_tcgaueceallgene)){
  var.vec <- c(var.vec, var(df_tcgaueceallgene[, ii]))
  mean.vec <- c(mean.vec, mean(df_tcgaueceallgene[, ii]))
  var_bin.vec <- c(var_bin.vec, var(df_tcgaueceallgene[, ii] != 0))
  mean_bin.vec <- c(mean_bin.vec, mean(df_tcgaueceallgene[, ii] != 0))
}
tcgaueceallgene_gene_summary.df <- data.frame(gene = colnames(df_tcgaueceallgene), var = var.vec, mean = mean.vec, var_bin = var_bin.vec, mean_bin = mean_bin.vec)
tcgaueceallgene_gene_summary.df$is_in_59 <- tcgaueceallgene_gene_summary.df$gene %in% colnames(df_tcgauecePois)
tcgaueceallgene_gene_summary.df <- tcgaueceallgene_gene_summary.df[order(tcgaueceallgene_gene_summary.df$mean_bin, decreasing = T), ]
tcgaueceallgene_gene_summary.df$rank <- c(1:nrow(tcgaueceallgene_gene_summary.df))
tcgaueceallgene_gene_summary.df[tcgaueceallgene_gene_summary.df$is_in_59 , ]
tcgaueceallgene_gene_summary.df[1:100, ]

df_tcgauecefiltered <- df_tcgaueceallgene[, order(mean_bin.vec, decreasing = T)[1:68]]
#df_tcgauecefiltered <- df_tcgaueceallgene[, colnames(df_tcgaueceallgene) %in% genelist_nature2016]

df_tcgauecefilteredPois <- df_tcgauecefiltered
df_tcgauecefilteredPois_transformed <- poisson_transform(df_tcgauecefilteredPois)

res_allen_tcgauecefilteredPois <- allen_method(df_tcgauecefilteredPois_transformed, B = 100, n_lambda = 300, glm_family = "poisson")
res_my_tcgauecefilteredPois <- my_prior_stability_search_samesubsample(
  df_tcgauecefilteredPois_transformed, glm_family = "poisson", n_lambda = 100, B = 100,
  prior_a.vec = 1,
  prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgauecefilteredPois$summary$n_edge, res_allen_tcgauecefilteredPois$summary$instability,
     xlim = c(0,800), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgauecefilteredPois$summary$n_edge, res_my_tcgauecefilteredPois$summary$instability, col = "red")
abline(h=0.05)

adj_mat_my_tcgauecefilteredPois <- res_my_tcgauecefilteredPois$adj_mat
rownames(adj_mat_my_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
colnames(adj_mat_my_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
write.csv(adj_mat_my_tcgauecefilteredPois, "adj_mat_my_tcgauecefilteredPois.csv")

adj_mat_allen_tcgauecefilteredPois <- res_allen_tcgauecefilteredPois$adj_mat
rownames(adj_mat_allen_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
colnames(adj_mat_allen_tcgauecefilteredPois) <- colnames(df_tcgauecefilteredPois)
write.csv(adj_mat_allen_tcgauecefilteredPois, "adj_mat_allen_tcgauecefilteredPois.csv")

plot(NULL, xlim = c(0, 20), ylim = c(0, 15), xlab = "Degree", ylab = "Freq",
     main = "TCGA.UCEC distribution comparison LHGM vs LLGM")
for(ii in 1:1){
  tmp_df <- data.frame(table(apply(adj_mat_allen_tcgauecefilteredPois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 0:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")
  tmp_df <- data.frame(table(apply(adj_mat_my_tcgauecefilteredPois, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:20){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
}

dat <- rbind(data.frame(xx = apply(adj_mat_allen_tcgauecefilteredPois, 2, sum), yy = "LLGM"),
             data.frame(xx = apply(adj_mat_my_tcgauecefilteredPois, 2, sum), yy = "LHGM"))
ggplot(dat,aes(x=xx, fill = facotr(yy))) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="group", values=c("red","darkgray"),labels=c("a","b")) +
  xlab("degree") + ylab("#nodes")
ggplot(dat, aes(x=xx, fill = yy)) +
  geom_histogram(data=subset(dat,yy == "LHGM"),fill = "red", alpha = 0.2, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),fill = "blue", alpha = 0.2, binwidth = 1) +
  scale_fill_manual(name="Method",values=c("red","blue"),labels=c("LHGM","LLGM"))

ggplot(dat,aes(x=xx))+
  geom_histogram(data=subset(dat,yy == "LHGM"),aes(fill=yy),alpha=0.8, binwidth = 1)+
  geom_histogram(data=subset(dat,yy == "LLGM"),aes(fill=yy),alpha=0.2, binwidth = 1)+
  scale_fill_manual(name="method", values=c("pink","blue"),labels=c("LHGM","LLGM"))+
  xlab("degree") + ylab("#nodes")

ggplot(dat,aes(x=xx)) + geom_histogram(aes(fill = yy), binwidth = 1) +
  facet_grid(yy ~ ., scales = "fixed") +
  theme(legend.position = "none",
        axis.text=element_text(size=16),
        axis.title=element_text(size=20),
        strip.text.y = element_text(size = 18, face = "bold.italic")
        ) +
  xlab("degree") + ylab("#nodes")

#---- Real-word data: TCGA Colon ----

# read data: TCGAColon
df_tcgacolon <- read.csv("data\\TCGAcolon_bmi_gene_mut_Elsie.txt", header = T,sep = "\t", stringsAsFactors = F, row.names = 1)
df_tcgacolon <- data.frame(t(df_tcgacolon))

df_tcgacolonPois <- df_tcgacolon

for(ii in 1:ncol(df_tcgacolon)){
  df_tcgacolon[, ii] <- as.integer(df_tcgacolon[, ii] != 0)
}

var.vec <- c()
for(ii in 1:ncol(df_tcgacolon)){
  var.vec <- c(var.vec, var(df_tcgacolon[, ii]))
}
df_tcgacolon <- df_tcgacolon[, order(var.vec, decreasing = T)[1:60]]

var.vec <- c()
for(ii in 1:ncol(df_tcgacolonPois)){
  var.vec <- c(var.vec, var(df_tcgacolonPois[, ii]))
}
df_tcgacolonPois <- df_tcgacolonPois[, order(var.vec, decreasing = T)[1:60]]


res_allen_tcgacolon <- allen_method(df_tcgacolon, B = 100, n_lambda = 300)
res_my_tcgacolon <- my_prior_stability_search_samesubsample(
  df_tcgacolon, n_lambda = 100, B = 100,
  prior_a.vec = 1,
  prior_b.vec = c(5, 10, 15, 20:35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgacolon$summary$n_edge, res_allen_tcgacolon$summary$instability,
     xlim = c(0,400), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgacolon$summary$n_edge, res_my_tcgacolon$summary$instability, col = "red")
abline(h=0.05)

adj_mat_my_tcgacolon <- res_my_tcgacolon$adj_mat
rownames(adj_mat_my_tcgacolon) <- colnames(df_tcgacolon)
colnames(adj_mat_my_tcgacolon) <- colnames(df_tcgacolon)
write.csv(adj_mat_my_tcgacolon, "adj_mat_my_tcgacolon.csv")

adj_mat_allen_tcgacolon <- res2_allen_tcgacolon$adj_mat
rownames(adj_mat_allen_tcgacolon) <- colnames(df_tcgacolon)
colnames(adj_mat_allen_tcgacolon) <- colnames(df_tcgacolon)
write.csv(adj_mat_allen_tcgacolon, "adj_mat_allen_tcgacolon.csv")

# Poisson

df_tcgacolonPois_transformed <- poisson_transform(df_tcgacolonPois)

res_allen_tcgacolonPois <- allen_method(df_tcgacolonPois_transformed, B = 100, n_lambda = 300, glm_family = "poisson")
res_my_tcgacolonPois <- my_prior_stability_search_samesubsample(
  df_tcgacolonPois_transformed, glm_family = "poisson", n_lambda = 100, B = 100,
  prior_a.vec = 1,
  prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

plot(res_allen_tcgacolonPois$summary$n_edge, res_allen_tcgacolonPois$summary$instability,
     xlim = c(0,400), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_tcgacolonPois$summary$n_edge, res_my_tcgacolonPois$summary$instability, col = "red")
abline(h=0.05)

adj_mat_my_tcgacolonPois <- res_my_tcgacolonPois$adj_mat
rownames(adj_mat_my_tcgacolonPois) <- colnames(df_tcgacolonPois)
colnames(adj_mat_my_tcgacolonPois) <- colnames(df_tcgacolonPois)
write.csv(adj_mat_my_tcgacolonPois, "adj_mat_my_tcgacolonPois.csv")

adj_mat_allen_tcgacolonPois <- res_allen_tcgacolonPois$adj_mat
rownames(adj_mat_allen_tcgacolonPois) <- colnames(df_tcgacolonPois)
colnames(adj_mat_allen_tcgacolonPois) <- colnames(df_tcgacolonPois)
write.csv(adj_mat_allen_tcgacolonPois, "adj_mat_allen_tcgacolonPois.csv")


#---- Simulation: 6 clu (mvp) unequal noise unequal signal ----

n_node_perclu <- 6
n_clu <- 10
n_sample <- 200
A_truth_mvp6clu <- matrix(rep(0, (n_node_perclu * n_clu) ^ 2), nrow = n_node_perclu * n_clu)
for(ii in 1:n_clu){
  for(jj in 1:n_node_perclu){
    for(kk in 1:n_node_perclu){
      A_truth_mvp6clu[(ii - 1) * n_node_perclu + jj, (ii - 1) * n_node_perclu + kk] <- 1
      A_truth_mvp6clu[(ii - 1) * n_node_perclu + kk, (ii - 1) * n_node_perclu + jj] <- 1
    }
  }
}
diag(A_truth_mvp6clu) <- 0

#   10 times
res_allen_mvp6clu.list <- list()
res_my_mvp6clu.list <- list()
for(cc in 1:10){
  cat("*********** sim mvp6clu, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6clu <- A_truth_mvp6clu * 0.2
  diag(corr_mvp6clu) <- 1

  df_sim_mvp6clu <- data.frame(genPoisNor(n_sample, 0, n_node_perclu * n_clu, corr_mvp6clu, rep(1:n_clu, each = n_node_perclu), c(), c()))
  for(ii in 1:(n_clu * n_node_perclu)){
    df_sim_mvp6clu[, ii] <- df_sim_mvp6clu[, ii] + rpois(n_sample, ceiling(ii / n_node_perclu))
  }

  res_allen_mvp6clu.list[[cc]] <- allen_method(df_sim_mvp6clu, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6clu.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6clu, glm_family = "poisson",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6clu.list[[1]]$summary$n_edge, res_allen_mvp6clu.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6clu.list[[1]]$summary$n_edge, res_my_mvp6clu.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6clu distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6clu.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6clu.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6clu.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6clu.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6clu)){
    for(kk in 1:nrow(A_truth_mvp6clu)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6clu[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6clu.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6clu.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_mvp6clu <- df_accu
df_accu_mvp6clu
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6clu, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvp6clu, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvp6clu.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvp6clu.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")

#---- Simulation: 6 clu (mvp) unequal signal equal noise ----

n_node_perclu <- 6
n_clu <- 10
n_sample <- 200
A_truth_mvp6cluUSEN <- matrix(rep(0, (n_node_perclu * n_clu) ^ 2), nrow = n_node_perclu * n_clu)
for(ii in 1:n_clu){
  for(jj in 1:n_node_perclu){
    for(kk in 1:n_node_perclu){
      A_truth_mvp6cluUSEN[(ii - 1) * n_node_perclu + jj, (ii - 1) * n_node_perclu + kk] <- 1
      A_truth_mvp6cluUSEN[(ii - 1) * n_node_perclu + kk, (ii - 1) * n_node_perclu + jj] <- 1
    }
  }
}
diag(A_truth_mvp6cluUSEN) <- 0

#   10 times
res_allen_mvp6cluUSEN.list <- list()
res_my_mvp6cluUSEN.list <- list()
for(cc in 1:10){
  cat("*********** sim mvp6cluUSEN, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6cluUSEN <- A_truth_mvp6cluUSEN * 0.2
  diag(corr_mvp6cluUSEN) <- 1

  df_sim_mvp6cluUSEN <- data.frame(genPoisNor(n_sample, 0, n_node_perclu * n_clu, corr_mvp6cluUSEN, rep(1:n_clu, each = n_node_perclu), c(), c()))
  for(ii in 1:(n_clu * n_node_perclu)){
    df_sim_mvp6cluUSEN[, ii] <- df_sim_mvp6cluUSEN[, ii] + rpois(n_sample, 5)
  }

  res_allen_mvp6cluUSEN.list[[cc]] <- allen_method(df_sim_mvp6cluUSEN, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6cluUSEN.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6cluUSEN, glm_family = "poisson",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6cluUSEN.list[[1]]$summary$n_edge, res_allen_mvp6cluUSEN.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6cluUSEN.list[[1]]$summary$n_edge, res_my_mvp6cluUSEN.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6cluUSEN distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6cluUSEN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6cluUSEN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6cluUSEN.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6cluUSEN.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6cluUSEN)){
    for(kk in 1:nrow(A_truth_mvp6cluUSEN)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6cluUSEN[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6cluUSEN.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6cluUSEN.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
roc(y_true, y_pred_prob_allen, plot = T, legacy.axes = T, percent = T,
    xlab = "False Positive Perc", ylab = "True Positive Perc", col = "#377eb8", lwd = 4, print.auc = T)
plot.roc(y_true, y_pred_prob_my, percent = T, col = "#4daf4a", lwd = 4, print.auc = T, add = T, print.auc.y=40)

df_accu_mvp6cluUSEN <- df_accu
df_accu_mvp6cluUSEN
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6cluUSEN, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvp6cluUSEN, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvp6cluUSEN.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvp6cluUSEN.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")

#---- Simulation: 6 clu (mvp) equal signal unequal noise ----

n_node_perclu <- 6
n_clu <- 10
n_sample <- 200
A_truth_mvp6cluESUN <- matrix(rep(0, (n_node_perclu * n_clu) ^ 2), nrow = n_node_perclu * n_clu)
for(ii in 1:n_clu){
  for(jj in 1:n_node_perclu){
    for(kk in 1:n_node_perclu){
      A_truth_mvp6cluESUN[(ii - 1) * n_node_perclu + jj, (ii - 1) * n_node_perclu + kk] <- 1
      A_truth_mvp6cluESUN[(ii - 1) * n_node_perclu + kk, (ii - 1) * n_node_perclu + jj] <- 1
    }
  }
}
diag(A_truth_mvp6cluESUN) <- 0

#   10 times
res_allen_mvp6cluESUN.list <- list()
res_my_mvp6cluESUN.list <- list()
for(cc in 1:10){
  cat("*********** sim mvp6cluESUN, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6cluESUN <- A_truth_mvp6cluESUN * 0.2
  diag(corr_mvp6cluESUN) <- 1

  df_sim_mvp6cluESUN <- data.frame(genPoisNor(n_sample, 0, n_node_perclu * n_clu, corr_mvp6cluESUN, rep(5, n_clu * n_node_perclu), c(), c()))
  for(ii in 1:(n_clu * n_node_perclu)){
    df_sim_mvp6cluESUN[, ii] <- df_sim_mvp6cluESUN[, ii] + rpois(n_sample, ceiling(ii / n_node_perclu))
  }

  res_allen_mvp6cluESUN.list[[cc]] <- allen_method(df_sim_mvp6cluESUN, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6cluESUN.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6cluESUN, glm_family = "poisson",
                                                                           n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                           prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6cluESUN.list[[1]]$summary$n_edge, res_allen_mvp6cluESUN.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6cluESUN.list[[1]]$summary$n_edge, res_my_mvp6cluESUN.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6cluESUN distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6cluESUN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6cluESUN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6cluESUN.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6cluESUN.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6cluESUN)){
    for(kk in 1:nrow(A_truth_mvp6cluESUN)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6cluESUN[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6cluESUN.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6cluESUN.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
roc(y_true, y_pred_prob_allen, plot = T, legacy.axes = T, percent = T,
    xlab = "False Positive Perc", ylab = "True Positive Perc", col = "#377eb8", lwd = 4, print.auc = T)
plot.roc(y_true, y_pred_prob_my, percent = T, col = "#4daf4a", lwd = 4, print.auc = T, add = T, print.auc.y=40)

df_accu_mvp6cluESUN <- df_accu
df_accu_mvp6cluESUN
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6cluESUN, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvp6cluESUN, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvp6cluESUN.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvp6cluESUN.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")


#---- Simulation: 6 clu (mvp) equal signal equal noise ----

n_node_perclu <- 6
n_clu <- 10
n_sample <- 200
A_truth_mvp6clueqhigh <- matrix(rep(0, (n_node_perclu * n_clu) ^ 2), nrow = n_node_perclu * n_clu)
for(ii in 1:n_clu){
  for(jj in 1:n_node_perclu){
    for(kk in 1:n_node_perclu){
      A_truth_mvp6clueqhigh[(ii - 1) * n_node_perclu + jj, (ii - 1) * n_node_perclu + kk] <- 1
      A_truth_mvp6clueqhigh[(ii - 1) * n_node_perclu + kk, (ii - 1) * n_node_perclu + jj] <- 1
    }
  }
}
diag(A_truth_mvp6clueqhigh) <- 0

#   10 times
res_allen_mvp6clueqhigh.list <- list()
res_my_mvp6clueqhigh.list <- list()
for(cc in 1:10){
  cat("*********** sim mvp6clueqhigh, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6clueqhigh <- A_truth_mvp6clueqhigh * 0.2
  diag(corr_mvp6clueqhigh) <- 1

  df_sim_mvp6clueqhigh <- data.frame(genPoisNor(n_sample, 0, n_node_perclu * n_clu, corr_mvp6clueqhigh, rep(5, n_node_perclu * n_clu), c(), c()))
  for(ii in 1:(n_clu * n_node_perclu)){
    df_sim_mvp6clueqhigh[, ii] <- df_sim_mvp6clueqhigh[, ii] + rpois(n_sample, 5)
  }

  res_allen_mvp6clueqhigh.list[[cc]] <- allen_method(df_sim_mvp6clueqhigh, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6clueqhigh.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6clueqhigh, glm_family = "poisson",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6clueqhigh.list[[1]]$summary$n_edge, res_allen_mvp6clueqhigh.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6clueqhigh.list[[1]]$summary$n_edge, res_my_mvp6clueqhigh.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6clueqhigh distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6clueqhigh.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6clueqhigh.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6clueqhigh.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6clueqhigh.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6clueqhigh)){
    for(kk in 1:nrow(A_truth_mvp6clueqhigh)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6clueqhigh[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6clueqhigh.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6clueqhigh.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_mvp6clueqhigh <- df_accu
df_accu_mvp6clueqhigh
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6clueqhigh, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvp6clueqhigh, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvp6clueqhigh.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvp6clueqhigh.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")


#---- Simulation: 6 hub (mvp) unequal noise unequal signal ----

n_neighbor_hub <- 9
n_hub <- 6
n_sample <- 200
A_truth_mvp6hub <- matrix(rep(0, (n_hub * (n_neighbor_hub + 1)) ^ 2), nrow = n_hub * (n_neighbor_hub + 1))
for(ii in c(1:n_hub)){
  for(jj in c(1:n_neighbor_hub)){
    A_truth_mvp6hub[n_hub*n_neighbor_hub + ii, n_neighbor_hub * (ii - 1) + jj] <- 1
    A_truth_mvp6hub[n_neighbor_hub * (ii - 1) + jj, n_hub*n_neighbor_hub + ii] <- 1
  }
}
diag(A_truth_mvp6hub) <- 0

#   10 times
res_allen_mvp6hub.list <- list()
res_my_mvp6hub.list <- list()
for(cc in c(1:10)){
  cat("*********** sim mvp6hub, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6hub <- A_truth_mvp6hub * 0.3
  diag(corr_mvp6hub) <- 1

  df_sim_mvp6hub <- data.frame(genPoisNor(n_sample, 0, n_node_perclu * n_clu, corr_mvp6hub,
                                          #rep(1, n_node_perclu * n_clu),
                                          c(rep(1:n_hub, each = n_neighbor_hub), 1:n_hub),
                                          c(), c()))
  for(ii in c(1:n_hub)){
    for(jj in c(1:n_neighbor_hub)){
      df_sim_mvp6hub[, n_neighbor_hub * (ii - 1) + jj] <- df_sim_mvp6hub[, n_neighbor_hub * (ii - 1) + jj] + rpois(n_sample, ii)
    }
  }
  for(ii in c(1:n_hub)){
    df_sim_mvp6hub[, n_hub * n_neighbor_hub + ii] <- df_sim_mvp6hub[, n_hub * n_neighbor_hub + ii] + rpois(n_sample, ii)
  }

  res_allen_mvp6hub.list[[cc]] <- allen_method(df_sim_mvp6hub, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6hub.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6hub, glm_family = "poisson",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6hub.list[[1]]$summary$n_edge, res_allen_mvp6hub.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6hub.list[[1]]$summary$n_edge, res_my_mvp6hub.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_mvp6hub), as.vector(prob_predicting(res_allen_mvp6hub.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_mvp6hub), as.vector(prob_predicting(res_my_mvp6hub.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6hub distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6hub.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6hub.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6hub.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6hub.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6hub)){
    for(kk in 1:nrow(A_truth_mvp6hub)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6hub[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6hub.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6hub.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)

}
df_accu_mvp6hub <- df_accu
df_accu_mvp6hub
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6hub, "res_sim\\df_accu_6hubLogit.csv")
write.csv(A_truth_mvp6hub, "res_sim\\adj_mat_truth_6hubLogit.csv")
write.csv(res_my_mvp6hub.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6hubLogit.csv")
write.csv(res_allen_mvp6hub.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6hubLogit.csv")


#---- Simulation: 6 hub (mvp) unequal signal equal noise ----

n_neighbor_hub <- 9
n_hub <- 6
n_sample <- 200
A_truth_mvp6hubeqhigh <- matrix(rep(0, (n_hub * (n_neighbor_hub + 1)) ^ 2), nrow = n_hub * (n_neighbor_hub + 1))
for(ii in c(1:n_hub)){
  for(jj in c(1:n_neighbor_hub)){
    A_truth_mvp6hubeqhigh[n_hub*n_neighbor_hub + ii, n_neighbor_hub * (ii - 1) + jj] <- 1
    A_truth_mvp6hubeqhigh[n_neighbor_hub * (ii - 1) + jj, n_hub*n_neighbor_hub + ii] <- 1
  }
}
diag(A_truth_mvp6hubeqhigh) <- 0

#   10 times
res_allen_mvp6hubeqhigh.list <- list()
res_my_mvp6hubeqhigh.list <- list()
for(cc in c(1:1)){
  cat("*********** sim mvp6hubeqhigh, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6hubeqhigh <- A_truth_mvp6hubeqhigh * 0.3
  diag(corr_mvp6hubeqhigh) <- 1

  df_sim_mvp6hubeqhigh <- data.frame(genPoisNor(n_sample, 0, n_hub*(n_neighbor_hub + 1), corr_mvp6hubeqhigh,
                                          #rep(5, n_hub*(n_neighbor_hub + 1)),
                                          c(rep(1:n_hub, each = n_neighbor_hub), 1:n_hub),
                                          c(), c()))
  for(ii in 1:ncol(df_sim_mvp6hubeqhigh)){
    df_sim_mvp6hubeqhigh[, ii] <- df_sim_mvp6hubeqhigh[, ii] + rpois(n_sample, 5)
  }


  res_allen_mvp6hubeqhigh.list[[cc]] <- allen_method(df_sim_mvp6hubeqhigh, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6hubeqhigh.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6hubeqhigh, glm_family = "poisson",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6hubeqhigh.list[[1]]$summary$n_edge, res_allen_mvp6hubeqhigh.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6hubeqhigh.list[[1]]$summary$n_edge, res_my_mvp6hubeqhigh.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_mvp6hubeqhigh), as.vector(prob_predicting(res_allen_mvp6hubeqhigh.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_mvp6hubeqhigh), as.vector(prob_predicting(res_my_mvp6hubeqhigh.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6hubeqhigh distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6hubeqhigh.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6hubeqhigh.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6hubeqhigh.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6hubeqhigh.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6hubeqhigh)){
    for(kk in 1:nrow(A_truth_mvp6hubeqhigh)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6hubeqhigh[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6hubeqhigh.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6hubeqhigh.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)

}
df_accu_mvp6hubeqhigh <- df_accu
df_accu_mvp6hubeqhigh
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6hubeqhigh, "res_sim\\df_accu_6hubLogit.csv")
write.csv(A_truth_mvp6hubeqhigh, "res_sim\\adj_mat_truth_6hubLogit.csv")
write.csv(res_my_mvp6hubeqhigh.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6hubLogit.csv")
write.csv(res_allen_mvp6hubeqhigh.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6hubLogit.csv")


#---- Simulation: 6 hub (mvp) equal signal unequal noise ----

n_neighbor_hub <- 9
n_hub <- 6
n_sample <- 200
A_truth_mvp6hubESUN <- matrix(rep(0, (n_hub * (n_neighbor_hub + 1)) ^ 2), nrow = n_hub * (n_neighbor_hub + 1))
for(ii in c(1:n_hub)){
  for(jj in c(1:n_neighbor_hub)){
    A_truth_mvp6hubESUN[n_hub*n_neighbor_hub + ii, n_neighbor_hub * (ii - 1) + jj] <- 1
    A_truth_mvp6hubESUN[n_neighbor_hub * (ii - 1) + jj, n_hub*n_neighbor_hub + ii] <- 1
  }
}
diag(A_truth_mvp6hubESUN) <- 0

#   10 times
res_allen_mvp6hubESUN.list <- list()
res_my_mvp6hubESUN.list <- list()
for(cc in 6:10){
  cat("*********** sim mvp6hubESUN, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6hubESUN <- A_truth_mvp6hubESUN * 0.3
  diag(corr_mvp6hubESUN) <- 1

  df_sim_mvp6hubESUN <- data.frame(genPoisNor(n_sample, 0, n_hub*(n_neighbor_hub + 1), corr_mvp6hubESUN,
                                                rep(5, n_hub*(n_neighbor_hub + 1)),
                                                #c(rep(1:n_hub, each = n_neighbor_hub), 1:n_hub),
                                                c(), c()))
  for(ii in 1:ncol(df_sim_mvp6hubESUN)){
    df_sim_mvp6hubESUN[, ii] <- df_sim_mvp6hubESUN[, ii] + rpois(n_sample, 5)
  }
  for(ii in c(1:n_hub)){
    for(jj in c(1:n_neighbor_hub)){
      df_sim_mvp6hubESUN[, n_neighbor_hub * (ii - 1) + jj] <- df_sim_mvp6hubESUN[, n_neighbor_hub * (ii - 1) + jj] + rpois(n_sample, ii)
    }
    df_sim_mvp6hubESUN[, n_hub*n_neighbor_hub + ii] <- df_sim_mvp6hubESUN[, n_hub*n_neighbor_hub + ii] + rpois(n_sample, ii)
  }

  res_allen_mvp6hubESUN.list[[cc]] <- allen_method(df_sim_mvp6hubESUN, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6hubESUN.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6hubESUN, glm_family = "poisson",
                                                                             n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                             prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6hubESUN.list[[1]]$summary$n_edge, res_allen_mvp6hubESUN.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6hubESUN.list[[1]]$summary$n_edge, res_my_mvp6hubESUN.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_mvp6hubESUN), as.vector(prob_predicting(res_allen_mvp6hubESUN.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_mvp6hubESUN), as.vector(prob_predicting(res_my_mvp6hubESUN.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6hubESUN distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6hubESUN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6hubESUN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6hubESUN.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6hubESUN.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6hubESUN)){
    for(kk in 1:nrow(A_truth_mvp6hubESUN)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6hubESUN[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6hubESUN.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6hubESUN.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)

}
df_accu_mvp6hubESUN <- df_accu
df_accu_mvp6hubESUN
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6hubESUN, "res_sim\\df_accu_6hubLogit.csv")
write.csv(A_truth_mvp6hubESUN, "res_sim\\adj_mat_truth_6hubLogit.csv")
write.csv(res_my_mvp6hubESUN.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6hubLogit.csv")
write.csv(res_allen_mvp6hubESUN.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6hubLogit.csv")



#---- Simulation: 6 hub (mvp) equal signal equal noise ----

n_neighbor_hub <- 9
n_hub <- 6
n_sample <- 200
A_truth_mvp6hubESEN <- matrix(rep(0, (n_hub * (n_neighbor_hub + 1)) ^ 2), nrow = n_hub * (n_neighbor_hub + 1))
for(ii in c(1:n_hub)){
  for(jj in c(1:n_neighbor_hub)){
    A_truth_mvp6hubESEN[n_hub*n_neighbor_hub + ii, n_neighbor_hub * (ii - 1) + jj] <- 1
    A_truth_mvp6hubESEN[n_neighbor_hub * (ii - 1) + jj, n_hub*n_neighbor_hub + ii] <- 1
  }
}
diag(A_truth_mvp6hubESEN) <- 0

#   10 times
res_allen_mvp6hubESEN.list <- list()
res_my_mvp6hubESEN.list <- list()
for(cc in 1:1){
  cat("*********** sim mvp6hubESEN, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6hubESEN <- A_truth_mvp6hubESEN * 0.3
  diag(corr_mvp6hubESEN) <- 1

  df_sim_mvp6hubESEN <- data.frame(genPoisNor(n_sample, 0, n_hub*(n_neighbor_hub + 1), corr_mvp6hubESEN,
                                                rep(5, n_hub*(n_neighbor_hub + 1)),
                                                #c(rep(1:n_hub, each = n_neighbor_hub), 1:n_hub),
                                                c(), c()))
  for(ii in 1:ncol(df_sim_mvp6hubESEN)){
    df_sim_mvp6hubESEN[, ii] <- df_sim_mvp6hubESEN[, ii] + rpois(n_sample, 5)
  }


  res_allen_mvp6hubESEN.list[[cc]] <- allen_method(df_sim_mvp6hubESEN, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvp6hubESEN.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvp6hubESEN, glm_family = "poisson",
                                                                             n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                             prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvp6hubESEN.list[[1]]$summary$n_edge, res_allen_mvp6hubESEN.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvp6hubESEN.list[[1]]$summary$n_edge, res_my_mvp6hubESEN.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_mvp6hubESEN), as.vector(prob_predicting(res_allen_mvp6hubESEN.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_mvp6hubESEN), as.vector(prob_predicting(res_my_mvp6hubESEN.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvp6hubESEN distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvp6hubESEN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvp6hubESEN.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvp6hubESEN.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvp6hubESEN.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvp6hubESEN)){
    for(kk in 1:nrow(A_truth_mvp6hubESEN)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvp6hubESEN[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvp6hubESEN.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvp6hubESEN.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)

}
df_accu_mvp6hubESEN <- df_accu
df_accu_mvp6hubESEN
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvp6hubESEN, "res_sim\\df_accu_6hubLogit.csv")
write.csv(A_truth_mvp6hubESEN, "res_sim\\adj_mat_truth_6hubLogit.csv")
write.csv(res_my_mvp6hubESEN.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6hubLogit.csv")
write.csv(res_allen_mvp6hubESEN.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6hubLogit.csv")



#---- Simulation: random link (mvp) ----

n_sample <- 200
p <- 60
perc_link <- 0.1
A_truth_mvpranlink <- matrix(as.integer(runif(p * p) < perc_link), nrow = p)
A_truth_mvpranlink <- A_truth_mvpranlink + t(A_truth_mvpranlink)
A_truth_mvpranlink[A_truth_mvpranlink != 0] <- 1
diag(A_truth_mvpranlink) <- 0

#   10 times
res_allen_mvpranlink.list <- list()
res_my_mvpranlink.list <- list()
for(cc in c(1:1)){
  cat("*********** sim mvpranlink, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvpranlink <- A_truth_mvpranlink * 0.1
  diag(corr_mvpranlink) <- 1

  df_sim_mvpranlink <- data.frame(genPoisNor(n_sample, 0, p, corr_mvpranlink, ceiling(1:p / 10), c(), c()))

  for(ii in 1:p){
    df_sim_mvpranlink[, ii] <- df_sim_mvpranlink[, ii] + rpois(n_sample, ceiling(ii / 10))
  }

  res_allen_mvpranlink.list[[cc]] <- allen_method(df_sim_mvpranlink, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_mvpranlink.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvpranlink, glm_family = "poisson",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvpranlink.list[[1]]$summary$n_edge, res_allen_mvpranlink.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvpranlink.list[[1]]$summary$n_edge, res_my_mvpranlink.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_mvpranlink), as.vector(prob_predicting(res_allen_mvpranlink.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_mvpranlink), as.vector(prob_predicting(res_my_mvpranlink.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvpranlink distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvpranlink.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvpranlink.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  for(jj in 1:nrow(A_truth_mvpranlink)){
    for(kk in 1:nrow(A_truth_mvpranlink)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvpranlink[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvpranlink.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvpranlink.list[[ii]]$adj_mat[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(as.vector(A_truth_mvpranlink), as.vector(prob_predicting(res_my_mvpranlink.list[[ii]]$adj_mat.list)))
  df_accu[ii, 15] <- auc(as.vector(A_truth_mvpranlink), as.vector(prob_predicting(res_allen_mvpranlink.list[[ii]]$adj_mat.list)))
}
df_accu_mvpranlink <- df_accu
df_accu_mvpranlink
summary(df_accu)
write.csv(df_accu_mvpranlink, "res_sim\\df_accu_ranlinkLogit.csv")
write.csv(A_truth_mvpranlink, "res_sim\\adj_mat_truth_ranlinkLogit.csv")
write.csv(res_my_mvpranlink.list[[1]]$adj_mat, "res_sim\\adj_mat_my_ranlinkLogit.csv")
write.csv(res_allen_mvpranlink.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_ranlinkLogit.csv")


#---- Simulation: scale free ----

n_sample <- 200
p <- 60
set.seed(434)
A_truth_scalefree <- as.matrix(as_adjacency_matrix(sample_pa(p, power = 0.2, directed = F)))
table(apply(A_truth_scalefree, 2, sum))
qqplot(qbeta(ppoints(60), 1, 3), apply(A_truth_scalefree, 2, sum))


#   10 times
res_allen_scalefree.list <- list()
res_my_scalefree.list <- list()
for(cc in c(1:1)){
  cat("*********** sim scalefree, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_scalefree <- A_truth_scalefree * 0.1
  diag(corr_scalefree) <- 1

  df_sim_scalefree <- data.frame(genPoisNor(n_sample, 0, p, corr_scalefree, rep(5, p), c(), c()))

  for(ii in 1:p){
    df_sim_scalefree[, ii] <- df_sim_scalefree[, ii] + rpois(n_sample, 1)
  }

  res_allen_scalefree.list[[cc]] <- allen_method(df_sim_scalefree, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_scalefree.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_scalefree, glm_family = "poisson",
                                                                          n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                          prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_scalefree.list[[1]]$summary$n_edge, res_allen_scalefree.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_scalefree.list[[1]]$summary$n_edge, res_my_scalefree.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_scalefree), as.vector(prob_predicting(res_allen_scalefree.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_scalefree), as.vector(prob_predicting(res_my_scalefree.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.scalefree distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_scalefree.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_scalefree.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_scalefree.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_scalefree.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_scalefree)){
    for(kk in 1:nrow(A_truth_scalefree)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_scalefree[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_scalefree.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_scalefree.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_scalefree <- df_accu
df_accu_scalefree
summary(df_accu)
write.csv(df_accu_scalefree, "res_sim\\df_accu_mvpranlinkLogit.csv")
write.csv(A_truth_scalefree, "res_sim\\adj_mat_truth_mvpranlinkLogit.csv")
write.csv(res_my_scalefree.list[[1]]$adj_mat, "res_sim\\adj_mat_my_mvpranlinkLogit.csv")
write.csv(res_allen_scalefree.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_mvpranlinkLogit.csv")


#---- Simulation: Consistency of Adaptive Searching (hub 3*3) ----

n_neighbor_hub_small <- 3
n_hub_small <- 3
n_sample <- 200
prior.b <- 3
A_truth_mvp6hubeqhigh_small <- matrix(rep(0, (n_hub_small * (n_neighbor_hub_small + 1)) ^ 2), nrow = n_hub_small * (n_neighbor_hub_small + 1))
for(ii in c(1:n_hub_small)){
  for(jj in c(1:n_neighbor_hub_small)){
    A_truth_mvp6hubeqhigh_small[n_hub_small*n_neighbor_hub_small + ii, n_neighbor_hub_small * (ii - 1) + jj] <- 1
    A_truth_mvp6hubeqhigh_small[n_neighbor_hub_small * (ii - 1) + jj, n_hub_small*n_neighbor_hub_small + ii] <- 1
  }
}
diag(A_truth_mvp6hubeqhigh_small) <- 0

#   10 times
edges_brute <- c()
edges_adaptive <- c()
for(cc in c(1:10)){
  cat("*********** sim Consistency of Adaptive Searching (hub 4*3), epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvp6hubeqhigh_small <- A_truth_mvp6hubeqhigh_small * 0.3
  diag(corr_mvp6hubeqhigh_small) <- 1

  df_sim_mvp6hubeqhigh_small <- data.frame(genPoisNor(n_sample, 0, (n_neighbor_hub_small + 1) * n_hub_small, corr_mvp6hubeqhigh_small,
                                                #rep(1, (n_neighbor_hub_small + 1) * n_hub_small),
                                                c(rep(1:n_hub_small, each = n_neighbor_hub_small), 1:n_hub_small),
                                                c(), c()))
  for(ii in 1:ncol(df_sim_mvp6hubeqhigh_small)){
    df_sim_mvp6hubeqhigh_small[, ii] <- df_sim_mvp6hubeqhigh_small[, ii] + rpois(n_sample, 6)
  }

  for(ii in 1:ncol(df_sim_mvp6hubeqhigh_small)){
    tmp_X <- data.matrix(df_sim_mvp6hubeqhigh_small[, -ii])
    tmp_y <- data.matrix(df_sim_mvp6hubeqhigh_small[, ii])
    #best_model_given_prior <- bruteforce_model_given_prior(tmp_X, tmp_y, prior.a = 1, prior.b = 10, glm_family = "poisson")
    res_brute <- bruteforce_model_given_prior(
      X = tmp_X, y = tmp_y,
      prior.a = 1, prior.b = prior.b, glm_family = "poisson")
    edges_brute <- c(edges_brute, as.integer(res_brute$edges))

    res <- best_model_given_prior(
      X = tmp_X, y = tmp_y, glm_family = "poisson",
      n_lambda = 300, prior.a = 1, prior.b = prior.b)
    edges_adaptive <- c(edges_adaptive, res$edges)
    print(c(res_brute$p.connection, res$p.connection, mean(res_brute$edges == res$edges), res$lambda_final))
  }

}

mean(edges_brute != edges_adaptive) * 100


#---- Simulation: reproduce Allen's 3-hub network ----

n_sample <- 200
p <- 50

#   10 times
res_allen_allen3hub.list <- list()
res_my_allen3hub.list <- list()
for(cc in 1:1){
  cat("*********** sim allen3hub, epoch", cc, "**************\n")
  set.seed(434 + cc)

  tmp <- XMRF.Sim(n = n_sample, p = p, model = "LPGM", graph.type = "hub")
  df_sim_allen3hub <- t(tmp$X)
  A_truth_allen3hub <- tmp$B

  res_allen_allen3hub.list[[cc]] <- allen_method(df_sim_allen3hub, glm_family = "poisson", B = 100, n_lambda = 300)
  res_my_allen3hub.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_allen3hub, glm_family = "poisson",
                                                                           n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                           prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_allen3hub.list[[1]]$summary$n_edge, res_allen_allen3hub.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_allen3hub.list[[1]]$summary$n_edge, res_my_allen3hub.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

plot(roc(as.vector(A_truth_allen3hub), as.vector(prob_predicting(res_allen_allen3hub.list[[1]]$adj_mat.list))))
plot(roc(as.vector(A_truth_allen3hub), as.vector(prob_predicting(res_my_allen3hub.list[[1]]$adj_mat.list))))

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.allen3hub distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_allen3hub.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_allen3hub.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_allen3hub.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_allen3hub.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_allen3hub)){
    for(kk in 1:nrow(A_truth_allen3hub)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_allen3hub[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_allen3hub.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_allen3hub.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)

}
df_accu_allen3hub <- df_accu
df_accu_allen3hub
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_allen3hub, "res_sim\\df_accu_6hubLogit.csv")
write.csv(A_truth_allen3hub, "res_sim\\adj_mat_truth_6hubLogit.csv")
write.csv(res_my_allen3hub.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6hubLogit.csv")
write.csv(res_allen_allen3hub.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6hubLogit.csv")


#---- plot for paper ----

df_accu.list <- list(clu_ESEN = df_accu_mvp6clueqhigh,
                     clu_USUN = df_accu_mvp6clu,
                     clu_ESUN = df_accu_mvp6cluESUN,
                     clu_USEN = df_accu_mvp6cluUSEN,
                     hub_ESEN = df_accu_mvp6hubESEN,
                     hub_USUN = df_accu_mvp6hub,
                     hub_ESUN = df_accu_mvp6hubESUN,
                     hub_USEN = df_accu_mvp6hubeqhigh)

df_accu_cb <- data.frame()
for(key in names(df_accu.list)){
  df_accu <- df_accu.list[[key]]
  df_accu$type <- key
  df_accu_cb <- rbind(df_accu_cb, df_accu)
}
write.csv(df_accu_cb, "df_accu_cb.csv")


#---- Simulation: 6 clu (MVN) unequal noise unequal signal ----

n_node_perclu <- 6
n_clu <- 10
n_sample <- 200
A_truth_mvn6cluUNUS <- matrix(rep(0, (n_node_perclu * n_clu) ^ 2), nrow = n_node_perclu * n_clu)
for(ii in 1:n_clu){
  for(jj in 1:n_node_perclu){
    for(kk in 1:n_node_perclu){
      A_truth_mvn6cluUNUS[(ii - 1) * n_node_perclu + jj, (ii - 1) * n_node_perclu + kk] <- 1
      A_truth_mvn6cluUNUS[(ii - 1) * n_node_perclu + kk, (ii - 1) * n_node_perclu + jj] <- 1
    }
  }
}
diag(A_truth_mvn6cluUNUS) <- 0

#   10 times
res_allen_mvn6cluUNUS.list <- list()
res_my_mvn6cluUNUS.list <- list()
for(cc in 1:10){
  cat("*********** sim mvn6cluUNUS, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvn6cluUNUS <- A_truth_mvn6cluUNUS * 0.2
  diag(corr_mvn6cluUNUS) <- 1
  mu_mvn6cluUNUS <- rep(1:n_clu, each = n_node_perclu)
  for(ii in 1:nrow(corr_mvn6cluUNUS)){
    for(jj in 1:nrow(corr_mvn6cluUNUS)){
      corr_mvn6cluUNUS[ii, jj] <- corr_mvn6cluUNUS[ii, jj] * sqrt(mu_mvn6cluUNUS[ii] * mu_mvn6cluUNUS[jj])
    }
  }

  df_sim_mvn6cluUNUS <- data.frame(mvrnorm(n_sample, mu_mvn6cluUNUS, corr_mvn6cluUNUS))

  for(ii in 1:(n_clu * n_node_perclu)){
    df_sim_mvn6cluUNUS[, ii] <- df_sim_mvn6cluUNUS[, ii] + rnorm(n_sample, ceiling(ii / n_node_perclu), sqrt(ceiling(ii / n_node_perclu)))
  }

  res_allen_mvn6cluUNUS.list[[cc]] <- allen_method(df_sim_mvn6cluUNUS, glm_family = "gaussian", B = 100, n_lambda = 300)
  res_my_mvn6cluUNUS.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvn6cluUNUS, glm_family = "gaussian",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvn6cluUNUS.list[[1]]$summary$n_edge, res_allen_mvn6cluUNUS.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvn6cluUNUS.list[[1]]$summary$n_edge, res_my_mvn6cluUNUS.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvn6cluUNUS distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvn6cluUNUS.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvn6cluUNUS.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvn6cluUNUS.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvn6cluUNUS.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvn6cluUNUS)){
    for(kk in 1:nrow(A_truth_mvn6cluUNUS)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvn6cluUNUS[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvn6cluUNUS.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvn6cluUNUS.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_mvn6cluUNUS <- df_accu
df_accu_mvn6cluUNUS
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvn6cluUNUS, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvn6cluUNUS, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvn6cluUNUS.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvn6cluUNUS.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")

#---- Simulation: 6 clu (MVN) equal noise equal signal ----

n_node_perclu <- 6
n_clu <- 10
n_sample <- 200
A_truth_mvn6cluENES <- matrix(rep(0, (n_node_perclu * n_clu) ^ 2), nrow = n_node_perclu * n_clu)
for(ii in 1:n_clu){
  for(jj in 1:n_node_perclu){
    for(kk in 1:n_node_perclu){
      A_truth_mvn6cluENES[(ii - 1) * n_node_perclu + jj, (ii - 1) * n_node_perclu + kk] <- 1
      A_truth_mvn6cluENES[(ii - 1) * n_node_perclu + kk, (ii - 1) * n_node_perclu + jj] <- 1
    }
  }
}
diag(A_truth_mvn6cluENES) <- 0

#   10 times
res_allen_mvn6cluENES.list <- list()
res_my_mvn6cluENES.list <- list()
for(cc in 1:10){
  cat("*********** sim mvn6cluENES, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvn6cluENES <- A_truth_mvn6cluENES * 0.2
  diag(corr_mvn6cluENES) <- 1
  mu_mvn6cluENES <- rep(5, n_node_perclu * n_clu)
  for(ii in 1:nrow(corr_mvn6cluENES)){
    for(jj in 1:nrow(corr_mvn6cluENES)){
      corr_mvn6cluENES[ii, jj] <- corr_mvn6cluENES[ii, jj] * 5
    }
  }

  df_sim_mvn6cluENES <- data.frame(mvrnorm(n_sample, mu_mvn6cluENES, corr_mvn6cluENES))

  for(ii in 1:(n_clu * n_node_perclu)){
    df_sim_mvn6cluENES[, ii] <- df_sim_mvn6cluENES[, ii] + rnorm(n_sample, 5, sqrt(5))
  }

  res_allen_mvn6cluENES.list[[cc]] <- allen_method(df_sim_mvn6cluENES, glm_family = "gaussian", B = 100, n_lambda = 300)
  res_my_mvn6cluENES.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvn6cluENES, glm_family = "gaussian",
                                                                       n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                       prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvn6cluENES.list[[1]]$summary$n_edge, res_allen_mvn6cluENES.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvn6cluENES.list[[1]]$summary$n_edge, res_my_mvn6cluENES.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvn6cluENES distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvn6cluENES.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvn6cluENES.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvn6cluENES.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvn6cluENES.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvn6cluENES)){
    for(kk in 1:nrow(A_truth_mvn6cluENES)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvn6cluENES[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvn6cluENES.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvn6cluENES.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_mvn6cluENES <- df_accu
df_accu_mvn6cluENES
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvn6cluENES, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvn6cluENES, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvn6cluENES.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvn6cluENES.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")

#---- Simulation: 6 hub (MVN) unequal noise unequal signal ----


n_neighbor_hub <- 9
n_hub <- 6
n_sample <- 200
A_truth_mvn6hubUNUS <- matrix(rep(0, (n_hub * (n_neighbor_hub + 1)) ^ 2), nrow = n_hub * (n_neighbor_hub + 1))
for(ii in c(1:n_hub)){
  for(jj in c(1:n_neighbor_hub)){
    A_truth_mvn6hubUNUS[n_hub*n_neighbor_hub + ii, n_neighbor_hub * (ii - 1) + jj] <- 1
    A_truth_mvn6hubUNUS[n_neighbor_hub * (ii - 1) + jj, n_hub*n_neighbor_hub + ii] <- 1
  }
}
diag(A_truth_mvn6hubUNUS) <- 0

#   10 times
res_allen_mvn6hubUNUS.list <- list()
res_my_mvn6hubUNUS.list <- list()
for(cc in 1:10){
  cat("*********** sim mvn6hubUNUS, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvn6hubUNUS <- A_truth_mvn6hubUNUS * 0.2
  diag(corr_mvn6hubUNUS) <- 1
  mu_mvn6hubUNUS <- c(rep(1:n_hub, each = n_neighbor_hub), 1:n_hub)
  for(ii in 1:nrow(corr_mvn6hubUNUS)){
    for(jj in 1:nrow(corr_mvn6hubUNUS)){
      corr_mvn6hubUNUS[ii, jj] <- corr_mvn6hubUNUS[ii, jj] * sqrt(mu_mvn6hubUNUS[ii] * mu_mvn6hubUNUS[jj])
    }
  }

  df_sim_mvn6hubUNUS <- data.frame(mvrnorm(n_sample, mu_mvn6hubUNUS, corr_mvn6hubUNUS))

  for(ii in c(1:n_hub)){
    for(jj in c(1:n_neighbor_hub)){
      df_sim_mvn6hubUNUS[, n_neighbor_hub * (ii - 1) + jj] <- df_sim_mvn6hubUNUS[, n_neighbor_hub * (ii - 1) + jj] + rnorm(n_sample, ii, sqrt(ii))
    }
  }
  for(ii in c(1:n_hub)){
    df_sim_mvn6hubUNUS[, n_hub * n_neighbor_hub + ii] <- df_sim_mvn6hubUNUS[, n_hub * n_neighbor_hub + ii] + rnorm(n_sample, ii, sqrt(ii))
  }

  res_allen_mvn6hubUNUS.list[[cc]] <- allen_method(df_sim_mvn6hubUNUS, glm_family = "gaussian", B = 100, n_lambda = 300)
  res_my_mvn6hubUNUS.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvn6hubUNUS, glm_family = "gaussian",
                                                                           n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                           prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvn6hubUNUS.list[[1]]$summary$n_edge, res_allen_mvn6hubUNUS.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvn6hubUNUS.list[[1]]$summary$n_edge, res_my_mvn6hubUNUS.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvn6hubUNUS distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvn6hubUNUS.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvn6hubUNUS.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvn6hubUNUS.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvn6hubUNUS.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvn6hubUNUS)){
    for(kk in 1:nrow(A_truth_mvn6hubUNUS)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvn6hubUNUS[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvn6hubUNUS.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvn6hubUNUS.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_mvn6hubUNUS <- df_accu
df_accu_mvn6hubUNUS
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvn6hubUNUS, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvn6hubUNUS, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvn6hubUNUS.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvn6hubUNUS.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")


#---- Simulation: 6 hub (MVN) unequal noise unequal signal ----


n_neighbor_hub <- 9
n_hub <- 6
n_sample <- 200
A_truth_mvn6hubENES <- matrix(rep(0, (n_hub * (n_neighbor_hub + 1)) ^ 2), nrow = n_hub * (n_neighbor_hub + 1))
for(ii in c(1:n_hub)){
  for(jj in c(1:n_neighbor_hub)){
    A_truth_mvn6hubENES[n_hub*n_neighbor_hub + ii, n_neighbor_hub * (ii - 1) + jj] <- 1
    A_truth_mvn6hubENES[n_neighbor_hub * (ii - 1) + jj, n_hub*n_neighbor_hub + ii] <- 1
  }
}
diag(A_truth_mvn6hubENES) <- 0

#   10 times
res_allen_mvn6hubENES.list <- list()
res_my_mvn6hubENES.list <- list()
for(cc in 1:10){
  cat("*********** sim mvn6hubENES, epoch", cc, "**************\n")
  set.seed(434 + cc)

  corr_mvn6hubENES <- A_truth_mvn6hubENES * 0.2
  diag(corr_mvn6hubENES) <- 1
  mu_mvn6hubENES <- rep(5, n_hub * (n_neighbor_hub + 1))
  for(ii in 1:nrow(corr_mvn6hubENES)){
    for(jj in 1:nrow(corr_mvn6hubENES)){
      corr_mvn6hubENES[ii, jj] <- corr_mvn6hubENES[ii, jj] * sqrt(mu_mvn6hubENES[ii] * mu_mvn6hubENES[jj])
    }
  }

  df_sim_mvn6hubENES <- data.frame(mvrnorm(n_sample, mu_mvn6hubENES, corr_mvn6hubENES))

  for(ii in c(1:n_hub)){
    for(jj in c(1:n_neighbor_hub)){
      df_sim_mvn6hubENES[, n_neighbor_hub * (ii - 1) + jj] <- df_sim_mvn6hubENES[, n_neighbor_hub * (ii - 1) + jj] + rnorm(n_sample, 5, sqrt(5))
    }
  }
  for(ii in c(1:n_hub)){
    df_sim_mvn6hubENES[, n_hub * n_neighbor_hub + ii] <- df_sim_mvn6hubENES[, n_hub * n_neighbor_hub + ii] + rnorm(n_sample, 5, sqrt(5))
  }

  res_allen_mvn6hubENES.list[[cc]] <- allen_method(df_sim_mvn6hubENES, glm_family = "gaussian", B = 100, n_lambda = 300)
  res_my_mvn6hubENES.list[[cc]] <- my_prior_stability_search_samesubsample(df_sim_mvn6hubENES, glm_family = "gaussian",
                                                                           n_lambda = 100, B = 100, prior_a.vec = 1,
                                                                           prior_b.vec = c(2:30, 35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))
}

plot(res_allen_mvn6hubENES.list[[1]]$summary$n_edge, res_allen_mvn6hubENES.list[[1]]$summary$instability,
     xlim = c(0, 500), ylim = c(0, 0.1),
     main = "Instability Comparison (Allen's vs Ours)",
     xlab = "number of edges in total graph",
     ylab = "instability (2 * average variance of all slots of edges)")
points(res_my_mvn6hubENES.list[[1]]$summary$n_edge, res_my_mvn6hubENES.list[[1]]$summary$instability, col = "red")
abline(h = 0.05)

df_accu <- data.frame(matrix(nrow = 0, ncol = 15))
colnames(df_accu) <- c("epoch", "err_rate_my", "err_rate_allen",
                       "weighted_err_my", "weighted_err_allen",
                       "sensitivity_my", "sensitivity_allen",
                       "specificity_my", "specificity_allen",
                       "false_positive_my", "false_positive_allen",
                       "false_negative_my", "false_negative_allen",
                       "AUC_my", "AUC_allen")
plot(NULL, xlim = c(0, 17), ylim = c(0, 25), xlab = "degree", ylab = "Freq",
     main = "sim.mvn6hubENES distribution comparison LHGM vs LLGM")
for(ii in 1:10){
  tmp_df <- data.frame(table(apply(res_my_mvn6hubENES.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "red", type = "l")
  tmp_df <- data.frame(table(apply(res_allen_mvn6hubENES.list[[ii]]$adj_mat, 2, sum)))
  tmp_df[, 1] <- as.numeric(as.character(tmp_df[, 1]))
  for(jj in 1:19){
    if(!(jj %in% tmp_df[,1])){
      tmp_df[nrow(tmp_df) + 1, ]<- c(jj, 0)
    }
  }
  tmp_df <- tmp_df[order(tmp_df[, 1]), ]
  points(tmp_df[, 1], tmp_df[, 2] + rnorm(nrow(tmp_df), 0, 0.1), col = "blue", type = "l")

  y_true <- c()
  y_pred_my <- c()
  y_pred_allen <- c()
  y_pred_prob_my <- c()
  y_pred_prob_allen <- c()
  mat_prob_pred_my <- prob_predicting(res_my_mvn6hubENES.list[[ii]]$adj_mat.list)
  mat_prob_pred_allen <- prob_predicting(res_allen_mvn6hubENES.list[[ii]]$adj_mat.list)
  for(jj in 1:nrow(A_truth_mvn6hubENES)){
    for(kk in 1:nrow(A_truth_mvn6hubENES)){
      if(jj <= kk){
        next
      }
      y_true <- c(y_true, A_truth_mvn6hubENES[jj, kk])
      y_pred_my <- c(y_pred_my, res_my_mvn6hubENES.list[[ii]]$adj_mat[jj, kk])
      y_pred_allen <- c(y_pred_allen, res_allen_mvn6hubENES.list[[ii]]$adj_mat[jj, kk])
      y_pred_prob_my <- c(y_pred_prob_my, mat_prob_pred_my[jj, kk])
      y_pred_prob_allen <- c(y_pred_prob_allen, mat_prob_pred_allen[jj, kk])
    }
  }

  cat("epoch", ii, "\n")
  confmat_my <- table(y_true, y_pred_my)
  confmat_allen <- table(y_true, y_pred_allen)
  print(confmat_my)
  print(confmat_allen)
  cat("err_rate_my =", (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my), " ")
  cat("false_positive_my =", confmat_my[1, 2] / sum(confmat_my[1, ]), " ")
  cat("false_negative_my =", confmat_my[2, 1] / sum(confmat_my[2, ]), "\n")
  cat("err_rate_allen =", (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen), " ")
  cat("false_positive_allen =", confmat_allen[1, 2] / sum(confmat_allen[1, ]), " ")
  cat("false_negative_allen =", confmat_allen[2, 1] / sum(confmat_allen[2, ]), "\n")
  df_accu[ii, 1] <- ii
  df_accu[ii, 2] <- (confmat_my[1, 2] + confmat_my[2, 1]) / sum(confmat_my)
  df_accu[ii, 3] <- (confmat_allen[1, 2] + confmat_allen[2, 1]) / sum(confmat_allen)
  df_accu[ii, 10] <- confmat_my[1, 2] / sum(confmat_my[1, ])

  df_accu[ii, 6] <- confmat_my[2, 2] / sum(confmat_my[2, ])
  df_accu[ii, 8] <- confmat_my[2, 2] / sum(confmat_my[, 2])
  df_accu[ii, 7] <- confmat_allen[2, 2] / sum(confmat_allen[2, ])
  df_accu[ii, 9] <- confmat_allen[2, 2] / sum(confmat_allen[, 2])

  df_accu[ii, 11] <- confmat_allen[1, 2] / sum(confmat_allen[1, ])
  df_accu[ii, 12] <- confmat_my[2, 1] / sum(confmat_my[2, ])
  df_accu[ii, 13] <- confmat_allen[2, 1] / sum(confmat_allen[2, ])
  df_accu[ii, 4] <- (df_accu[ii, 10] + df_accu[ii, 12]) / 2
  df_accu[ii, 5] <- (df_accu[ii, 11] + df_accu[ii, 13]) / 2

  df_accu[ii, 14] <- auc(y_true, y_pred_prob_my)
  df_accu[ii, 15] <- auc(y_true, y_pred_prob_allen)
}
df_accu_mvn6hubENES <- df_accu
df_accu_mvn6hubENES
print(xtable(df_accu[, c("epoch", "weighted_err_my", "weighted_err_allen", "sensitivity_my",
                         "sensitivity_allen", "specificity_my", "specificity_allen", "AUC_my", "AUC_allen")],
             type = "latex"), file = "tmp.tex")
summary(df_accu)
write.csv(df_accu_mvn6hubENES, "res_sim\\df_accu_6cluLogit.csv")
write.csv(A_truth_mvn6hubENES, "res_sim\\adj_mat_truth_6cluLogit.csv")
write.csv(res_my_mvn6hubENES.list[[1]]$adj_mat, "res_sim\\adj_mat_my_6cluLogit.csv")
write.csv(res_allen_mvn6hubENES.list[[1]]$adj_mat, "res_sim\\adj_mat_allen_6cluLogit.csv")


#---- Finished ----




