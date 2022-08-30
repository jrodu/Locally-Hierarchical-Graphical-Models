

library(glmnet)
library(glm2)
library(MASS)


lhgm <- function(df, prior_a.vec, prior_b.vec, glm_family = "binomial", theta = NULL,
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
    cat("subsample", b, "/", B, "finished", "\n")
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
      fit <- best_model_given_prior(X = data.matrix(df[, -i]), y = data.matrix(df[, i]), prior_weight = prior_weight,
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
      #print("lambda vec:")
      #print(lambda_final.vec)
      cat("prior.a = ", prior.a, ", prior.b = ", prior.b, ", instability = ", instability, "\n")
    }
  }
  
  df_summary <- data.frame(prior_a = prior_a.vec, prior_b = prior_b.vec,
                           instability = instability.vec,
                           n_edge = n_edge.vec)
  
  tmp <- df_summary
  tmp$index <- c(1: nrow(tmp))
  tmp <- tmp[tmp$instability <= 0.1, ]
  if(is_sparse_graph){
    tmp <- tmp[tmp$n_edge < p * (p - 1) / 4, ]
  }
  tmp <- tmp[order(tmp$n_edge, decreasing = T), ]
  tmp_index <- tmp[1, "index"]
  prior_a.opt <- tmp[1, "prior_a"]
  prior_b.opt <- tmp[1, "prior_b"]
  adj_mat <- adj_mat.list[[tmp_index]]
  
  res <- list(summary = df_summary, 
              adj_mat.list = adj_mat.list,
              adj_mat = adj_mat,
              prior_a.opt = prior_a.opt,
              prior_b.opt = prior_b.opt)
  return(res)
}


instability_measure <- function(L_adj_mat.array) {
  mean_adj_mat <- apply(L_adj_mat.array, c(1, 2), mean)
  p <- nrow(mean_adj_mat)
  tot_instability <- sum(mean_adj_mat * (1 - mean_adj_mat)) / (p * (p - 1)) * 2
  return(tot_instability)
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


#---- Test: Real-word data TCGA UECE (DELETE before packaging) ----

df_tcgauece <- read.csv("endometrial_cancer.csv", header = T,sep = "\t")
df_tcgauece <- df_tcgauece[, -1]
df_tcgauecePois <- df_tcgauece

for(ii in 1:ncol(df_tcgauece)){
  df_tcgauece[, ii] <- as.integer(df_tcgauece[, ii] != 0)
}

res_my_tcgauece <- lhgm(df_tcgauece,  n_lambda = 100, B = 100, prior_a.vec = 1,
                        prior_b.vec = c(5, 10:35, 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))

df_tcgauecePois_transformed <- poisson_transform(df_tcgauecePois)

res_my_tcgauecePois <- lhgm(df_tcgauecePois_transformed, glm_family = "poisson", n_lambda = 100, B = 100,
                            prior_a.vec = 1, prior_b.vec = c(5, seq(10, 35, 0.5), 40, 45, 50, 55, 60, 65, 70, 80, 90, 100))


