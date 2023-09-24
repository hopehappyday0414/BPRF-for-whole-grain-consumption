library(scales)
#source("./src/utils/egger_functions.R")

extract_data_info <- function(signal_model,
                              linear_model,
                              ref_covs,
                              alt_covs,
                              exp_quantiles = c(0.15, 0.85),
                              exp_bounds = NULL,
                              num_points = 150L,
                              normalize_to = NULL) {
  # extract dataframe
  df_signal <- signal_model$data$to_df()
  df_linear <- linear_model$data$to_df()
  df_signal$data_id <- 1:nrow(df_signal)
  common_cols <- intersect(colnames(df_signal), colnames(df_linear))
  df <- merge(df_signal, df_linear, by = common_cols, all.x = TRUE)
  df <- df[!duplicated(df$data_id),]
  df <- df[order(df$data_id),]
  
  # get exposure information
  ref_exp <- df[ref_covs]
  alt_exp <- df[alt_covs]
  df$ref_mid <- rowMeans(ref_exp)
  df$alt_mid <- rowMeans(alt_exp)
  
  exp_limits <- c(min(min(ref_exp), min(alt_exp)), max(max(ref_exp), max(alt_exp)))
  pred_exp <- seq(exp_limits[1], exp_limits[2], length = 100L)
  if (is.null(NULL)) {
    exp_bounds = c(quantile(df$ref_mid, 0.15), quantile(df$alt_mid, 0.85))
  }
  
  # create temperary prediction
  covs <- data.frame(
    c(sapply(ref_covs, function(x) rep(exp_limits[1], length.out = 100L), simplify = FALSE, USE.NAMES = TRUE),
      sapply(alt_covs, function(x) pred_exp, simplify = FALSE, USE.NAMES = TRUE))
  )
  data <- mrtool$MRData()
  data$load_df(covs, col_covs=c(ref_covs, alt_covs))
  pred <- signal_model$predict(data)
  
  # extract normalization exposure
  if (is.null(NULL)) {
    normalize_to <- exp_limits[1]
    normalize_value <- pred[1]
  } else if (normalize_to == "min") {
    normalize_to <- pred_exp[which.min(pred)]
    normalize_value <- min(pred)
  } else {
    normalize_value <- approx(pred_exp, pred, xout=normalize_to)$y
  }
  
  # create prediction of log relative risk
  pred_lrr <- pred - normalize_value
  
  # add outlier column
  w <- t(do.call(rbind, 
                 lapply(1:length(signal_model$sub_models), 
                        function(i){signal_model$sub_models[[i]]$w_soln}))
  ) %*% signal_model$weights
  df<-cbind(df,w)
  df$outlier <- ifelse(df$w<0.1, TRUE, FALSE)
  
  # add signal column
  df$signal = signal_model$predict(signal_model$data)
  
  # create prediction at reference and alternative mid points
  covs <- data.frame(
    c(sapply(ref_covs, function(x) rep(normalize_to, length.out = nrow(df)), simplify = FALSE, USE.NAMES = TRUE),
      sapply(alt_covs, function(x) df$ref_mid, simplify = FALSE, USE.NAMES = TRUE))
  )
  data <- mrtool$MRData()
  data$load_df(covs, col_covs=c(ref_covs, alt_covs))
  df$ref_pred <- signal_model$predict(data)
  df$alt_pred <- df$ref_pred + df$obs
  df$residual <- df$obs - df$signal
  df$residual_se <- sqrt(df$obs_se^2 + df$signal^2*linear_model$gamma_soln[[1]])
  
  list(
    df = df,
    ref_covs = ref_covs,
    alt_covs = alt_covs,
    exp_limits = exp_limits,
    exp_bounds = exp_bounds,
    exp_quantiles = exp_quantiles,
    pred_exp = pred_exp,
    pred_lrr = pred_lrr,
    normalize_to = normalize_to,
    normalize_value = normalize_value
  )
}


get_df_fill <- function(df) {
  # get number of filled data points
  absr_rank <- rep(0L, length=nrow(df))
  absr_rank[order(abs(df$residual))] <- seq(1, nrow(df))
  sort_index <- order(df$residual, decreasing = mean(df$residual/df$residual_se) > 0)
  num_fill <- nrow(df) - absr_rank[tail(sort_index, n=1)]
  
  # get the fill-dataframe
  df_fill <- df[sort_index[1:num_fill],]
  df_fill$study_id <- paste0('fill_', df_fill$study_id)
  df_fill$residual <- - df_fill$residual
  df_fill$obs <- df_fill$signal + df_fill$residual
  df_fill$alt_pred <- df_fill$obs + df_fill$ref_pred
  df_fill
}

plot_residual <- function(df, title){
  residualplot<-paste(ro_pair, "residual.tiff")
  tiff(filename =  residualplot,
       width =1000, height = 850, units = "px", pointsize = 12, res=100)
  residual <- df$residual
  obs_se <- df$residual_se
  max_obs_se <- quantile(obs_se, 0.99)
  fill_index <- rownames(df)[grepl('fill', df$study_id, fixed= TRUE)]
  # funnel plot
  plot(residual, obs_se, pch=19, col=alpha('gray', 1.0), 
       ylim=c(max_obs_se,0), xlim=c(-2*max_obs_se, 2*max_obs_se),
       yaxs='i', xlab="residual", ylab="ln relative risk", main=title)
  if (length(fill_index) > 0){
    points(df[fill_index,'residual'], df[fill_index,'residual_se'], 
           col=alpha('#008080',1.0), pch=19,)
  }
  
  if (sum(df$outlier) > 0) {
    outlier_index <- rownames(df)[df$outlier==TRUE]
    points(df[outlier_index, 'residual'], df[outlier_index, 'residual_se'],
           col=alpha('red', 1.8), pch=4,lwd=4)
  }
  x <- c(0, -1.96*max_obs_se, 1.96*max_obs_se)
  y <- c(0, max_obs_se, max_obs_se)
  polygon(x,y,col=alpha('#B0E0E6', 0.4))
  lines(c(0, 0), c(0, max_obs_se), lty='dashed')
  lines(c(0, -1.96*max_obs_se), c(0, max_obs_se), col='#87CEFA')
  lines(c(0, 1.96*max_obs_se), c(0, max_obs_se), col='#87CEFA')
  dev.off()
}

get_gamma_sd <- function(model){
  gamma <- model$gamma_soln
  gamma_fisher <- model$lt$get_gamma_fisher(gamma)
  return(1/sqrt(gamma_fisher[1,1]))
}

get_beta_sd <- function(model){
  model_specs <- mrtool$core$other_sampling$extract_simple_lme_specs(model)
  beta_hessian <- mrtool$core$other_sampling$extract_simple_lme_hessian(model_specs)
  beta_sd <- sqrt(diag(solve(beta_hessian)))
  names(beta_sd) <- model$cov_names
  return(beta_sd[["signal"]])
}


summarize_model <- function(data_info,
                            uncertainty_info,
                            linear_model,
                            signal_model,
                            egger_model,
                            egger_model_all,
                            uncertainty_info_fill = NULL,
                            linear_model_fill = NULL) {
  num_fill <- sum(grepl('fill', data_info$df$study_id, fixed=TRUE))
  summary <- list(
    ro_pair = data_info$ro_pair,
    num_fill = num_fill, 
    gamma = linear_model$gamma_soln[[1]], 
    gamma_sd = get_gamma_sd(linear_model),
    score = uncertainty_info$score
  )
  
  if (!is.null(linear_model_fill)){
    summary$gamma_adjusted <- linear_model_fill$gamma_soln[[1]]
    summary$gamma_sd_adjusted <- get_gamma_sd(linear_model_fill)
    summary$score_adjusted <- uncertainty_info_fill$score
  }
  
  for (name in names(egger_model)){
    summary[[paste0("egger_",name)]] <- egger_model[[name]]
  }
  
  for (name in names(egger_model_all)){
    summary[[paste0("egger_all_",name)]] <- egger_model_all[[name]]
  }
  
  data.frame(summary)
}

get_draws <- function(data_info, linear_model, num_draws = 1000L) {
  # extract info
  beta <- linear_model$beta_soln
  names(beta) <- linear_model$cov_names
  beta <- beta[["signal"]]
  gamma <- linear_model$gamma_soln[[1]]
  beta_sd <- get_beta_sd(linear_model)
  gamma_sd <- get_gamma_sd(linear_model)
  
  outer_fe_sd <- sqrt(beta_sd^2 + gamma + 2*gamma_sd)
  beta_samples <- rnorm(num_draws, mean = beta, sd = outer_fe_sd)
  
  draws <- as.data.frame(
    cbind(data_info$pred_exp, outer(data_info$pred_lrr, beta_samples))
  )
  names(draws) <- c("exposure", sapply(1:num_draws, function(i) paste0("draw_", i)))
  draws
}

egger_regression <- function(residual, residual_sd, one_sided = TRUE) {
  weighted_residual <- residual/residual_sd
  r_mean <- mean(weighted_residual)
  r_sd <- 1/sqrt(length(weighted_residual))
  r_pval <- get_pval(r_mean, r_sd, one_sided = one_sided)
  list(mean = r_mean, sd = r_sd, pval = r_pval)
}

get_pval <- function(beta, beta_sd, one_sided = FALSE) {
  zscore <- abs(beta/beta_sd)
  if (one_sided) {
    pval <- 1 - pnorm(zscore)
  } else {
    pval <- 2*(1 - pnorm(zscore))
  }
  pval
}

get_cov_names <- function(signal_model) {
  cov_model <- signal_model$sub_models[[1]]$cov_models[[1]]
  list(alt_covs = cov_model$alt_cov,
       ref_covs = cov_model$ref_cov)
}

get_risk_limits <- function(signal_model) {
  cov_names <- get_cov_names(signal_model)
  risk_data <- signal_model$data$get_covs(unlist(cov_names))
  c(min(risk_data), max(risk_data))
}

get_signal <- function(signal_model, risk) {
  cov_names <- get_cov_names(signal_model)
  risk_limits <- get_risk_limits(signal_model)
  df_covs <- data.frame(
    c(sapply(cov_names$ref_covs, function(x) rep(risk_limits[1], length.out = length(risk)),
             simplify = FALSE, USE.NAMES = TRUE),
      sapply(cov_names$alt_covs, function(x) risk,
             simplify = FALSE, USE.NAMES = TRUE))
  )
  data <- mrtool$MRData()
  data$load_df(df_covs, col_covs=unlist(cov_names))
  signal_model$predict(data)
}

get_beta <- function(linear_model) {
  beta <- linear_model$beta_soln
  names(beta) <- linear_model$cov_names
  specs <- mrtool$core$other_sampling$extract_simple_lme_specs(linear_model)
  beta_hessian <- mrtool$core$other_sampling$extract_simple_lme_hessian(specs)
  beta_sd <- 1/sqrt(diag(beta_hessian))
  names(beta_sd) <- linear_model$cov_names
  c(beta["signal"], beta_sd["signal"])
}

get_gamma <- function(linear_model) {
  gamma <- linear_model$gamma_soln[[1]]
  gamma_fisher <- linear_model$lt$get_gamma_fisher(linear_model$gamma_soln)
  gamma_sd <- 1/sqrt(diag(gamma_fisher)[[1]])
  c(gamma, gamma_sd)
}

get_soln <- function(linear_model) {
  list(
    beta_soln = get_beta(linear_model),
    gamma_soln = get_gamma(linear_model)
  )
}

get_ln_rr_draws <- function(signal_model,
                            linear_model,
                            risk,
                            num_draws = 150L,
                            normalize_to_tmrel = FALSE,
                            include_re = TRUE) {
  # set seed inside function
  set.seed(1234)
  
  signal <- get_signal(signal_model, risk)
  re_signal <- signal
  soln <- get_soln(linear_model)
  
  fe_samples <- rnorm(num_draws, mean=soln$beta[1], sd=soln$beta[2])
  re_samples <- rnorm(num_draws, mean=0, sd=sqrt(soln$gamma[1] + 2*soln$gamma[2]))
  
  draws <- outer(signal, fe_samples)
  if (include_re) {
    draws <- draws + outer(re_signal, re_samples)
  }
  
  if (normalize_to_tmrel) {
    tmrel_index <- which.min(signal)
    draws <- apply(draws, 2, function(x) x - x[tmrel_index])
  }
  
  df <- as.data.frame(cbind(risk, draws))
  names(df) <- c("risk", sapply(1:num_draws, function(i) paste0("draw_", i)))
  return(df)
}


summarize_draws <- function(data){
  
  df <- as.data.table(copy(data))
  draw_cols <- colnames(df)[grepl("draw_", colnames(df))]
  
  df[, mean := apply(.SD, 1, mean), .SDcols = draw_cols]
  df[, upper := apply(.SD, 1, quantile, 0.975), .SDcols = draw_cols]
  df[, lower := apply(.SD, 1, quantile, 0.025), .SDcols = draw_cols]
  
  df[, (draw_cols) := NULL]
  return(df)
  
}
