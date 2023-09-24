
get_uncertainty_info <- function(data_info, linear_model) {
  # extract info
  beta <- linear_model$beta_soln
  names(beta) <- linear_model$cov_names
  beta <- beta[["signal"]]
  gamma <- linear_model$gamma_soln[[1]]
  beta_sd <- get_beta_sd(linear_model)
  gamma_sd <- get_gamma_sd(linear_model)
  
  # compute uncertainty of fixed effects
  inner_fe_sd <- sqrt(beta_sd^2 + gamma)
  outer_fe_sd <- sqrt(beta_sd^2 + gamma + 2*gamma_sd)
  inner_betas <- c(qnorm(0.05, mean=beta, sd=inner_fe_sd),
                   qnorm(0.95, mean=beta, sd=inner_fe_sd))
  outer_betas <- c(qnorm(0.05, mean=beta, sd=outer_fe_sd),
                   qnorm(0.95, mean=beta, sd=outer_fe_sd))
  
  df <- data_info$df
  
  # compute the lower and upper draws
  pred_lrr <- data_info$pred_lrr
  inner_draws <- rbind(inner_betas[1]*pred_lrr, inner_betas[2]*pred_lrr) 
  outer_draws <- rbind(outer_betas[1]*pred_lrr, outer_betas[2]*pred_lrr)
  # compute and score
  pred_exp <- data_info$pred_exp
  exp_bounds <- data_info$exp_bounds
  index <- (pred_exp >= exp_bounds[1]) & (pred_exp <= exp_bounds[2])
  sign_score <- sign(mean(pred_lrr[index]))
  score <- min(rowMeans(outer_draws[,index])*sign_score)
  list(
    score = score,
    inner_draws = inner_draws,
    outer_draws = outer_draws
  )
}


plot_BPRFmodel <- function(data_info,draws,
                           uncertainty_info,
                           linear_model,
                           signal_model,
                           uncertainty_info_fill = NULL,
                           linear_model_fill = NULL) {
  evidence_BPRF<-paste(ro_pair,"evidence_BPRF.tiff")
  tiff(filename = evidence_BPRF,
       width =1000, height = 850, units = "px", pointsize = 12,res=100)
  #main <- paste0(ro_pair, ": score=", "0.095")
  main <- paste0(ro_pair, ": score=", round(uncertainty_info$score,3))
  pred_exp <- data_info$pred_exp
  pred_lrr <- data_info$pred_lrr
  inner_draws <- uncertainty_info$inner_draws
  outer_draws <- uncertainty_info$outer_draws
  xrange <- max(pred_exp) - min(pred_exp)
  draw_mean <- apply(draws, 1, function(x) mean(x))
  #draw_lower_fe <- apply(draws, 1, function(x) quantile(x, probs=.025))
  #draw_lower <- apply(draws, 1, function(x) quantile(x, probs=.05))
  draw_upper <- apply(draws, 1, function(x) quantile(x, probs=.95))
  # plot data
  plot(df$alt_mid, df$alt_pred,
       col = alpha('gray', 2.0), cex = 0.15/df$obs_se,
       xlim = c(min(pred_exp) - xrange/20, max(pred_exp) + xrange/20),
       ylim = c(-1.2,0.5),
       xlab = "exposure", ylab = "ln_relative_risk", pch = 19)
  df_outlier <- df[df$outlier,]
  points(df_outlier$alt_mid, df_outlier$alt_pred,
         cex = 0.15/df_outlier$obs_se, col = alpha('red',2.0), pch = 4, lwd=6)
  
  # plot prediction
  lines(pred_exp, pred_lrr, col="#008080", lwd=2.5)
  lines(draws$exposure, draw_upper,col='#AA0C2D', lwd=2.5)
  # plot uncertainties
  polygon(c(pred_exp, rev(pred_exp)), c(inner_draws[1,], rev(inner_draws[2,])),
          col=alpha("#69b3a2", 0.2), border=FALSE)
  polygon(c(pred_exp, rev(pred_exp)), c(outer_draws[1,], rev(outer_draws[2,])),
          col=alpha("#69b3a2", 0.3), border=FALSE)
  
  # plot filled model
  if (!is.null(linear_model_fill)){
    lines(pred_exp, uncertainty_info_fill$outer_draws[1,], lty='dashed', col=alpha('gray', 0.9))
    lines(pred_exp, uncertainty_info_fill$outer_draws[2,], lty='dashed', col=alpha('gray', 0.9))
    df_fill <- df[grepl('fill', df$study_id, fixed=TRUE),]
    points(df_fill$alt_mid, df_fill$alt_pred,
           cex = 0.1/df_fill$obs_se, pch=1, col = alpha("#008080", 0.8))
    #main <- paste0(main, ", score_fill= 0.095" )
    main <- paste0(main, ", score_fill= ", round(uncertainty_info_fill$score, 3))
  }
  # plot bounds
  abline(v = data_info$exp_bounds[1], lty = 'dashed', lwd = 1, col = 'black')
  abline(v = data_info$exp_bounds[2], lty = 'dashed', lwd = 1, col = 'black')
  # plot 0 line
  abline(h = 0.0, lwd = 1, col = 'black')
  title(main = main) 
  dev.off()
}

plot_BPRFmodel_rr <- function(data_info,draws,
                              uncertainty_info,
                              linear_model,
                              signal_model,
                              uncertainty_info_fill = NULL,
                              linear_model_fill = NULL) {
  evidence_BPRF<-paste(ro_pair,"evidence_BPRF_rr.tiff")
  tiff(filename = evidence_BPRF,
       width =1000, height = 850, units = "px", pointsize = 12,res=100)
  #main <- paste0(ro_pair, ": score=", "0.095")
  main <- paste0(ro_pair, ": score=", round(uncertainty_info$score,3))
  pred_exp <- data_info$pred_exp
  pred_lrr <- data_info$pred_lrr
  inner_draws <- uncertainty_info$inner_draws
  outer_draws <- uncertainty_info$outer_draws
  xrange <- max(pred_exp) - min(pred_exp)
  draw_mean <- apply(draws, 1, function(x) mean(x))
  #draw_lower_fe <- apply(draws, 1, function(x) quantile(x, probs=.025))
  #draw_lower <- apply(draws, 1, function(x) quantile(x, probs=.05))
  draw_upper <- apply(draws, 1, function(x) quantile(x, probs=.95))
  # plot data
  plot(df$alt_mid, exp(df$alt_pred),
       col = alpha('gray', 2.0), cex = 0.15/df$obs_se,
       xlim = c(min(pred_exp) - xrange/20, max(pred_exp) + xrange/20),
       ylim = c(0.0,1.3),
       xlab = "exposure", ylab = "Relative risk", pch = 19)
  df_outlier <- df[df$outlier,]
  points(df_outlier$alt_mid, exp(df_outlier$alt_pred),
         cex = 0.15/df_outlier$obs_se, col = alpha('red', 2.0), pch = 4, lwd=6)
  
  # plot prediction
  lines(pred_exp, exp(pred_lrr), col="#008080", lwd=2.5)
  lines(draws$exposure, exp(draw_upper),col='#AA0C2D', lwd=2.5)
  # plot uncertainties
  polygon(c(pred_exp, rev(pred_exp)), c(exp(inner_draws[1,]), rev(exp(inner_draws[2,]))),
          col=alpha("#69b3a2", 0.2), border=FALSE)
  polygon(c(pred_exp, rev(pred_exp)), c(exp(outer_draws[1,]), rev(exp(outer_draws[2,]))),
          col=alpha("#69b3a2", 0.3), border=FALSE)
  
  # plot filled model
  if (!is.null(linear_model_fill)){
    lines(pred_exp, exp(uncertainty_info_fill$outer_draws[1,]), lty='dashed', col=alpha('gray', 0.9))
    lines(pred_exp, exp(uncertainty_info_fill$outer_draws[2,]), lty='dashed', col=alpha('gray', 0.9))
    df_fill <- df[grepl('fill', df$study_id, fixed=TRUE),]
    points(df_fill$alt_mid, exp(df_fill$alt_pred),
           cex = 0.1/df_fill$obs_se, pch=1, col = alpha("#008080", 0.8))
    #main <- paste0(main, ", score_fill= 0.095" )
    main <- paste0(main, ", score_fill=", round(uncertainty_info_fill$score, 3))
  }
  
  # plot bounds
  abline(v = exp(data_info$exp_bounds[1]), lty = 'dashed', lwd = 1, col = 'black')
  abline(v = exp(data_info$exp_bounds[2]), lty = 'dashed', lwd = 1, col = 'black')
  # plot 0 line
  abline(h = 1.0, lwd = 1, col = 'black')
  title(main = main) 
  dev.off()
}