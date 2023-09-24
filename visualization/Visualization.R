#-------------------------------------------------------------------------------------------------------------------------------------------------
# 05_publication_bias and plots.R
#-------------------------------------------------------------------------------------------------------------------------------------------------
ro_pair<-"T2D"
library(data.table)
library(reticulate)
library(ggplot2)
np <- import("numpy")

linear_model_path <- paste0("/home/hopehappy/BPRF_wholegrain/04_mixed_effects_pkl_files/", ro_pair, "_linear_model", "_non_trimmed.pkl")
signal_model_path <- paste0("/home/hopehappy/BPRF_wholegrain/01_template_pkl_files/", ro_pair, "_signal_model", "_non_trimmed.pkl")
ref_covs <- c("a_0", "a_1")
alt_covs <- c("b_0", "b_1")

### Load model objects
linear_model <- py_load_object(filename = linear_model_path, pickle = "dill")
signal_model <- py_load_object(filename = signal_model_path, pickle = "dill")

linear_model$predict
data_info <- extract_data_info(signal_model,
                               linear_model,
                               ref_covs = ref_covs,
                               alt_covs = alt_covs)

data_info$ro_pair <- "T2D"
df <- data_info$df

### Detect publication bias
df_no_outlier <- df[!df$outlier,]
egger_model_all <- egger_regression(df$residual, df$residual_se)
egger_model <- egger_regression(df_no_outlier$residual, df_no_outlier$residual_se)
has_pub_bias <- egger_model$pval < 0.05

### Adjust for publication bias
if (has_pub_bias) {
  df_fill <- get_df_fill(df[!df$outlier,])
  num_fill <- nrow(df_fill)
} else {
  num_fill <- 0
}

# fill the data if needed and refit the model
if (num_fill > 0) {
  df <- rbind(df, df_fill)
  data_info$df <- df
  df2<-df[!df$outlier,]
  
  # refit the model
  data <- mrtool$MRData()
  data$load_df(
    data = df,
    col_obs = 'obs',
    col_obs_se = 'obs_se',
    col_study_id = 'study_id',
    col_covs = as.list(linear_model$cov_names)
  )
  linear_model_fill <- mrtool$MRBRT(data, cov_models=linear_model$cov_models)
  linear_model_fill$fit_model()
} else {
  linear_model_fill <- NULL
}

### Extract scores
uncertainty_info <- get_uncertainty_info(data_info, linear_model)
if (is.null(linear_model_fill)) {
  uncertainty_info_fill <- NULL
} else {
  uncertainty_info_fill <- get_uncertainty_info(data_info, linear_model_fill)
}

### Output diagnostics
# figures
draws <- get_draws(data_info, linear_model)
###risk is a numeric sequence
risk <- 0:400
df_rr <- get_ln_rr_draws(signal_model,
                         linear_model,
                         risk,
                         # normalize_to_tmrel = is_j_shaped,
                         num_draws = 150L,
                         normalize_to_tmrel = FALSE)

# visual check draws
plot_BPRFmodel(data_info,draws,
               uncertainty_info,
               linear_model,
               signal_model,
               uncertainty_info_fill,
               linear_model_fill)

title <- paste0(ro_pair, ": egger_mean=", round(egger_model$mean, 3),
                ", egger_sd=", round(egger_model$sd,3), ", egger_pval=", 
                round(egger_model$pval, 3))
funnel_plot<-plot_residual(df, title)

plot_BPRFmodel_rr(data_info,draws,
                  uncertainty_info,
                  linear_model,
                  signal_model,
                  uncertainty_info_fill,
                  linear_model_fill)

# summary
summary <- summarize_model(data_info,
                           uncertainty_info,
                           linear_model,
                           signal_model,
                           egger_model,
                           egger_model_all,
                           uncertainty_info_fill,
                           linear_model_fill)

write.csv(summary, paste0(ro_pair,"_summary.csv"), row.names=FALSE)