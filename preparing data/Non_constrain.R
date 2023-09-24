###Data cleaning###
rm(list = ls())
library(reticulate)
library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
WORK_DIR <- "/home/hopehappy/BPRF_wholegrain"
setwd(WORK_DIR)

source(paste0(WORK_DIR, "/src/utils/prep_diet_data_function.R"))
source(paste0(WORK_DIR, "/src/utils/qsub_function.R"))
OUT_DIR <- c("/home/hopehappy/BPRF_wholegrain")
INPUT_DATA_DIR <- c("/mnt/c/Users/ctbhe/Desktop/T2D/BPRF_OUT")

####mrtool module import####
mrtool<-import("mrtool")
evidence_score <- import("mrtool.evidence_score.scorelator")
limetr<-import("limetr")
utils<-import("mrtool.core.utils")
knots<-import("mrtool.core.model")
core<-import("mrtool.core")
np <- import("numpy")

ro_pair<-"CRC"
data <- readRDS(paste0(INPUT_DATA_DIR, "/00_prepped_data/", ro_pair, ".RDS"))
df <- data$df
###1 calculate signal model####
#setting global variable
OBS_VAR <- "ln_effect"
OBS_SE_VAR <- "ln_se"
STUDY_ID_VAR <- "nid"
REF_EXPOSURE_COLS <- c("a_0", "a_1")
ALT_EXPOSURE_COLS <- c("b_0", "b_1")
cov_names <- c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)

# Specify all the columns you need for your application
mrdata <- mrtool$MRData()
mrdata$load_df(
  data =df, 
  col_obs = OBS_VAR,
  col_obs_se = OBS_SE_VAR, 
  col_study_id = STUDY_ID_VAR, 
  col_covs = as.list(cov_names)
)

#monotonicity <- 'increasing'
monotonicity <- NULL
#monotonicity <- c("decreasing")
degree <- 3L
N_I_KNOTS <- 3L
PRIOR_VAR_RSLOPE = 1e-6
PRIOR_VAR_MAXDER <- 1e-4

ensemble_cov_model <- mrtool$LogCovModel(
  alt_cov = ALT_EXPOSURE_COLS,
  ref_cov = REF_EXPOSURE_COLS,
  use_spline = TRUE,
  use_re = FALSE,
  spline_degree = degree,
  spline_knots_type = 'domain',
  spline_r_linear = TRUE,
  prior_spline_funval_uniform = array(c(-1 + 1e-6, 19)),
  prior_spline_num_constraint_points = 150L,
  spline_knots = array(seq(0, 1, length.out = N_I_KNOTS + 2)),
  prior_spline_maxder_gaussian = cbind(rbind(rep(0, N_I_KNOTS),
                                             rep(sqrt(PRIOR_VAR_MAXDER), N_I_KNOTS)), c(0, sqrt(PRIOR_VAR_RSLOPE))),
  prior_spline_monotonicity = monotonicity,
  prior_spline_der2val_gaussian_domain = array(c(0.0, 1.0))
)

# Create knot samples
knots <- import("mrtool.core.model")
knots_samples <- knots$create_knots_samples(
  data = mrdata, l_zero = TRUE, num_splines = 50L, 
  num_knots = 5L, width_pct = 0.2,
  alt_cov_names = ALT_EXPOSURE_COLS,
  ref_cov_names = REF_EXPOSURE_COLS)
# Ensemble model with exposure only 
signal_model <- mrtool$MRBeRT(mrdata,
                              ensemble_cov_model=ensemble_cov_model,
                              ensemble_knots=knots_samples)

signal_model$fit_model(inner_print_level=5L, inner_max_iter=200L,outer_step_size=200L, outer_max_iter=100L)

# create "new covariates" for later use

signal <- signal_model$predict(mrdata, predict_for_study=FALSE)
plot(exp(signal), type="l" )
plot(signal,type="l")

###test shape######
#df_signal_pred <- data.table(a_0 = 0,
#                             a_1 = 0,
#                             b_0 = 0:400,
#                             b_1 = 0:400)
#data_signal_pred <- mrtool$MRData()
#data_signal_pred$load_df(
#  data = df_signal_pred,
#  col_covs = as.list(names(df_signal_pred))
#)
#signal_plot <- signal_model$predict(data_signal_pred, predict_for_study=FALSE)

#plot(exp(signal_plot), type="l" )
#evidence_BPRF<-paste(ro_pair,"_check.tiff")
#  tiff(filename = evidence_BPRF,
#     width =1000, height = 850, units = "px", pointsize = 12,res=100)
#plot(df_signal_pred$b_1,signal_plot,
#    ylim = c(-1.0,0.0),
#      xlab = "exposure", ylab = "ln_relative_risk", pch = 19)
#dev.off()


# Extract weights of data point
w <- t(do.call(rbind, 
               lapply(1:length(signal_model$sub_models), 
                      function(i){signal_model$sub_models[[i]]$w_soln}))
) %*% signal_model$weights
#check draws
dev.new
with(df,plot(b_1,ln_effect))
with(df, lines(b_1, signal))

df_data <-  mrdata$to_df()

# Assign signal to data for use in later stage

df_data$signal <-  signal

# Drop data trimmed
df_data <-  df_data[w >= 0.1,]

# Save data and model
py_save_object(object = signal_model, 
               filename = paste0(OUT_DIR, "/01_template_pkl_files/",ro_pair, "_signal_model", "_trimmed.pkl"), 
               pickle = "dill")

out <- append(data, list(df_data=df_data))

saveRDS(out, paste0(OUT_DIR, "/01_template_models/",ro_pair, "_signal_model", "_trimmed.RDS"))

# 02_covariate_selection.R
#-------------------------------------------------------------------------------------------------------------------------------------------------
#
###covariates selected and create draws####
###ro_pair is aim of analysis 

library(dplyr)
np<-import("numpy")
data <- readRDS(paste0(OUT_DIR, "/01_template_models/",ro_pair, "_signal_model", "_trimmed.RDS"))
df_data <- data$df_data
mrdata <- mrtool$MRData()
mrdata$load_df(
  data = df_data, 
  col_obs = c('obs'),
  col_obs_se = c('obs_se'), 
  col_study_id = c('study_id'), 
  col_covs = as.list(c("signal"))
)

# Fit Linear Cov model with signal
cov_models <- list(mrtool$LinearCovModel(
  alt_cov = "signal",
  use_re = TRUE, 
  prior_beta_uniform=array(c(1.0, 1.0))
))

# No trimming
model <- mrtool$MRBRT(
  data = mrdata,
  cov_models = cov_models,
  inlier_pct = 1.0
)

model$fit_model(inner_print_level=5L, inner_max_iter=200L, 
                outer_step_size=200L, outer_max_iter=100L) 

# Sample betas to use as priors for covariate selection.
sampling <- import("mrtool.core.other_sampling")
beta_samples <- sampling$sample_simple_lme_beta(1000L, model)
beta_std <- sd(beta_samples)


# Save data and model
py_save_object(object = model, 
               filename = paste0(OUT_DIR, "/02_loglinear_pkl_files/",ro_pair, "_trimmed.pkl"), 
               pickle = "dill")

out <- append(data, list(beta_std = beta_std))
saveRDS(out, paste0(OUT_DIR, "/02_loglinear_models/",ro_pair, "_trimmed.RDS"))

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 03_covariate_selection.R
#------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)
BIAS_COVARIATES_AS_INTX <- TRUE
BETA_PRIOR_MULTIPLIER = 0.1

COV_FINDER_CONFIG = list(
  pre_selected_covs = list("exposure_linear"), 
  num_samples = 1000L,
  power_range = list(-4, 4), 
  power_step_size = 0.05,
  laplace_threshold = 1e-5,
  inlier_pct = 1.0, #since we trim in stage 1
  bias_zero = TRUE
)

# Read data
data <- readRDS(paste0(OUT_DIR, "/02_loglinear_models/",ro_pair, "_trimmed.RDS"))
df_data <- data$df_data
df_tmp <- data$df
df_tmp <- df_tmp[as.numeric(rownames(df_data)),]

cv_excl <- names(df_tmp)[sapply(df_tmp, setequal, 0) & grepl('cv_',names(df_tmp))]
cov_names <- c("exposure_linear", data$x_covs[!data$x_covs %in% eval(cv_excl)])

# Delete in future development; skip covariate selection for now. WHY? 
candidate_covs <- cov_names[!cov_names == "exposure_linear"]

# Interaction with signal
if (BIAS_COVARIATES_AS_INTX){
  for (cov in candidate_covs) df_data[, cov] <- df_data$signal * df_tmp[, cov]
}

# Change the name of signal to exposure_linear, since some
# underlying code deal with column name `exposure_linear`
df_data$exposure_linear <- df_data$signal
mrdata <- mrtool$MRData()

mrdata$load_df(
  data = df_data, 
  col_obs = c('obs'),
  col_obs_se = c('obs_se'), 
  col_study_id = c('study_id'), 
  col_covs = as.list(cov_names)
)

loglinear_model <- readRDS(paste0(OUT_DIR, "/02_loglinear_models/",ro_pair, "_trimmed.RDS"))

# Beta prior from first loglinear model results.
beta_gprior_std <- loglinear_model$beta_std
covfinder <- do.call(
  mrtool$CovFinder,
  c(COV_FINDER_CONFIG, 
    list(
      data = mrdata, 
      covs = as.list(candidate_covs)),
    beta_gprior_std = BETA_PRIOR_MULTIPLIER * beta_gprior_std
  )
)

#covfinder: There is no covariates to select, will return the pre-selected covariates.
covfinder$select_covs(verbose = TRUE)

selected_covs <- covfinder$selected_covs
selected_covs

# Save data and selected covariates
out <- append(data, list(df_cov_selection=df_data, selected_covs=selected_covs))
saveRDS(out, paste0(OUT_DIR, "/03_covariate_selection_models/", ro_pair, "_trimmed.RDS"))

#-------------------------------------------------------------------------------------------------------------------------------------------------
# 04_mixed_effects_models.R
#-------------------------------------------------------------------------------------------------------------------------------------------------
library(dplyr)

# Extract selected covariates
data <- readRDS(paste0(OUT_DIR, "/03_covariate_selection_models/", ro_pair,"_trimmed.RDS"))
##df_data is calculated  by signal model and filter w<0.1 (trimmed)
##df_tmp is raw data and not trimmed

df_data <- data$df_data
df_tmp <- data$df

# Only keep rows that are not trimmed
df_tmp <- df_tmp[as.numeric(rownames(df_data)),]

cov_names <- data$selected_covs
bias_covs <- cov_names[!cov_names == "exposure_linear"]

# Add interaction
for (cov in bias_covs) df_data[, cov] <- df_data$signal * df_tmp[, cov]

# Selected bias covariates plus signal
covs <- c("signal", bias_covs)

mrdata <- mrtool$MRData()
mrdata$load_df(
  df_data,
  col_obs = c('obs'),
  col_obs_se = c('obs_se'), 
  col_study_id = c('study_id'),
  col_covs=as.list(covs)
)

#loglinear_model <- readRDS(paste0(OUT_DIR, "/02_loglinear_models/", ro_pair, ".RDS"))

# Beta prior from first loglinear model results.
# beta_std drived from BPRF data

#beta_gprior_std <- loglinear_model$beta_std

# Combine cov models
cov_models <- list()
for (cov in bias_covs) cov_models <- append(cov_models, 
                                            list(
                                              do.call(
                                                mrtool$LinearCovModel, 
                                                list( prior_beta_gaussian=array(c(0, BETA_PRIOR_MULTIPLIER * data$beta_std)),
                                                      use_re = F,
                                                      alt_cov=cov
                                                )
                                              )
                                            )
)

# Mixed effects model
cov_models <- append(cov_models, mrtool$LinearCovModel('signal', 
                                                       use_re=TRUE,
                                                       prior_beta_uniform=array(c(1.0, 1.0))
)
)


model <- mrtool$MRBRT(
  data=mrdata,
  cov_models = cov_models,
  inlier_pct = 1.0
)

model$fit_model(inner_print_level=5L, inner_max_iter=200L, 
                outer_step_size=200L, outer_max_iter=100L)

# Load signal model and data in Stage 1
signal_model <- py_load_object(filename=paste0(OUT_DIR , "/01_template_pkl_files/",ro_pair, "_signal_model" , "_trimmed.pkl"), 
                               pickle = "dill")

orig_data <- readRDS(paste0(OUT_DIR, "/01_template_models/",ro_pair, "_signal_model", "_trimmed.RDS"))
###df orignal data without trimmed and filter
df <- orig_data$df

# This should be provided by the user
NUM_POINTS <- 200L
exposure_lower <- min(df[,c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)])
exposure_upper <- max(df[,c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS)])
#exposure <- seq(0, 400, length.out=NUM_POINTS)
exposure <- seq(exposure_lower, exposure_upper, length.out=NUM_POINTS)
min_cov <- rep(exposure_lower, NUM_POINTS)

# Deal with Sarah's data
if ('a_0' %in% REF_EXPOSURE_COLS){
  df_signal_pred <- data.frame(a_0=min_cov, a_1=min_cov, b_0=exposure, b_1=exposure)
} else {
  df_signal_pred <- data.frame(a_0=min_cov, b_0=exposure)
}

# Predict using signal model and gridded exposure
data_signal_pred <- mrtool$MRData()
data_signal_pred$load_df(
  df_signal_pred,
  col_covs = as.list(c(REF_EXPOSURE_COLS, ALT_EXPOSURE_COLS))
)
signal_pred <- signal_model$predict(data_signal_pred)

df_final_pred <- data.table(signal=signal_pred)

# add bias covs in the dataset (do not need this?)
#df_final_pred[, (bias_covs) := 0]

data_final_pred <- mrtool$MRData()
data_final_pred$load_df(
  df_final_pred,
  col_covs = as.list(c("signal"))
)

# create draws and prediction
set.seed(123)
sampling <- import("mrtool.core.other_sampling")
num_samples <- 1000L
beta_samples <- sampling$sample_simple_lme_beta(num_samples, model)
gamma_samples <- rep(model$gamma_soln, num_samples) * matrix(1, num_samples)



#data_final_pred$load_df(
#  df_final_pred,
#  col_covs = as.list(c("signal", bias_covs))
#)


curve <- model$predict(data_final_pred)
draws <- model$create_draws(
  data_final_pred,
  beta_samples=beta_samples,
  gamma_samples=gamma_samples
)


# Save model
py_save_object(object = model, 
               filename = paste0(OUT_DIR, "/04_mixed_effects_pkl_files/", ro_pair, "_linear_model", "_trimmed.pkl"), 
               pickle = "dill")