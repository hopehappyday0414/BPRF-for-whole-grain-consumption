#-----------------------------------------------------------------
# Purpose: cleaning the newly extracted whole grains data for mr-brt es
#-------------------------------------------------------------

rm(list=ls())
library(tidyverse)
library(data.table)
library(readxl)
library(stringr)

new_rr_dir <- c("/mnt/c/Users/ctbhe/Desktop/T2D/BPRF")
new_out_dir<-c("/mnt/c/Users/ctbhe/Desktop/T2D/BPRF_OUT")
new_rr_files <- list.files(new_rr_dir, full.names = T)
new_rr_files <- new_rr_files[!grepl("~$",new_rr_files, fixed = T) & !grepl("template|map",new_rr_files)]
new_rr_data <- lapply(new_rr_files,read_excel) %>% rbindlist(use.names = T, fill = T) %>% data.table
outcomes <- new_rr_data[, ro_pair] %>% unique %>% sort


new_rr_data_sub <- new_rr_data[ro_pair=="Stroke",]
#-----------------
## clean the data
#-----------------

confounders <- names(new_rr_data_sub)[grepl('confounders_', names(new_rr_data_sub)) & !grepl('_other', names(new_rr_data_sub))]

# change the data to numeric
int_vars <- c("age_start", "age_end", "rep_population", confounders,"pooled cohort")

num_vars <- c("age_mean", "age_sd", "percent_male", "follow_up", "effect_size", "lower", "upper", "custom_exp_level_lower", "custom_exp_level_upper",
              "custom_unexp_level_lower", "custom_unexp_level_upper")

new_rr_data_sub[, (int_vars) := lapply(.SD, str_trim), .SDcols=int_vars]
new_rr_data_sub[, (num_vars) := lapply(.SD, str_trim), .SDcols=num_vars]

new_rr_data_sub[, (int_vars) := lapply(.SD, as.integer), .SDcols=int_vars]
new_rr_data_sub[, (num_vars) := lapply(.SD, as.numeric), .SDcols=num_vars]

# remove former smoking exposure
#new_rr_data_sub <- new_rr_data_sub[!grepl("former", exp_temporality, ignore.case = T)]

#---------------------------
# log-transform the effects
#---------------------------

new_rr_data_sub[, effect_size_log := log(effect_size)]
new_rr_data_sub[, se_effect_log := abs((log(upper)-log(lower))/3.92)] 
new_rr_data_sub[, se_effect_log] %>% summary

# remove missing SE
new_rr_data_sub <- new_rr_data_sub[!is.na(effect_size_log)]
new_rr_data_sub <- new_rr_data_sub[!is.na(se_effect_log)]

# creating bias covariates
new_rr_data_sub[,cv_subpopulation := ifelse(rep_population==1, 0, 1)]
new_rr_data_sub[,cv_exposure_selfreport := ifelse(exp_method_1=="Self-report",1,0)]
new_rr_data_sub[,cv_outcome_selfreport := ifelse(outcome_assess_1=="Self-report",1,0)]
new_rr_data_sub[,cv_exposure_study := ifelse(exp_assess_period=="only at baseline",1,0)]

if(all(is.na(new_rr_data_sub[, age_start]))){
  new_rr_data_sub[!is.na(age_mean),cv_older := ifelse(age_mean>=50,1,0)]
} else {
  new_rr_data_sub[!is.na(age_start),cv_older := ifelse(age_start>=50,1,0)]
}

# new_rr_data_sub[,cv_non_smoker := ifelse(unexp_def=="never smokers",0,1)]
new_rr_data_sub[,cv_risk_measure := ifelse(effect_size_measure=="Relative risk (RR)" | effect_size_measure=="Hazard ratio (HR)",0,1)]

cvs <- names(new_rr_data_sub)[grepl('cv_', names(new_rr_data_sub))]

# for confounders, replace NA with 0
for (j in c(confounders)){
  set(new_rr_data_sub,which(is.na(new_rr_data_sub[[j]])),j,0)
}

# for cv covariates, replacing NA to 1 if necessary (very rare)
for (j in c(cvs)){
  set(new_rr_data_sub,which(is.na(new_rr_data_sub[[j]])),j,1)
}

# replace missing age_start and age_end with the median value (this is fine for non-cvd outcomes since age is not relavent. May be problematic for CVDs though)
set(new_rr_data_sub, which(is.na(new_rr_data_sub[["age_start"]])),"age_start", median(new_rr_data_sub[["age_start"]], na.rm = T))
set(new_rr_data_sub, which(is.na(new_rr_data_sub[["age_end"]])),"age_end", median(new_rr_data_sub[["age_end"]], na.rm = T))
new_rr_data_sub[is.na(percent_male), percent_male := 0.5]

# check the exposure upper missing values
new_rr_data_sub[, .(custom_exp_level_lower, custom_exp_level_upper, custom_unexp_level_lower, custom_unexp_level_upper)]
new_rr_data_sub[is.na(custom_exp_level_lower),]

#  new_rr_data_sub[is.na(custom_unexp_level_lower) & is.na(custom_unexp_level_upper) & !is.na(custom_unexp_level_mid), c("custom_unexp_level_lower", "custom_unexp_level_upper") := unexp_level_mid]

# change NA to 0 for unexposed group
for (j in c("custom_unexp_level_lower", "custom_unexp_level_upper")){
  set(new_rr_data_sub,which(is.na(new_rr_data_sub[[j]])),j,0)
}

# change NA in exp_upper to 1.5*lower if missing
new_rr_data_sub[is.na(custom_exp_level_upper), custom_exp_level_upper := 1.5*custom_exp_level_lower]

# count numbers of confounders
new_rr_data_sub[, adjustment := apply(.SD, 1, sum), .SDcols=confounders]

# count confounders in confounde_other
new_rr_data_sub[, confounders_other := gsub(",", ";", confounders_other)]

if(nrow(new_rr_data_sub[!is.na(confounders_other)])!=0){
  new_rr_data_sub[!is.na(confounders_other), num := unlist(lapply(str_split(confounders_other, ";"), length))]
}

new_rr_data_sub[is.na(confounders_other), num := 0]
new_rr_data_sub[, adjustment:= adjustment + num]

# creating cascading dummies
new_rr_data_sub[ ,cv_adj := as.numeric()]
new_rr_data_sub[,adj_energy := 0]
new_rr_data_sub[,adj_age_sex := 0]
new_rr_data_sub[confounders_age==1 & confounders_sex==1 , adj_age_sex := 1]
new_rr_data_sub[confounders_caloric==1, adj_energy := 1]

# cv_adj=2 if age or sex is not adjusted
new_rr_data_sub[adj_age_sex==0, cv_adj := 2]

# cv_adj=1 if only age and sex are adjusted
new_rr_data_sub[adj_age_sex==1 & adj_energy == 0, cv_adj := 1]

# cv_adj=1 if age+sex + caloric and less than 3 covariates are adjusted
new_rr_data_sub[adj_age_sex==1 & adj_energy == 1 & adjustment>3 , cv_adj := 0]


# check whether there is missing in cv_adj
message("there is ", nrow(new_rr_data_sub[is.na(cv_adj)]), " missing values in cv_adj")



# remove cv_adj
# new_rr_data_sub[, cv_adj := NULL]

# change the exposed and unexposed name
new_rr_data_sub <- setnames(new_rr_data_sub, old = c("custom_exp_level_lower", "custom_exp_level_upper", "custom_unexp_level_lower", "custom_unexp_level_upper"), 
                            new = c("b_0", "b_1", "a_0", "a_1"))

new_rr_data_sub <- setnames(new_rr_data_sub, old = c("effect_size_log", "se_effect_log"), 
                            new = c("ln_effect", "ln_se"))

# add mean exposure
new_rr_data_sub[, mean_exp := (b_0+b_1)/2]


# check data before saving
bias_covs <- names(new_rr_data_sub)[grepl('cv_', names(new_rr_data_sub))]

#cvd_ro <- c("ihd", "stroke", "afib_and_flutter", "peripheral_artery_disease", "aortic_aneurism")


incl_vars <- c("nid", "ln_effect", "ln_se", "b_0","b_1","a_0","a_1", "mean_exp", 'percent_male',"follow_up", bias_covs)
nesc_vars <- c("nid", "ln_effect", "ln_se", "b_0","b_1","a_0","a_1", "mean_exp", 'percent_male',"follow_up", bias_covs)


#keep the validation of input data
if(nrow(new_rr_data_sub[ln_se<0,])>0) message(paste0("ln_se < 0 for ", nrow(new_rr_data_sub[ln_se<0,]), " row" ))
if(nrow(new_rr_data_sub[lower>upper,])>0) message(paste0("lower and upper swapped for ", nrow(new_rr_data_sub[lower>upper,]), " rows" ))


if(new_rr_data_sub[, nesc_vars, with=F] %>% is.na %>% any){
  message("there are missing values in the required variables... please check")
} else {
  message("there are no missing values! Good to go!")
}

# keep usable variables
new_data <- new_rr_data_sub[,incl_vars, with=F]


#-------------------------------------------------------------------------------------------------------------------------------------------------
# save the .RDS file
#-------------------------------------------------------------------------------------------------------------------------------------------------
REF_EXPOSURE_COLS <- c("a_0", "a_1")
ALT_EXPOSURE_COLS <- c("b_0", "b_1")
obs_var <- "ln_effect"; obs_se_var <- "ln_se"
ref_vars <- c("a_0", "a_1"); alt_vars <- c("b_0", "b_1")
allow_ref_gt_alt = FALSE
drop_x_covs <- drop_z_covs <- keep_x_covs <- keep_z_covs <- NA
study_id_var <- "nid"
verbose <- TRUE
dt_sub_new<-new_rr_data_sub

# NOTE: these variables cannot be missing! ##
no_missing_vars <- c("nid", "ln_effect", "ln_se", alt_vars, ref_vars, 'percent_male', "follow_up", bias_covs)

# check the missingness in the no_missing_vars
if(dt_sub_new[, no_missing_vars, with=F] %>% is.na %>% any){
  message("there are missing values in the required variables... please check")
} else {
  message("there are no missing values! Good to go!")
}

#remove the cv_exposure_selfreport
#dt_sub_new[, cv_exposure_selfreport := NULL]

write.csv(dt_sub_new, file.path(new_out_dir, paste0(ro_pair, '.csv')), row.names = F)


# NOTE: this function will drop any incomplete rows. No missing values is allow.  
out <- prep_diet_data(ro_pair, obs_var, obs_se_var, ref_vars, alt_vars, allow_ref_gt_alt = FALSE,
                      study_id_var = study_id_var,
                      drop_x_covs = NA, keep_x_covs = NA, drop_z_covs = NA, keep_z_covs = NA,
                      diet_dir = new_out_dir,
                      verbose = TRUE)
print(out)
message(paste0('saving cleaned data..'))
saveRDS(out, paste0(new_out_dir, "/00_prepped_data/", ro_pair, ".RDS"))
message(paste0('saved cleaned data to: ', paste0(new_out_dir, "/00_prepped_data/",ro_pair, ".RDS")))