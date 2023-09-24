# =================================================================================================
# all_csvs.R
#=================================================================================================
# process results value

OUT_DIR<- "/home/hopehappy/BPRF_wholegrain/Results_csv/"


pkl_dir <- "/home/hopehappy/BPRF_wholegrain"

rr_files<-c("CRC","IHD","Stroke","T2D")


for (ro in rr_files) {
  print(ro)
  signal_model_path <- paste0(pkl_dir,"/01_template_pkl_files/", ro, "_signal_model", "_trimmed.pkl")
  linear_model_path <- paste0(pkl_dir, "/04_mixed_effects_pkl_files/", ro, "_linear_model", "_trimmed.pkl")
  
  signal_model <- py_load_object(filename = signal_model_path)
  linear_model <- py_load_object(filename = linear_model_path, pickle = "dill")
  
  # specify risk, you need to input the exposures that you want to predict

  risk <- seq(from=0, to=400, by=2.67)

  get_draws
  set.seed(123)
  np$random$seed(as.integer(321))
  
  ###save conservative estimation####
  df_re <- get_ln_rr_draws(signal_model,
                           linear_model,
                           risk,
                           num_draws = 170L,
                           normalize_to_tmrel = FALSE,
                           include_re = T)
  
  write.csv(df_re, paste0(OUT_DIR,ro,"_re.csv"), row.names = F)
  
  ###save contraditional estimation####
  set.seed(123)
  np$random$seed(as.integer(321))
  df_fe <- get_ln_rr_draws(signal_model,
                           linear_model,
                           risk,
                           num_draws = 170L,
                           normalize_to_tmrel = FALSE,
                           include_re = F)
  
  write.csv(df_fe, paste0(OUT_DIR,ro,"_fe.csv"), row.names = F)
  
  # visual check draws
  draws <- df_re[, 2:ncol(df_re)]
  draw_mean <- apply(draws, 1, function(x) mean(x))
  draw_lower <- apply(draws, 1, function(x) quantile(x, probs=.025))
  draw_upper <- apply(draws, 1, function(x) quantile(x, probs=.975))
  ###save_RR conservative value
  summary_rr <- cbind(risk, round(exp(draw_mean), digits = 4), round(exp(draw_lower), digits = 4), round(exp(draw_upper), digits = 4))
  write.csv(summary_rr, paste0(OUT_DIR,ro,'_summary_re.csv'), row.names = F)
  
  
  ###save_RR contraditional value
  # traditional 95% CI
  draws_fe <- df_fe[, 2:ncol(df_fe)]
  draw_mean_fe <- apply(draws_fe, 1, function(x) mean(x))
  draw_lower_fe <- apply(draws_fe, 1, function(x) quantile(x, probs=.025))
  draw_upper_fe <- apply(draws_fe, 1, function(x) quantile(x, probs=.975))
  
  summary_rr_fe <- cbind(risk, round(exp(draw_mean_fe), digits = 4), round(exp(draw_lower_fe), digits = 4), round(exp(draw_upper_fe), digits = 4))
  write.csv(summary_rr_fe, paste0(OUT_DIR,ro,'_summary_fe.csv'), row.names = F)
}