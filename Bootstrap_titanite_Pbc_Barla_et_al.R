###############################################
# Script for Barla et al. (submitted) paper
# in Chemical Geology
##############################################

# libraries
library(IsoplotR)
library("readxl")
library(dplyr)
library(plotly)

# set working directory
setwd("your/working/directory")

# import data and some try
df_in <- read.csv("your/csv/file/directory", header = T, sep=";")
# plot(df_in$Final207_206, df_in$Final238_206)
f207206 <- df_in$Final207_206
df_new <- cbind(df_in, f207206)

# add another for loop to iterate through subsample size
for (samp in seq(5, 50, by = 5))
{
  print("Run for subsample size : ")
  print(samp)
  
  # declare empty variables like lists
  val = 0
  mean_76_val_list <- c()
  mean_76_sd_list <- c()
  mean_76_ic_list <- c()
  mean_86_val_list <- c()
  mean_86_sd_list <- c()
  mean_86_ic_list <- c()
  li_age_list <- c()
  ui_val_list <- c()
  li_age_unc_list <- c()
  age_mswd_list <- c()
  chi_r_list <- c()
  mins_76_val_list <- c()
  maxs_76_val_list <- c()
  mins_86_val_list <- c()
  maxs_86_val_list <- c()
  vars_86_list <- c()
  vars_76_list <- c()
  cor_age_list <- c()
  cor_age_unc_list <- c()
  wm_cor_86age_list <- c()
  wm_cor_86age_unc_list <- c()
  f207_list <- c()
  f207std_list <- c()
  f207min_list <- c()
  f207max_list <- c()
  f207pct25_list <- c()
  f207pct75_list <- c()
  f207pct5_list <- c()
  f207pct95_list <- c()
  
  # start iterations for bootstrap (here 200 iterations)
  for (val in 1: 200)
  {
    df_samp <- sample_n(df_in, samp)
    
    #basic stats for 207-206 ratio using weighted mean function
    df_samp_76 <- df_samp[ , c("Final207_206", "Final207_206_Int2SE")]
    
    mean_76 = IsoplotR::weightedmean(df_samp_76, )
    mean_76_res <- mean_76["mean"]
    
    mean76_val <- mean_76_res[[1]][[1]]
    mean76_sd <- mean_76_res[[1]][[2]]
    var76 = var(df_samp$"Final207_206")
    
    # add el in list
    mean_76_val_list <- append(mean_76_val_list, mean76_val)
    mean_76_sd_list <- append(mean_76_sd_list, mean76_sd)
    min_76 <- min(df_samp$"Final207_206")
    max_76 <- max(df_samp$"Final207_206")
    mins_76_val_list <- append(mins_76_val_list, min_76)
    maxs_76_val_list <- append(maxs_76_val_list, max_76)
    vars_76_list <- append(vars_76_list, var76)
    
    #basic stats for 238-206 ratio using weighted mean function
    df_samp_86 <- df_samp[ , c("Final238_206", "Final238_206_Int2SE")]
    mean_86 = IsoplotR::weightedmean(df_samp_86, )
    mean_86_res <- mean_86["mean"]
    mean86_val <- mean_86_res[[1]][[1]]
    mean86_sd <- mean_86_res[[1]][[2]]
    var86 = var(df_samp$"Final238_206")
    
    # add el in list
    mean_86_val_list <- append(mean_86_val_list, mean86_val[1])
    mean_86_sd_list <- append(mean_86_sd_list, mean86_sd)
    
    min_86 <- min(df_samp$"Final238_206")
    max_86 <- max(df_samp$"Final238_206")
    mins_86_val_list <- append(mins_86_val_list, min_86)
    maxs_86_val_list <- append(maxs_86_val_list, max_86)
    vars_86_list <- append(vars_86_list, var86)
    
    # calculate low intercept age using isoplotR
    df_upb <- IsoplotR::read.data(df_samp, method = "U-Pb", format=2)
    age <- IsoplotR::age(df_upb, type=3)

    li_age <- age["par"][[1]][[1]]
    ui_val <- age["par"][[1]][[2]]
    li_age_unc <- age["err"][[1]][2]
    age_mswd <- age["mswd"][[1]]
    chi_r <- age["value"][[1]]
    
    # add el in list
    li_age_list <- append(li_age_list, li_age)
    ui_val_list <- append(ui_val_list, ui_val)
    li_age_unc_list <- append(li_age_unc_list, li_age_unc)
    age_mswd_list <- append(age_mswd_list, age_mswd)
    chi_r_list <- append(chi_r_list, chi_r)
    
    # calc corr age as concordia
    pb_corr <- IsoplotR::Pb0corr(df_upb, option = 2, omit4c = NULL)
    age_corr2 <- IsoplotR::age(pb_corr, type=2)
    cor_age <- age_corr2["age"][[1]][[1]]
    cor_age_unc <- age_corr2["age"][[1]][[2]]
    cor_age_list <- append(cor_age_list, cor_age)
    cor_age_unc_list <- append(cor_age_unc_list, cor_age_unc)
    
    df_orig <- pb_corr$x.raw
    df_orig2 = as.data.frame(df_orig)
    df_corr <- pb_corr$x
    df_corr2 = as.data.frame(df_corr)
    df_corr2$"pbc" <- ui_val
    
    # calc corr age as weighted mean 238-206
    wm_86_age <- IsoplotR::weightedmean(pb_corr, detect.outliers = TRUE, plot = TRUE, type=2)
    wm_cor_86age <- wm_86_age["mean"][[1]][[1]]
    wm_cor_86age_unc <- wm_86_age["mean"][[1]][[2]]
    
    wm_cor_86age_list <- append(wm_cor_86age_list, wm_cor_86age)
    wm_cor_86age_unc_list <- append(wm_cor_86age_unc_list, wm_cor_86age_unc)
    
    df_corr2$"f207" <- (df_orig2$"Pb207Pb206" - df_corr2$"Pb207Pb206") / (df_corr2$"pbc" - df_corr2$"Pb207Pb206")
    df_corr2$"f206" <- df_corr2$"f206" * df_corr2$"pbc" / df_orig2$"Pb207Pb206"
    f206_min <- min(df_corr2$"f207")
    f206_max <- max(df_corr2$"f207")
    f207_perct75 <- quantile(df_corr2$"f207", probs = 0.75)
    f207_perct25 <- quantile(df_corr2$"f207", probs = 0.25)
    f207_perct95 <- quantile(df_corr2$"f207", probs = 0.95)
    f207_perct5 <- quantile(df_corr2$"f207", probs = 0.05)
    f207_mean_pct <- 100*mean(df_corr2$"f207")
    f207_stdv <- sd(df_corr2$"f207")
    
    # add el in list
    f207_list <- append(f207_list, f207_mean_pct)
    f207std_list <- append(f207std_list, f207_stdv)
    f207min_list <- append(f207min_list, f207_min)
    f207max_list <- append(f207max_list, f207_max)
    f207pct75_list <- append(f207pct75_list, f207_perct75)
    f207pct25_list <- append(f207pct25_list, f207_perct25)
    f207pct95_list <- append(f207pct95_list, f207_perct95)
    f207pct5_list <- append(f207pct5_list, f207_perct5)
  }
  
  print(mean_76_sd_list)
  
  file_out = paste('bra_out_', toString(samp), 'samp.csv', sep = '_')
  sav_path = paste("D:/_PAPERS/__Submitted/Titanite_UPb-TREE_dev_2019/bootstrap_project2/", file_out)
  df_out <- data.frame(li_age_list, li_age_unc_list, age_mswd_list, chi_r_list, ui_val_list, mean_86_val_list, mean_86_sd_list, mean_76_val_list, mean_76_sd_list, mins_86_val_list, maxs_86_val_list, mins_76_val_list, maxs_76_val_list, vars_86_list, vars_76_list, cor_age_list, cor_age_unc_list, wm_cor_86age_list, wm_cor_86age_unc_list, f207_list, f207std_list, f207min_list, f207max_list, f207pct75_list, f207pct25_list, f207pct5_list, f207pct95_list)
  
  print(df_out)
  write.csv(df_out,sav_path, row.names = FALSE)
  
}
