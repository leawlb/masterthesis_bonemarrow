# this script calculates scoring variance for mmus replicates

#-------------------------------------------------------------------------------
# PREPARATION
#-------------------------------------------------------------------------------

# load libraries, source code
library(tidyverse)
library(rlist)

setwd("/home/l012t/Interspecies_BM")

source(file = "code/filepath.R")
source(file = "code/colors.R")
source(file = "interactome/functions/cpi_functions_analysis.R")

#-------------------------------------------------------------------------------

# load metrics objects (only ctpints required)
ctpints_yng1 <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng1_cpi_ctpints"))

ctpints_yng2 <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng2_cpi_ctpints"))

ctpints_yng3 <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng3_cpi_ctpints"))

ctpints_old2 <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old2_cpi_ctpints"))

ctpints_old3 <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old3_cpi_ctpints"))

#-------------------------------------------------------------------------------

ctpints_list <- list(
  ctpints_yng1,
  ctpints_yng2,
  ctpints_yng3,
  ctpints_old2,
  ctpints_old3)

names <- c(
  "mmus_yng1_ctpints",
  "mmus_yng2_ctpints",
  "mmus_yng3_ctpints",
  "mmus_old2_ctpints",
  "mmus_old3_ctpints")

names(ctpints_list) <- names

#-------------------------------------------------------------------------------

# make a DF to store the score percentiles 
score_df_list <- lapply(ctpints_list, function(ctpints){
  
  # start with the first element of the "separate" list inctpints
  score_df <- ctpints[[2]][[1]]
  score_df$score_percentile <- vector(length = nrow(score_df))
  score_percentile <- ecdf(score_df$score)
  
  for(j in 1:nrow(score_df)){
    score_df$score_percentile[j] <- score_percentile(score_df$score[j])
  }
  
  # now repeat for all other cell types starting from that
  for(i in 2:length(ctpints[[2]])){
    
    temp_df <- ctpints[[2]][[i]]
    temp_df$score_percentile <- vector(length = nrow(temp_df))
    score_percentile <- ecdf(temp_df$score)
    
    for(j in 1:nrow(temp_df)){
      temp_df$score_percentile[j] <- score_percentile(temp_df$score[j])
    }
    
    score_df <- rbind(score_df, temp_df)
  }
  
  # add metadata
  score_df$emitter <- factor(score_df$emitter, levels = names(col_emi_red))
  score_df$receiver <- factor(score_df$receiver, levels = names(col_rec_red))
  score_df$score_percentile <- score_df$score_percentile * 100
  
  return(score_df)
})

#-------------------------------------------------------------------------------

# calculate score variance 
score_df_list <- lapply(score_df_list, function(scdf){
  
  # calculate the variance of LRI scores for each interaction per cell type 
  scdf$emi_var <- vector(length = nrow(scdf))
  scdf$rec_var <- vector(length = nrow(scdf))
  
  # for each emitter
  for(i in scdf$emitter){
    
    # for each group (each LRI) calculate variance (but still within emitter)
    temp_df <- scdf[scdf$emitter == i,]
    temp_df <- group_by(temp_df, interaction)
    var_df <- temp_df %>% 
      summarise(var_int = var(score_percentile))
    
    # fill into score_df
    temp_df$emi_var <- var_df$var_int[
      match(temp_df$interaction, var_df$interaction)]
    scdf[scdf$emitter == i,] <- temp_df
    
  }
  
  # for each receiver
  for(i in scdf$receiver){
    
    temp_df <- scdf[scdf$receiver == i,]
    temp_df <- group_by(temp_df, interaction)
    var_df <- temp_df %>% 
      summarise(var_int = var(score_percentile))
    
    temp_df$rec_var <- var_df$var_int[
      match(temp_df$interaction, var_df$interaction)]
    
    scdf[scdf$receiver == i,] <- temp_df
    
  }
  # NAs are generated when there is an interaction detected only in one ctp/ct
  
  return(scdf)
})

names(score_df_list) <- names

#-------------------------------------------------------------------------------

# save it 
for(i in 1:length(score_df_list)){
  saveRDS(score_df_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(score_df_list)[i], "_score_df"))
}
