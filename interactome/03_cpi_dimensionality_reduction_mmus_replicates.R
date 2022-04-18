# this script performs dimensionality reduction on mmus replicate CPI objects

#-------------------------------------------------------------------------------
# PREPARATION
#-------------------------------------------------------------------------------

# load libraries, source code
library(tidyverse)
library(rlist)

setwd("/home/l012t/Interspecies_BM")

source(file = "code/filepath.R")
source(file = "interactome/functions/cpi_functions_analysis.R")

#-------------------------------------------------------------------------------

# load cpi objects
mmus_yng1_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_young1_cpi"))

mmus_yng2_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_young2_cpi"))

mmus_yng3_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_young3_cpi"))

mmus_old2_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old2_cpi"))

mmus_old3_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old3_cpi"))

#-------------------------------------------------------------------------------

cpi_list <- list(mmus_yng1_cpi,
                 mmus_yng2_cpi,
                 mmus_yng3_cpi,
                 mmus_old2_cpi,
                 mmus_old3_cpi)

names <- c("mmus_yng1_cpi",
           "mmus_yng2_cpi",
           "mmus_yng3_cpi",
           "mmus_old2_cpi",
           "mmus_old3_cpi")

names(cpi_list) <- names

# remove NAs or PCA functions won't work
cpi_list <- lapply(cpi_list, function(x){
  x$Score[is.na(x$Score)] <- 0
  return(x)
})

#-------------------------------------------------------------------------------

# PCA on all objects separately (no merging required)

cpi_list <- lapply(cpi_list, function(cpi){
  
  # perform pca
  cpi_pca <- get_pca(cpi)
  
  # collect info to df 
  pca_df <- as.data.frame(cpi_pca$ind$coord[,1:3])
  # add info on the percent of variation explained by each PC for top 5 PCs
  pca_df <- rbind(pca_df, c(cpi_pca$eig[1,2],
                            cpi_pca$eig[2,2],
                            cpi_pca$eig[3,2],
                            cpi_pca$eig[4,2],
                            cpi_pca$eig[5,2]))
  rownames(pca_df)[nrow(pca_df)] <- "Percentage of Variation"
  
  cpi$PCA <- pca_df
  return(cpi)
})

for(i in 1:length(cpi_list)){
  saveRDS(cpi_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list)[i], "_pca"))
}
