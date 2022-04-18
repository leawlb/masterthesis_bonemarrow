# this script performs dimensionality reduction on mmus CPI of different conds.

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

# load cpi objects of each condition
mmus_yng_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi"))

mmus_old_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi"))

mmus_yng_cpi_raw <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi_raw"))

mmus_old_cpi_raw <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi_raw"))

mmus_yng_cpi_nlv <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi_nlv"))

mmus_old_cpi_nlv <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi_nlv"))

mmus_yng_cpi_nct <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi_nct"))

mmus_old_cpi_nct <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi_nct"))

mmus_yng_cpi_nrm <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi_nrm"))

mmus_old_cpi_nrm <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi_nrm"))

mmus_yng_cpi_nds <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi_nds"))

mmus_old_cpi_nds <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi_nds"))

mmus_yng_cpi_rev <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_yng_cpi_rev"))

mmus_old_cpi_rev <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_old_cpi_rev"))

#-------------------------------------------------------------------------------

# pre-processing
cpi_list <- list(mmus_yng_cpi,
                 mmus_old_cpi,
                 mmus_yng_cpi_raw,
                 mmus_old_cpi_raw,
                 mmus_yng_cpi_nlv,
                 mmus_old_cpi_nlv,
                 mmus_yng_cpi_nct,
                 mmus_old_cpi_nct,
                 mmus_yng_cpi_nrm,
                 mmus_old_cpi_nrm,
                 mmus_yng_cpi_nds,
                 mmus_old_cpi_nds,
                 mmus_yng_cpi_rev,
                 mmus_old_cpi_rev)

names <- c("mmus_yng_cpi",
           "mmus_old_cpi",
           "mmus_yng_cpi_raw",
           "mmus_old_cpi_raw",
           "mmus_yng_cpi_nlv",
           "mmus_old_cpi_nlv",
           "mmus_yng_cpi_nct",
           "mmus_old_cpi_nct",
           "mmus_yng_cpi_nrm",
           "mmus_old_cpi_nrm",
           "mmus_yng_cpi_nds",
           "mmus_old_cpi_nds",
           "mmus_yng_cpi_rev",
           "mmus_old_cpi_rev")

names(cpi_list) <- names

# remove NAs or PCA functions won't work
cpi_list <- lapply(cpi_list, function(x){
  x$Score[is.na(x$Score)] <- 0
  return(x)
})

#-------------------------------------------------------------------------------

# add suffix to allow merging of otherwise identical colnames (ctp)

cpi_list_tomerge <- cpi_list

for(i in 1:length(cpi_list_tomerge)){
  
  # add suffix to ctp names 
  suffix <- paste0("_", gsub("_cpi", "", names(cpi_list)[i]))
  colnames(cpi_list_tomerge[[i]]$Score) <- paste(
    colnames(cpi_list_tomerge[[i]]$Score), suffix, sep = "")
}

#-------------------------------------------------------------------------------

# manually merge both age groups per condition

# all pre-processing steps 
mmus_cpi <- cpi_list_tomerge[[1]]

mmus_cpi$Score <- cbind(
  cpi_list_tomerge[[1]]$Score, cpi_list_tomerge[[2]]$Score)
mmus_cpi$Receptorrank <- cbind(
  cpi_list_tomerge[[1]]$Receptorrank, cpi_list_tomerge[[2]]$Receptorrank)
mmus_cpi$Ligandrank <- cbind(
  cpi_list_tomerge[[1]]$Ligandrank, cpi_list_tomerge[[2]]$Ligandrank)
mmus_cpi$Celltypes <- rbind(
  cpi_list_tomerge[[1]]$Celltypes, cpi_list_tomerge[[2]]$Celltypes)

# no pre-processing steps
mmus_cpi_raw <- cpi_list_tomerge[[3]]

mmus_cpi_raw$Score <- cbind(
  cpi_list_tomerge[[3]]$Score, cpi_list_tomerge[[4]]$Score)
mmus_cpi_raw$Receptorrank <- cbind(
  cpi_list_tomerge[[3]]$Receptorrank, cpi_list_tomerge[[4]]$Receptorrank)
mmus_cpi_raw$Ligandrank <- cbind(
  cpi_list_tomerge[[3]]$Ligandrank, cpi_list_tomerge[[4]]$Ligandrank)
mmus_cpi_raw$Celltypes <- rbind(
  cpi_list_tomerge[[3]]$Celltypes, cpi_list_tomerge[[4]]$Celltypes)

# no levelling
mmus_cpi_nlv <- cpi_list_tomerge[[5]]

mmus_cpi_nlv$Score <- cbind(
  cpi_list_tomerge[[5]]$Score, cpi_list_tomerge[[6]]$Score)
mmus_cpi_nlv$Receptorrank <- cbind(
  cpi_list_tomerge[[5]]$Receptorrank, cpi_list_tomerge[[6]]$Receptorrank)
mmus_cpi_nlv$Ligandrank <- cbind(
  cpi_list_tomerge[[5]]$Ligandrank, cpi_list_tomerge[[6]]$Ligandrank)
mmus_cpi_nlv$Celltypes <- rbind(
  cpi_list_tomerge[[5]]$Celltypes, cpi_list_tomerge[[6]]$Celltypes)

# no cut-off
mmus_cpi_nct <- cpi_list_tomerge[[7]]

mmus_cpi_nct$Score <- cbind(
  cpi_list_tomerge[[7]]$Score, cpi_list_tomerge[[8]]$Score)
mmus_cpi_nct$Receptorrank <- cbind(
  cpi_list_tomerge[[7]]$Receptorrank, cpi_list_tomerge[[8]]$Receptorrank)
mmus_cpi_nct$Ligandrank <- cbind(
  cpi_list_tomerge[[7]]$Ligandrank, cpi_list_tomerge[[8]]$Ligandrank)
mmus_cpi_nct$Celltypes <- rbind(
  cpi_list_tomerge[[7]]$Celltypes, cpi_list_tomerge[[8]]$Celltypes)

# no cell type removal
mmus_cpi_nrm <- cpi_list_tomerge[[9]]

mmus_cpi_nrm$Score <- cbind(
  cpi_list_tomerge[[9]]$Score, cpi_list_tomerge[[10]]$Score)
mmus_cpi_nrm$Receptorrank <- cbind(
  cpi_list_tomerge[[9]]$Receptorrank, cpi_list_tomerge[[10]]$Receptorrank)
mmus_cpi_nrm$Ligandrank <- cbind(
  cpi_list_tomerge[[9]]$Ligandrank, cpi_list_tomerge[[10]]$Ligandrank)
mmus_cpi_nrm$Celltypes <- rbind(
  cpi_list_tomerge[[9]]$Celltypes, cpi_list_tomerge[[10]]$Celltypes)

# no downsampling
mmus_cpi_nds <- cpi_list_tomerge[[11]]

mmus_cpi_nds$Score <- cbind(
  cpi_list_tomerge[[11]]$Score, cpi_list_tomerge[[12]]$Score)
mmus_cpi_nds$Receptorrank <- cbind(
  cpi_list_tomerge[[11]]$Receptorrank, cpi_list_tomerge[[12]]$Receptorrank)
mmus_cpi_nds$Ligandrank <- cbind(
  cpi_list_tomerge[[11]]$Ligandrank, cpi_list_tomerge[[12]]$Ligandrank)
mmus_cpi_nds$Celltypes <- rbind(
  cpi_list_tomerge[[11]]$Celltypes, cpi_list_tomerge[[12]]$Celltypes)

# reverse classification
mmus_cpi_rev <- cpi_list_tomerge[[13]]

mmus_cpi_rev$Score <- cbind(
  cpi_list_tomerge[[13]]$Score, cpi_list_tomerge[[14]]$Score)
mmus_cpi_rev$Receptorrank <- cbind(
  cpi_list_tomerge[[13]]$Receptorrank, cpi_list_tomerge[[14]]$Receptorrank)
mmus_cpi_rev$Ligandrank <- cbind(
  cpi_list_tomerge[[13]]$Ligandrank, cpi_list_tomerge[[14]]$Ligandrank)
mmus_cpi_rev$Celltypes <- rbind(
  cpi_list_tomerge[[13]]$Celltypes, cpi_list_tomerge[[14]]$Celltypes)

# put in new list
cpi_list_comp_age <- list(
  mmus_cpi,
  mmus_cpi_raw,
  mmus_cpi_nlv,
  mmus_cpi_nct,
  mmus_cpi_nrm,
  mmus_cpi_nds,
  mmus_cpi_rev
)

names(cpi_list_comp_age) <- c(
  "mmus_cpi",
  "mmus_cpi_raw",
  "mmus_cpi_nlv",
  "mmus_cpi_nct",
  "mmus_cpi_nrm",
  "mmus_cpi_nds",
  "mmus_cpi_rev"
)

#-------------------------------------------------------------------------------
# Dimensionality reduction
#-------------------------------------------------------------------------------

# use get_pca function from cpi_functions_analysis.R

#-------------------------------------------------------------------------------
# all objects separately

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

#-------------------------------------------------------------------------------
# to compare age groups within conditions

cpi_list_comp_age <- lapply(cpi_list_comp_age, function(cpi){
  
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

#-------------------------------------------------------------------------------

# save all 
for(i in 1:length(cpi_list)){
  saveRDS(cpi_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list)[i], "_pca"))
}

for(i in 1:length(cpi_list_comp_age)){
  saveRDS(cpi_list_comp_age[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_comp_age)[i], "_pca"))
}
