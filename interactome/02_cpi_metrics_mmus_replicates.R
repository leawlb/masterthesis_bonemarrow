# this script calculates important metrics for separate mmus replicates

#-------------------------------------------------------------------------------
# PREPARATION
#-------------------------------------------------------------------------------

species <- "mmus"

#-------------------------------------------------------------------------------

# load libraries, source code
library(SingleCellExperiment)
library(tidyverse)
library(rlist)

setwd("/home/l012t/Interspecies_BM")

source(file = "code/filepath.R")
source(file = "code/colors.R")
source(file = "interactome/functions/cpi_functions_analysis.R")

#-------------------------------------------------------------------------------
# load cpi objects of all replicates

mmus_yng1_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_young1_cpi"))

mmus_yng2_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_young2_cpi"))

mmus_yng3_cpi <- readRDS(file = paste0(
  filepath_lv1, "data/cci_files/mmus_young3_cpi"))

# mmus_old1 is missing a fraction

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

#-------------------------------------------------------------------------------
# calculate interactions per CTP

ctpints_list <- lapply(cpi_list, extract_ctp_info)
names(ctpints_list) <- names

# add info on age
for(i in 1:length(ctpints_list)){
  if(grepl("old", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$age <- "old"
  }else if(grepl("yng", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$age <- "young"
  }
}

# species
for(i in 1:length(ctpints_list)){
  ctpints_list[[i]][[1]]$species <- species
}

# factorize
ctpints_list <- lapply(ctpints_list, function(ctpints){
  
  ctpints$overview$receiver <- factor(ctpints$overview$receiver,
                                      levels = names(col_rec_red))
  ctpints$overview$emitter <- factor(ctpints$overview$emitter,
                                     levels = names(col_emi_red))
  ctpints$overview$age <- factor(ctpints$overview$age,
                                 levels = names(col_age))
  ctpints$overview$species <- factor(ctpints$overview$species,
                                     levels = names(col_spc))
  return(ctpints)
})

#-------------------------------------------------------------------------------
# calculate interactions per cell type

ctints_list <- lapply(ctpints_list, extract_ct_info)
names(ctints_list) <- names

# add info on age
for(i in 1:length(ctints_list)){
  if(grepl("old", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$age <- "old"
  }else if(grepl("yng", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$age <- "young"
  }
}

# species
for(i in 1:length(ctints_list)){
  ctints_list[[i]][[1]]$species <- species
}

# factorize
ctints_list <- lapply(ctints_list, function(ctints){
  
  ctints$overview$annotation <- factor(ctints$overview$annotation, 
                                       levels = names(col_ann))
  ctints$overview$age <- factor(ctints$overview$age,
                                levels = names(col_age))
  ctints$overview$species <- factor(ctints$overview$species,
                                    levels = names(col_spc))
  
  return(ctints)
})

#-------------------------------------------------------------------------------
# calculate nr of ligands/receptors per cell type

# this extracts all detected ligands or receptors per cell type
ctlrs_list <- lapply(ctints_list, function(ctints){
  ctlrs <- lapply(ctints[[2]], extract_lrs_info)
  return(ctlrs)
})

# this combines info from ctlrs_list and cpi_list to obtain an overview df
nrlrs_list <- extract_lrs_nrs(cpi = cpi_list, ctlrs = ctlrs_list)
names(nrlrs_list) <- names

# add info on age
for(i in 1:length(nrlrs_list)){
  if(grepl("old", names(nrlrs_list)[i])){
    nrlrs_list[[i]]$age <- "old"
  }else if(grepl("yng", names(nrlrs_list)[i])){
    nrlrs_list[[i]]$age <- "young"
  }
}

# species
for(i in 1:length(nrlrs_list)){
  nrlrs_list[[i]]$species <- species
}

# factorize
nrlrs_list <- lapply(nrlrs_list, function(nrlrs){
  
  nrlrs$annotation <- factor(nrlrs$annotation, levels = names(col_ann))
  nrlrs$age <- factor(nrlrs$age, levels = names(col_age))
  nrlrs$species <- factor(nrlrs$species, levels = names(col_spc))
  
  return(nrlrs)
})

#-------------------------------------------------------------------------------

# load all objects and subset them per replicate (same conditions)
mmus_yng_sce <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE_ds_ctrm"))
mmus_old_sce <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE_ds_ctrm"))

mmus_yng_sce1 <- mmus_yng_sce[,grep("mmus_young1", mmus_yng_sce$Sample)]
mmus_yng_sce2 <- mmus_yng_sce[,grep("mmus_young2", mmus_yng_sce$Sample)]
mmus_yng_sce3 <- mmus_yng_sce[,grep("mmus_young3", mmus_yng_sce$Sample)]

mmus_old_sce2 <- mmus_old_sce[,grep("mmus_old2", mmus_old_sce$Sample)]
mmus_old_sce3 <- mmus_old_sce[,grep("mmus_old3", mmus_old_sce$Sample)]

sce_list <- list(
  mmus_yng_sce1,
  mmus_yng_sce2,
  mmus_yng_sce3,
  
  mmus_old_sce2,
  mmus_old_sce3
)

nr_cells <- lapply(sce_list, function(sce){
  
  cts <- unique(sce$Identity)
  
  nr_cells <- vector()
  for(i in 1:length(cts)){
    nr_cells[i] <- length(which(sce$Identity == cts[i]))
  }
  
  names(nr_cells) <- cts
  return(nr_cells)
})

#-------------------------------------------------------------------------------

# add to the overview dfs
for(i in 1:length(sce_list)){
  
  cts <- names(nr_cells[[i]])
  
  # empty new vectors to fill for each list type
  ctpints_list[[i]]$overview$nr_cells_emi <- vector(
    length = nrow(ctpints_list[[i]]$overview))
  ctpints_list[[i]]$overview$nr_cells_rec <- vector(
    length = nrow(ctpints_list[[i]]$overview))
  
  ctints_list[[i]]$overview$nr_cells <- vector(
    length = nrow(ctints_list[[i]]$overview))
  
  nrlrs_list[[i]]$nr_cells <- vector(
    length = nrow(nrlrs_list[[i]]))
  
  for(j in 1:length(cts)){
    
    # proceed per list type
    if(cts[j] %in% ctpints_list[[i]]$overview$emitter){
      ctpints_list[[i]]$overview$nr_cells_emi[
        ctpints_list[[i]]$overview$emitter == cts[j]] <- nr_cells[[i]][
          grep(cts[j], names(nr_cells[[i]]))]
    }else if(cts[j] %in% ctpints_list[[i]]$overview$receiver){
      ctpints_list[[i]]$overview$nr_cells_rec[
        ctpints_list[[i]]$overview$receiver == cts[j]] <- nr_cells[[i]][
          grep(cts[j], names(nr_cells[[i]]))]
    }
    
    if(cts[j] %in% ctints_list[[i]]$overview$celltypes){
      ctints_list[[i]]$overview$nr_cells[
        ctints_list[[i]]$overview$celltypes == cts[j]] <- nr_cells[[i]][
          grep(cts[j], names(nr_cells[[i]]))]
    }
    
    if(cts[j] %in% nrlrs_list[[i]]$celltypes){
      nrlrs_list[[i]]$nr_cells[
        nrlrs_list[[i]]$celltypes == cts[j]] <- nr_cells[[i]][
          grep(cts[j], names(nr_cells[[i]]))]
    }
  }
}

#-------------------------------------------------------------------------------

for(i in 1:length(ctpints_list)){
  saveRDS(ctpints_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(ctpints_list)[i], "_ctpints"))
}

for(i in 1:length(ctints_list)){
  saveRDS(ctints_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(ctints_list)[i], "_ctints"))
}

for(i in 1:length(nrlrs_list)){
  saveRDS(nrlrs_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(nrlrs_list)[i], "_nrlrs"))
}
