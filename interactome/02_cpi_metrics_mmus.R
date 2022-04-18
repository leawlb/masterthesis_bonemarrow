# this script calculates important metrics for mmus for different CPI conditions

#-------------------------------------------------------------------------------
# PREPARATION
#-------------------------------------------------------------------------------

# set only_full to TRUE to only calculate for full cpi object
# set only_full to FALSE to only calculate for all conditions (takes long)

only_full <- FALSE
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

# pre-processing, put into list 

if(only_full == TRUE){
 
  cpi_list <- list(mmus_yng_cpi,
                   mmus_old_cpi)
  
  names <- c("mmus_yng_cpi",
             "mmus_old_cpi")
}else{

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
                   mmus_old_cpi_rev )

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
             "mmus_old_cpi_rev" )
}

names(cpi_list) <- names

#-------------------------------------------------------------------------------
# METRICS
#-------------------------------------------------------------------------------

# calculate interactions per CTP

ctpints_list <- list()
for(i in 1:length(cpi_list)){
  
  # in case of cutoff AND level == FALSE, method for "raw" objects must be used
  if(cpi_list[[i]]$Conditions["cutoff"] == FALSE){
    print(paste(names(cpi_list)[i], 
                "calculated by extract_ctp_info_raw", sep = " "))
    
    ctpints_list[[i]] <- extract_ctp_info_raw(
      cpi_list[[i]], cutoff = cpi_list[[i]]$Conditions["cutoff"])
    
  }else if(grepl("nrm", names(cpi_list)[[i]])){
    print(paste(names(cpi_list)[i], 
                "calculated by extract_ctp_info_raw", sep = " "))
    
    ctpints_list[[i]] <- extract_ctp_info_raw(
      cpi_list[[i]], cutoff = cpi_list[[i]]$Conditions["cutoff"])
  }
  else{
    
    # for these conditions  extract_ctp_info can be used which is much faster
    print(paste(names(cpi_list)[i], 
                "calculated by extract_ctp_info", sep = " "))
    
    ctpints_list[[i]] <- extract_ctp_info(cpi_list[[i]])
  }
}

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

# condition
for(i in 1:length(ctpints_list)){
  if(grepl("raw", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "no pre-processing steps"
  }else if(grepl("nlv", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "no levelling"
  }else if(grepl("nct", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "no cut-off"
  }else if(grepl("nrm", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "no cell type removal"
  }else if(grepl("nds", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "no downsampling"
  }else if(grepl("rev", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "reverse annotation"
  }else{
    ctpints_list[[i]][[1]]$condition <- "all pre-processing steps"
  }
}

# factorize cell types and age
ctpints_list <- lapply(ctpints_list, function(ctpints){
  
  ctpints$overview$receiver <- factor(
    ctpints$overview$receiver,
    levels = names(col_cts_all)[
      names(col_cts_all) %in% unique(ctpints$overview$receiver)])

  ctpints$overview$emitter <- factor(
    ctpints$overview$emitter,
    levels = names(col_cts_all)[
      names(col_cts_all) %in% unique(ctpints$overview$emitter)])

  ctpints$overview$age <- factor(ctpints$overview$age,
                                 levels = names(col_age))
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

# condition
for(i in 1:length(ctints_list)){
  if(grepl("raw", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$condition <- "no pre-processing steps"
  }else if(grepl("nlv", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$condition <- "no levelling"
  }else if(grepl("nct", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$condition <- "no cut-off"
  }else if(grepl("nrm", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$condition <- "no cell type removal"
  }else if(grepl("nds", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$condition <- "no downsampling"
  }else if(grepl("rev", names(ctints_list)[i])){
    ctints_list[[i]][[1]]$condition <- "reverse annotation"
  }else{
    ctints_list[[i]][[1]]$condition <- "all pre-processing steps"
  }
}

# factorize
ctints_list <- lapply(ctints_list, function(ctints){
  
  ctints$overview$celltypes <- factor(
    ctints$overview$celltypes,
    levels = names(col_cts_all)[
      names(col_cts_all) %in% unique(ctints$overview$celltypes)])
  
  ctints$overview$age <- factor(ctints$overview$age,
                                 levels = names(col_age))
  return(ctints)
})

#-------------------------------------------------------------------------------

# calculate nr of ligands/receptors = nr of lrs per cell type

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

# condition
for(i in 1:length(nrlrs_list)){
  if(grepl("raw", names(nrlrs_list)[i])){
    nrlrs_list[[i]]$condition <- "no pre-processing steps"
  }else if(grepl("nlv", names(nrlrs_list)[i])){
    nrlrs_list[[i]]$condition <- "no levelling"
  }else if(grepl("nct", names(nrlrs_list)[i])){
    nrlrs_list[[i]]$condition <- "no cut-off"
  }else if(grepl("nrm", names(nrlrs_list)[i])){
    nrlrs_list[[i]]$condition <- "no cell type removal"
  }else if(grepl("rev", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "reverse annotation"
  }else if(grepl("nds", names(ctpints_list)[i])){
    ctpints_list[[i]][[1]]$condition <- "no downsampling"
  }else{
    nrlrs_list[[i]]$condition <- "all pre-processing steps"
  }
}

# factorize
nrlrs_list <- lapply(nrlrs_list, function(nrlrs){
  
  nrlrs$annotation <- factor(nrlrs$annotation, levels = names(col_ann))
  nrlrs$age <- factor(nrlrs$age, levels = names(col_age))
  nrlrs$species <- factor(nrlrs$species, levels = names(col_spc))
  
  return(nrlrs)
})

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# add nr of cells per cell type from SCE objects

if(only_full == TRUE){
  
  # load and add corresponding SCE objects
  mmus_yng_sce_ds_ctrm <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE_ds_ctrm"))
  mmus_old_sce_ds_ctrm <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE_ds_ctrm"))
  
  sce_list <- list(
    "mmus_yng_cpi" = mmus_yng_sce_ds_ctrm,
    "mmus_old_cpi" = mmus_old_sce_ds_ctrm)
  
  names(sce_list) <- names
  
}else{
  
  # downsampled and removed cell types
  mmus_yng_sce_ds_ctrm <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE_ds_ctrm"))
  mmus_old_sce_ds_ctrm <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE_ds_ctrm"))
  
  # downsampled but no removed cell types
  mmus_yng_sce_ds <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE_ds"))
  mmus_old_sce_ds <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE_ds"))
  
  # not downsampled but removed cell types
  mmus_yng_sce_ctrm <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_yng_SCE_ctrm"))
  mmus_old_sce_ctrm <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE_ctrm"))
  
  # original objects
  mmus_yng_sce <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE"))
  mmus_old_sce <- readRDS(file = paste0(
    filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE"))
  
  sce_list <- list(
    "mmus_yng_cpi" = mmus_yng_sce_ds_ctrm,
    "mmus_old_cpi" = mmus_old_sce_ds_ctrm,
    "mmus_yng_cpi_raw" = mmus_yng_sce,
    "mmus_old_cpi_raw" = mmus_old_sce,
    "mmus_yng_cpi_nlv" = mmus_yng_sce_ds_ctrm,
    "mmus_old_cpi_nlv" = mmus_old_sce_ds_ctrm,
    "mmus_yng_cpi_nct" = mmus_yng_sce_ds_ctrm,
    "mmus_old_cpi_nct" = mmus_old_sce_ds_ctrm,
    "mmus_yng_cpi_nrm" = mmus_yng_sce_ds,
    "mmus_old_cpi_nrm" = mmus_old_sce_ds,
    "mmus_yng_cpi_nds" = mmus_yng_sce_ctrm,
    "mmus_old_cpi_nds" = mmus_old_sce_ctrm,
    "mmus_yng_cpi_rev" = mmus_yng_sce_ds_ctrm,
    "mmus_old_cpi_rev" = mmus_old_sce_ds_ctrm
  )
}

# extract nr of cells per cell type
nr_cells <- lapply(sce_list, function(sce){
  
  cts <- levels(sce$Identity)
  
  nr_cells <- vector()
  for(i in 1:length(cts)){
    nr_cells[i] <- length(which(sce$Identity == cts[i]))
  }
  
  names(nr_cells) <- cts
  return(nr_cells)
})

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

# save all of them
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
