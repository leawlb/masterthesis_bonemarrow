#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# this script calculates cell type pair interactions (CPI) from SCE objects

# running all of it at once takes very long and can cause memory problems

# SCE objects should have 
# -ct annotation in an colData row calles "Identity"
# -ct names that cannot grep each other and do not contain x, (, ), or &

# A ligand receptor database of suitable format is required

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# PRE-PROCESSING
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# load libraries, source code
library(SingleCellExperiment)
library(tidyverse)
library(rlist)
library(scuttle)

setwd("/home/l012t/Interspecies_BM")

source(file = "code/filepath.R")
source(file = "interactome/functions/cpi_functions_calculation.R")
# required for cts_ordered
source(file = "transcriptome/functions/sce_functions.R")

#-------------------------------------------------------------------------------

# load LRDB
lrdb_part <- list.load(file = paste(
  filepath_lv1, "data/databases/lrdb_part.rds", sep = ""))

lrdb_comp <- list.load(file = paste(
  filepath_lv1, "data/databases/lrdb_comp.rds", sep = ""))

#-------------------------------------------------------------------------------
# load one SCE object per condition to investigate

# BL6
mmus_yng_sce <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE"))
mmus_old_sce <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE"))

# CARO
#mcar_yng_sce <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mcar/mcar_young_SCE"))
#mcar_old_sce <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mcar/mcar_old_SCE"))

# SPRET
#mspr_yng_sce <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mspr/mspr_young_SCE"))
#mspr_old_sce <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mspr/mspr_old_SCE"))

# CAST
#mcas_yng_sce <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mcas/mcas_young_SCE"))
#mcas_old_sce <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mcas/mcas_old_SCE"))

#-------------------------------------------------------------------------------
# put SCE objects in lists 

sce_list <- list(
  mmus_yng_sce,
  mmus_old_sce#,
  #mcar_yng_sce,
  #mcar_old_sce,
  #mspr_yng_sce,
  #mspr_old_sce
  #,mcas_yng_sce,
  #mcas_old_sce
)

names <- c(
  "mmus_yng",
  "mmus_old"#,
  #"mcar_yng",
  #"mcar_old",
  #"mspr_yng",
  #"mspr_old"
  #,"mcas_yng",
  #"mcas_old"
)

species <- c(
  "mmus"#,
  #"mcar",
  #"mspr"#,
  #"mcas"
)

names(sce_list) <- names

#-------------------------------------------------------------------------------
# cell type vectors

# cell type vector of possible hematopoietic contaminants in stromal fraction
# these will be removed ONLY from the stromal fraction (optional)
hem_cts <- c("LMPPs", "Small pre-BCs", "Large pre-BCs", "Pro-BCs",
             "Mature B cells", "T cells", "NK cells", "Gran/Mono progs",
             "Monocyte progs", "Neutrophil progs", "Eo/Baso progs",
             "Monocytes", "Neutrophils", "Dendritic cells", "Ery/Mk progs",
             "Megakaryocyte progs", "Erythrocyte progs", "Erythroids",
             "Baso/Eo/Mast progs", "T cells/NK cells", "Monocyte-primed progs",
             "Dendritic-primed progs", "LT-HSCs", "Common myeloid progs",
             "Gran/Mono progs"
)

# cell type vector of cell types to remove from ALL fractions (optional)
# these are low-abundant cts or cts that are not interesting to the analysis
remove_cts <- c("Schwann cells", "Mature B cells", "T cells/NK cells", 
                "Pro-BCs", "Dendritic cells", "Monocytes", "Neutrophils", 
                "Erythroids", "Small pre-BCs", "Large pre-BCs", 
                "Immature B cells",  "Chondrocytes", "Myofibroblasts", 
                "Gran/Mono progs")

# a vector of ct names with corresponding annotation of emis or recs
# (REQUIRED)
cts_annotated <- c(
  "LT-HSCs" = "receiver",
  "ST-HSCs/MPPs" = "receiver",
  
  "LMPPs" = "receiver",
  
  "Small pre-BCs" = "receiver", 
  "Large pre-BCs" = "receiver",
  "Pro-BCs" = "receiver",
  "Immature B cells" = "receiver",
  "Mature B cells" = "receiver",
  "T cells/NK cells" = "receiver",
  
  "Common myeloid progs" = "receiver",
  "Gran/Mono progs" = "receiver",
  "Monocyte-primed progs" = "receiver",  
  "Monocyte progs" = "receiver",
  "Dendritic-primed progs" = "receiver",  
  "Neutrophil progs" = "receiver",
  "Baso/Eo/Mast progs" = "receiver",
  
  "Monocytes" = "receiver",
  "Dendritic cells" = "receiver",
  "Neutrophils" = "receiver",
  
  "Ery/Mk progs" = "receiver",
  "Megakaryocyte progs" = "receiver",
  "Erythrocyte progs" = "receiver",
  "Erythroids" = "receiver",
  
  "Ng2-MSCs" = "emitter",
  "Fibro/Chondro progs" = "emitter",
  "Myofibroblasts" = "emitter",
  "Arteriolar fibros" = "emitter",
  "Endosteal fibros" = "emitter",
  "Stromal fibros" = "emitter",
  
  "Adipo-CARs" = "emitter",
  "Osteo-CARs" = "emitter",
  
  "Arteriolar ECs" = "emitter",
  "Sinusoidal ECs" = "emitter",
  "Smooth muscle" = "emitter",
  
  "Chondrocytes" = "emitter",
  "Osteoblasts" = "emitter",
  
  "Schwann cells" = "emitter"
)

#-------------------------------------------------------------------------------
# remove hematopoietic contamination from stromal cells (optional)

sce_list <- lapply(sce_list, function(sce){
  
  # remove the hem_cts ONLY from cells that are also stromal
  hem_pos <- intersect(which(sce$Fraction == "stromal"),
                       which(sce$Identity %in% hem_cts))
  
  sce <- sce[,-hem_pos] 
  
  # print to be sure it's correct
  print(unique(sce$Identity[sce$Fraction == "stromal"]))
  #print(unique(sce$Identity))
  
  return(sce)
  
})

#-------------------------------------------------------------------------------
# adjust SCE objects 

sce_list_ds_ctrm <- sce_list
sce_list_ds <- sce_list
sce_list_ctrm <- sce_list

names(sce_list_ds_ctrm) <- names
names(sce_list_ds) <- names
names(sce_list_ctrm) <- names

# for each species, combine the young and old objects to sample down,
# then separate and save again in the same spot in the list

# sample down and remove cell types
for(i in species){
  
  # combine the young and old object of each species
  sce_yng_temp <- sce_list_ds_ctrm[[intersect(grep(species, 
                                                   names(sce_list_ds_ctrm)),
                                      grep("yng", names(sce_list_ds_ctrm)))]]
  sce_old_temp <- sce_list_ds_ctrm[[intersect(grep(species, 
                                                   names(sce_list_ds_ctrm)),
                                      grep("old", names(sce_list_ds_ctrm)))]]
  
  # combine and remove cell types
  sce_temp <- cbind(sce_yng_temp, sce_old_temp)
  sce_temp <- sce_temp[,-which(sce_temp$Identity %in% remove_cts)]
  
  # downsample
  assays(sce_temp)$downsampled <- downsampleBatches(assays(sce_temp)$counts, 
                                                    batch = sce_temp$Age)
  
  # separate again
  sce_temp_yng <- sce_temp[,sce_temp$Age == "young"]
  sce_temp_old <- sce_temp[,sce_temp$Age == "old"]
  
  # save in the same positions in the list the objects were taken out
  sce_list_ds_ctrm[[intersect(grep(species, names(sce_list_ds_ctrm)),
                      grep("yng", names(sce_list_ds_ctrm)))]] <- sce_temp_yng
  sce_list_ds_ctrm[[intersect(grep(species, names(sce_list_ds_ctrm)),
                      grep("old", names(sce_list_ds_ctrm)))]] <- sce_temp_old
}

# sample down but keep all cell types
for(i in species){
  
  sce_yng_temp <- sce_list_ds[[intersect(grep(species, names(sce_list_ds)),
                                      grep("yng", names(sce_list_ds)))]]
  sce_old_temp <- sce_list_ds[[intersect(grep(species, names(sce_list_ds)),
                                      grep("old", names(sce_list_ds)))]]
  
  sce_temp <- cbind(sce_yng_temp, sce_old_temp)
  assays(sce_temp)$downsampled <- downsampleBatches(assays(sce_temp)$counts, 
                                                    batch = sce_temp$Age)
  
  sce_temp_yng <- sce_temp[,sce_temp$Age == "young"]
  sce_temp_old <- sce_temp[,sce_temp$Age == "old"]
  
  sce_list_ds[[intersect(grep(species, names(sce_list_ds)),
                              grep("yng", names(sce_list_ds)))]] <- sce_temp_yng
  sce_list_ds[[intersect(grep(species, names(sce_list_ds)),
                              grep("old", names(sce_list_ds)))]] <- sce_temp_old
}

# no downsampling, but cell type removal
for(i in 1:length(sce_list_ctrm)){
  
  sce_list_ctrm[[i]] <- sce_list_ctrm[[i]][
    ,-which(sce_list_ctrm[[i]]$Identity %in% remove_cts)]
}

#-------------------------------------------------------------------------------
# remove empty factors

for(i in 1:length(sce_list_ds_ctrm)){
  sce_list_ds_ctrm[[i]]$Identity <- factor(
    sce_list_ds_ctrm[[i]]$Identity, 
    levels = cts_ordered[cts_ordered %in%
                           unique(sce_list_ds_ctrm[[i]]$Identity)])
}

for(i in 1:length(sce_list_ds)){
  sce_list_ds[[i]]$Identity <- factor(
    sce_list_ds[[i]]$Identity, 
    levels = cts_ordered[cts_ordered %in% unique(sce_list_ds[[i]]$Identity)])
}

for(i in 1:length(sce_list_ctrm)){
  sce_list_ctrm[[i]]$Identity <- factor(
    sce_list_ctrm[[i]]$Identity, 
    levels = cts_ordered[cts_ordered %in% unique(sce_list_ctrm[[i]]$Identity)])
}

for(i in 1:length(sce_list)){ # still important after hem_ct removal
  sce_list[[i]]$Identity <- factor(
    sce_list[[i]]$Identity, 
    levels = cts_ordered[cts_ordered %in% unique(sce_list[[i]]$Identity)])
}

# save all objects for later use 
for(i in 1:length(sce_list_ds_ctrm)){
  species <- sce_list_ds_ctrm[[i]]$Species[1]
  
  saveRDS(sce_list_ds_ctrm[[i]], file = paste0(
    filepath_lv1, "data/sce_objects/", species, "/", names[i], "_SCE_ds_ctrm"))
}

for(i in 1:length(sce_list_ds)){
  species <- sce_list_ds[[i]]$Species[1]
  
  saveRDS(sce_list_ds[[i]], file = paste0(
    filepath_lv1, "data/sce_objects/", species, "/", names[i], "_SCE_ds"))
}

for(i in 1:length(sce_list_ctrm)){
  species <- sce_list_ctrm[[i]]$Species[1]
  
  saveRDS(sce_list_ctrm[[i]], file = paste0(
    filepath_lv1, "data/sce_objects/", species, "/", names[i], "_SCE_ctrm"))
}

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# CPI CALCULATIONS
# this requires multiple functions from cpi_functions_calc.R
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# FULL CALCULATION FULL OBJECTS for analysis

cpi_list <- lapply(sce_list_ds_ctrm, function(sce){
  
  # track progress
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  # extracts an interaction matrix (im) of all LR genes and cells
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = remove_cts, 
                               annotation = cts_annotated,
                               assay = "downsampled")
  
  # interaction ranking (ir) and scoring of extracted LR genes and cells
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = TRUE,
                                top_level = 300)
  
  # transforms the ranked object into list of cell type pair interactomes (cpi)
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  # adds information on the CPI calc process, required for analysis functions
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = TRUE)
  
  # add metadata directly
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "all pre-processing steps"
  
  return(cpi_obj)
})

names(cpi_list) <- names

# order all objects alphabetically for direct comparisons
cpi_list <- lapply(cpi_list, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

# save the cpi objects 
for(i in 1:length(cpi_list)){
  
  saveRDS(cpi_list[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list)[i], "_cpi"))
}

# remove to make space for long run
rm(cpi_list)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# CPI CALCULATION SAMPLES (BATCHES)
# calculate cpi for separate batches with both fractions

# obtain all samples as separate objects 
sce_sample_list_batches <- lapply(sce_list_ds_ctrm, function(sce){
  
  # subset all SCE objects by Batch_exp and add to a nested list
  samples <- list()
  for(i in 1:length(unique(sce$Batch_exp))){
    samples[[i]] <- sce[,sce$Batch_exp == unique(sce$Batch_exp)[i]]
  }
  return(samples)
})

# "unnest" the list to obtain one entry per sample
sce_sample_list_batches <- unlist(sce_sample_list_batches)

# add the names of the samples for saving later on
names_keep <- vector()
j <- 1
for(i in 1:length(sce_sample_list_batches)){
  
  # for each Batch_exp, make unique sample name from the "Batch_exp" slot
  names(sce_sample_list_batches)[i] <- paste(
    sce_sample_list_batches[[i]]$Species[1], 
    sce_sample_list_batches[[i]]$Age[1],
    sce_sample_list_batches[[i]]$Batch_exp[1],  
    sep = "_")
  
  # to only use sce object batches with both str and hsc fractions
  if(length(unique(sce_sample_list_batches[[i]]$Fraction)) == 2){
    names_keep[j] <- names(sce_sample_list_batches)[i]
    j <- j +1
  }
}

# only keep the samples with both fractions
sce_sample_list_batches <- unlist(sce_sample_list_batches)[names_keep]
# control that all elements have both fractions
for(i in 1:length(sce_sample_list_batches)){
  print(unique(sce_sample_list_batches[[i]]$Fraction))
}

#-------------------------------------------------------------------------------

# same conditions as full object calculation above
cpi_sample_list_batches <- lapply(sce_sample_list_batches, function(sce){
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = remove_cts,
                               annotation = cts_annotated,
                               assay = "downsampled")
  
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = TRUE,
                                top_level = 300)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = TRUE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "all pre-processing steps"
  
  return(cpi_obj)
})

names(cpi_sample_list_batches) <- paste0(names(sce_sample_list_batches), "_cpi")

cpi_sample_list_batches <- lapply(cpi_sample_list_batches, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

# save the cpi objects
for(i in 1:length(cpi_sample_list_batches)){
  
  saveRDS(cpi_sample_list_batches[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_sample_list_batches)[i]))
}

rm(cpi_sample_list_batches)
rm(sce_sample_list_batches)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# CPI CALCULATION SAMPLES (REPLICATES)
# calculate cpi for separate replicates with both fractions

# obtain all samples as separate objects 
sce_sample_list_reps <- lapply(sce_list_ds_ctrm, function(sce){
  
  # subset all SCE objects by Sample 
  samples <- list()
  for(i in 1:length(unique(sce$Sample))){
    samples[[i]] <- sce[,sce$Sample == unique(sce$Sample)[i]]
  }
  return(samples)
})

sce_sample_list_reps <- unlist(sce_sample_list_reps)

for(i in 1:length(sce_sample_list_reps)){
  names(sce_sample_list_reps)[i] <- sce_sample_list_reps[[i]]$Sample[1]
}

#combine fractions (hsc and stromal) for each replicate
sce_sample_list_reps_merged <- list()
j <- 1
for(i in 1:length(sce_sample_list_reps)){
  
  # take each "hsc" object and if there is a corresponding "stromal" object,
  # add it
  if(grepl("hsc", names(sce_sample_list_reps)[i])){
    hsc <- sce_sample_list_reps[[i]]
    stromal_pos <- grep(gsub("_hsc","_stromal",names(sce_sample_list_reps)[i]), 
                        names(sce_sample_list_reps))
    if(length(stromal_pos) > 0){
      stromal <- sce_sample_list_reps[[stromal_pos]]
      merged <- cbind(hsc, stromal)
      sce_sample_list_reps_merged[[j]] <- merged
      j <- j + 1
    }
  }
}

# add names to merged
for(i in 1:length(sce_sample_list_reps_merged)){
  names(sce_sample_list_reps_merged)[i] <-  gsub(
    "_hsc", "", sce_sample_list_reps_merged[[i]]$Sample[1])
}

for(i in 1:length(sce_sample_list_reps_merged)){
  print(unique(sce_sample_list_reps_merged[[i]]$Fraction))
}

#-------------------------------------------------------------------------------

# same conditions as full object calculation above
cpi_sample_list_reps <- lapply(sce_sample_list_reps_merged, function(sce){

  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = remove_cts,
                               annotation = cts_annotated,
                               assay = "downsampled")
  
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = TRUE,
                                top_level = 300)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = TRUE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "all pre-processing steps"
  
  return(cpi_obj)
})

names(cpi_sample_list_reps) <- paste0(names(cpi_sample_list_reps), "_cpi")

cpi_sample_list_reps <- lapply(cpi_sample_list_reps, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

# save the cpi objects
for(i in 1:length(cpi_sample_list_reps)){
  
  saveRDS(cpi_sample_list_reps[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_sample_list_reps)[i]))
}

rm(cpi_sample_list_reps)
rm(sce_sample_list_reps_merged)
rm(sce_sample_list_reps)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# CONDITION CALCULATION 
# calculate CPI at different conditions

#-------------------------------------------------------------------------------

# raw (= no conditions, cutoff = 0, no cell types removed)
cpi_list_raw <- lapply(sce_list, function(sce){
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0, 
                               ct_remove = NULL, 
                               annotation = cts_annotated,
                               assay = "logcounts")
  
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = FALSE, 
                                level = FALSE)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = FALSE, 
                          level = FALSE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "no pre-processing steps"
  
  return(cpi_obj)
})

names(cpi_list_raw) <- paste0(names, "_cpi_raw")

cpi_list_raw <- lapply(cpi_list_raw, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

for(i in 1:length(cpi_list_raw)){
  saveRDS(cpi_list_raw[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_raw)[i]))
}

rm(cpi_list_raw)

#-------------------------------------------------------------------------------

# no levelling
cpi_list_nlv <- lapply(sce_list_ds_ctrm, function(sce){ 
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = remove_cts,
                               annotation = cts_annotated,
                               assay = "downsampled")
  
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = FALSE)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = FALSE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "no levelling"
  
  return(cpi_obj)
})

names(cpi_list_nlv) <- paste0(names, "_cpi_nlv")

cpi_list_nlv <- lapply(cpi_list_nlv, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

for(i in 1:length(cpi_list_nlv)){
  saveRDS(cpi_list_nlv[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_nlv)[i]))
}

rm(cpi_list_nlv)

#-------------------------------------------------------------------------------
# no cutoff

cpi_list_nct <- lapply(sce_list_ds_ctrm, function(sce){ 
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0, 
                               ct_remove = remove_cts,
                               annotation = cts_annotated,
                               assay = "downsampled")
  
  ir_obj <- interaction_ranking(im_obj, 
                                cutoff = FALSE, 
                                level = TRUE,
                                top_level = 800) # must be much higher
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = FALSE, 
                          level = TRUE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "no cut-off"
  
  return(cpi_obj)
})

names(cpi_list_nct) <- paste0(names, "_cpi_nct")

cpi_list_nct <- lapply(cpi_list_nct, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

for(i in 1:length(cpi_list_nct)){
  saveRDS(cpi_list_nct[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_nct)[i]))
}

rm(cpi_list_nct)

#-------------------------------------------------------------------------------
# no cell type removal
cpi_list_nrm <- lapply(sce_list_ds, function(sce){
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = NULL, 
                               annotation = cts_annotated,
                               assay = "downsampled")
  
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = TRUE,
                                top_level = 300)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = TRUE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "no cell type removal"
  
  return(cpi_obj)
})

names(cpi_list_nrm) <- paste0(names, "_cpi_nrm")

cpi_list_nrm <- lapply(cpi_list_nrm, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

for(i in 1:length(cpi_list_nrm)){
  saveRDS(cpi_list_nrm[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_nrm)[i]))
}

rm(cpi_list_nrm)

#-------------------------------------------------------------------------------
# no downsampling
cpi_list_nds <- lapply(sce_list_ctrm, function(sce){
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = NULL, 
                               annotation = cts_annotated,
                               assay = "logcounts")
  
  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = TRUE,
                                top_level = 300)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = TRUE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "no downsampling"
  
  return(cpi_obj)
})

names(cpi_list_nds) <- paste0(names, "_cpi_nds")

cpi_list_nds <- lapply(cpi_list_nds, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

for(i in 1:length(cpi_list_nds)){
  saveRDS(cpi_list_nds[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_nds)[i]))
}

rm(cpi_list_nds)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# REVERSE annotation
# calculate CPI with hematopoietic emitters and stromal receivers

#-------------------------------------------------------------------------------

cts_annotated_reverse <- c(
  "LT-HSCs" = "emitter",
  "ST-HSCs/MPPs" = "emitter",
  
  "LMPPs" = "emitter",
  
  "Small pre-BCs" = "emitter", 
  "Large pre-BCs" = "emitter",
  "Pro-BCs" = "emitter",
  "Immature B cells" = "emitter",
  "Mature B cells" = "emitter",
  "T cells/NK cells" = "emitter",
  
  "Common myeloid progs" = "emitter",
  "Gran/Mono progs" = "emitter",
  "Monocyte-primed progs" = "emitter",  
  "Monocyte progs" = "emitter",
  "Dendritic-primed progs" = "emitter",  
  "Neutrophil progs" = "emitter",
  "Baso/Eo/Mast progs" = "emitter",
  
  "Monocytes" = "emitter",
  "Dendritic cells" = "emitter",
  "Neutrophils" = "emitter",
  
  "Ery/Mk progs" = "emitter",
  "Megakaryocyte progs" = "emitter",
  "Erythrocyte progs" = "emitter",
  "Erythroids" = "emitter",
  
  "Ng2-MSCs" = "receiver",
  "Fibro/Chondro progs" = "receiver",
  "Myofibroblasts" = "receiver",
  "Arteriolar fibros" = "receiver",
  "Endosteal fibros" = "receiver",
  "Stromal fibros" = "receiver",
  
  "Adipo-CARs" = "receiver",
  "Osteo-CARs" = "receiver",
  
  "Arteriolar ECs" = "receiver",
  "Sinusoidal ECs" = "receiver",
  "Smooth muscle" = "receiver",
  
  "Chondrocytes" = "receiver",
  "Osteoblasts" = "receiver",
  
  "Schwann cells" = "receiver"
)

#-------------------------------------------------------------------------------

cpi_list_rev <- lapply(sce_list_ds_ctrm, function(sce){
  
  print(paste(sce$Species[1], sce$Age[1], sep = " "))
  
  im_obj <- interaction_matrix(sce, 
                               cutoff = 0.05, 
                               ct_remove = remove_cts, 
                               annotation = cts_annotated_reverse,
                               assay = "downsampled")

  ir_obj <- interaction_ranking(im_obj,
                                cutoff = TRUE, 
                                level = TRUE,
                                top_level = 300)
  
  cpi_obj <- interaction_list(ir_obj, lrdb_comp)
  cpi_obj$Conditions <- c(cutoff = TRUE, 
                          level = TRUE,
                          reverse = TRUE)
  
  cpi_obj$Celltypes$Age <- sce$Age[1]
  cpi_obj$Celltypes$Species <- sce$Species[1]
  cpi_obj$Celltypes$Condition <- "reversed annotation"
  
  return(cpi_obj)
})

names(cpi_list_rev) <- paste0(names, "_cpi_rev")

# order all objects in the same way
cpi_list_rev <- lapply(cpi_list_rev, function(cpi){
  
  alph_order_old <- order(colnames(cpi$Score))
  
  cpi$Score <- cpi$Score[,alph_order_old]
  cpi$Ligandrank <- cpi$Ligandrank[,alph_order_old]
  cpi$Receptorrank <- cpi$Receptorrank[,alph_order_old]
  cpi$Celltypes <- cpi$Celltypes[alph_order_old,]
  
  return(cpi)
})

for(i in 1:length(cpi_list_rev)){
  
  saveRDS(cpi_list_rev[[i]], file = paste0(
    filepath_lv1, "data/cci_files/", names(cpi_list_rev)[i]))
}
