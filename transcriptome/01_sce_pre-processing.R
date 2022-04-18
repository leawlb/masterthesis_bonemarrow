#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# This script uses pre-processed SCE objects per species with:
# -cell type annotations in "Identity" colData slot
# -reduced dimensions
# -high-quality cells
# -fastMNN batch correction

# and further prepares them by
# -factorising cell types
# -adding ENSMUS IDs
# -subsetting into objects per age group for separate analysis

#-------------------------------------------------------------------------------

# load libraries, source code
library(SingleCellExperiment)
library(rlist)

#setwd("/Users/leaw/Documents/repositories/Interspecies_BM/new_Scripts")
setwd("/home/l012t/Interspecies_BM")

source(file = "code/filepath.R")
source(file = "code/colors.R")
source(file = "transcriptome/functions/sce_functions.R")

#-------------------------------------------------------------------------------

# load insp_ref required for get_shared_ensids()
insp_ref <- list.load(file = paste0(
  filepath_lv1, "data/databases/insp_ref_shared_man.rds"))

#-------------------------------------------------------------------------------

# SCE objects
mmus_all <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_all_fastMNN_annotated"))

mcar_all <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mcar/mcar_all_fastMNN_annotated"))

mspr_all <- readRDS(file = paste0(
  filepath_lv1, "data/sce_objects/mspr/mspr_all_fastMNN_annotated"))

#mcas_all <- readRDS(file = paste0(
#  filepath_lv1, "data/sce_objects/mcas/mcas_all_fastMNN_annotated"))

#-------------------------------------------------------------------------------

# subset by age
mmus_yng <- mmus_all[,mmus_all$Age == "young"]
mmus_old <- mmus_all[,mmus_all$Age == "old"]

mcar_yng <- mcar_all[,mcar_all$Age == "young"]
mcar_old <- mcar_all[,mcar_all$Age == "old"]

mspr_yng <- mspr_all[,mspr_all$Age == "young"]
mspr_old <- mspr_all[,mspr_all$Age == "old"]

#mcas_yng <- mcas_all[,mcas_all$Age == "young"]
#mcas_old <- mcas_all[,mcas_all$Age == "old"]

#-------------------------------------------------------------------------------

# put SCE objects in list
sce_list <- list(
  mmus_yng,
  mmus_old,
  mcar_yng,
  mcar_old,
  mspr_yng,
  mspr_old
  #,mcas_yng,
  #mcas_old
)

#-------------------------------------------------------------------------------

# pre-process SCE objects

# this function adds info on transcripts shared between species and ensmus IDs 
sce_list <- lapply(sce_list, get_shared_ensids)

# this function changes ct names and factorises them for better visualisation
sce_list <- lapply(sce_list, change_cts)

#-------------------------------------------------------------------------------

# save the subset objects

# BL6
saveRDS(sce_list[[1]], file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_young_SCE"))
saveRDS(sce_list[[2]], file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_old_SCE"))

# CARO
saveRDS(sce_list[[3]], file = paste0(
  filepath_lv1, "data/sce_objects/mcar/mcar_young_SCE"))
saveRDS(sce_list[[4]], file = paste0(
  filepath_lv1, "data/sce_objects/mcar/mcar_old_SCE"))

# SPRET
saveRDS(sce_list[[5]], file = paste0(
  filepath_lv1, "data/sce_objects/mspr/mspr_young_SCE"))
saveRDS(sce_list[[6]], file = paste0(
  filepath_lv1, "data/sce_objects/mspr/mspr_old_SCE"))

# CAST
#saveRDS(sce_list_merged[[7]], file = paste0(
#  filepath_lv1, "data/sce_objects/mcas/mcas_young_SCE"))
#saveRDS(sce_list_merged[[8]], file = paste0(
#  filepath_lv1, "data/sce_objects/mcas/mcas_old_SCE"))

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# also pre-process and save the whole objects

sce_list_all <- list(
  mmus_all,
  mcar_all,
  mspr_all
  #,mcas_all
)

sce_list_all <- lapply(sce_list_all, get_shared_ensids)
sce_list_all <- lapply(sce_list_all, change_cts)

saveRDS(sce_list_all[[1]], file = paste0(
  filepath_lv1, "data/sce_objects/mmus/mmus_SCE"))

saveRDS(sce_list_all[[2]], file = paste0(
  filepath_lv1, "data/sce_objects/mcar/mcar_SCE"))

saveRDS(sce_list_all[[3]], file = paste0(
  filepath_lv1, "data/sce_objects/mspr/mspr_SCE"))
 
#saveRDS(sce_list_merged[[4]], file = paste0(
#  filepath_lv1, "data/sce_objects/mcas/mcas_SCE"))
