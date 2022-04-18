# functions for calculation of cell type interactomes
# for detailed explanations see cpi_calculation.Rmd

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
# FUNDAMENTAL FUNCTIONS (hidden functions)
# required for main functions below
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

add_annotation <- function(sce, annotation){
  
  sce$Annotation <- vector(mode = "character", length = ncol(sce))
  
  for(i in 1:length(annotation)){
    if(annotation[i] == "receiver"){
      sce$Annotation[sce$Identity == names(annotation)[i]] <- "receiver"
    }else if(annotation[i] == "emitter"){
      sce$Annotation[sce$Identity == names(annotation)[i]] <- "emitter"
    }
  }
  return(sce)
}

# adds required annotation to the SCE objects, required for datasheet

#-------------------------------------------------------------------------------

make_datasheet <- function(sce){

  datasheet <- data.frame(
    sample = colnames(sce),
    celltype = colData(sce)$Identity,
    annotation = colData(sce)$Annotation
  )
  return(datasheet)
}

# returns datasheet required for extract_matrix

#-------------------------------------------------------------------------------

expr_cells_perc <- function(sce, cts_list, cts, perm = FALSE, assay = assay){
  # cts defined in interaction_matrix()
  # perm to indicate if used for permutation or not
  
  pct_expr_cells_df <- data.frame(row.names = rownames(assays(sce)[[
    grep(assay, names(assays(sce)))]]))
  
  perc_list <- lapply(cts_list, function(x){
    
    perc <- vector()
    # subset the sce object to calculate ratios for each ct separately
    if(perm == FALSE){
      sce <- sce[,colData(sce)$Identity == x]
    }else if(perm == TRUE){
      sce <- sce[,colData(sce)$Annotation == x]
    }
    
    print(paste("........ ", x, ": ", ncol(sce), sep = "")) #track progress
    for(i in 1:nrow(pct_expr_cells_df)){
      # all cells of celltype x that express gene i
      nr_cell_expr <- length(which(assays(sce)[[
        grep(assay, names(assays(sce)))]][i,] > 0))
      # non-expressing cells are all remaining cells of celltype x
      nr_cell_nnex <- ncol(sce) - nr_cell_expr

      perc[i] <- nr_cell_expr/(nr_cell_expr + nr_cell_nnex)
    }
    return(perc)
  })
  
  for(i in 1:length(cts)){
    pct_expr_cells_df[,i] <- perc_list[[i]]
    colnames(pct_expr_cells_df)[i] <- cts[i]
  }
  return(pct_expr_cells_df)
}

# returns dataframe with genes in rows, cell types in cols, and percentage of 
# cells per cell type that express a certain gene for gene cutoff 

#-------------------------------------------------------------------------------

expr_cells_cutoff <- function(cts_list, pct_expr_cells, cutoff){
  # list of cell types
  # pct_expr_cells = percent of expressing cells df from expr_cells_perc()
  # cutoff of minimum percentage of expressing cells required
  
  expr_genes_pos <- which(
    pct_expr_cells[,colnames(pct_expr_cells) == cts_list] >= cutoff)
  return(expr_genes_pos)
}

# returns for each cell type (= col) a vector of gene positions that are 
# expressed in a given percentage of cells per cell type as defined by cutoff

#-------------------------------------------------------------------------------

# this function was adjusted from Adrien's function interactionmatrix()

extract_matrix <- function(counts, interactions, datasheet, expr_genes_pos){
  # counts matrix from SCE object 
  # interaction database
  # datasheet from make_datasheet()
  # gene positions list per cell types from expr_cells_cutoff()
  
  # prepare cts vector to iterate through celltypes separately by annotation
  cts_emitter <- unfactor(unique(datasheet[datasheet[,3] == "emitter", 2]))
  cts_receiver <- unfactor(unique(datasheet[datasheet[,3] == "receiver", 2]))
  
  # empty dataframe of correct dimensions
  interactionmatPlus <- matrix(ncol = ncol(counts), nrow = nrow(interactions))
  interactionmatPlus <- data.frame(interactionmatPlus)
  
  # use only expressed genes that remain after cutoff
  temp_expr <- expr_genes_pos
  for(i in 1:length(c(cts_emitter, cts_receiver))){
    pos <- expr_genes_pos[[i]]
    temp_expr[[i]] <- counts[pos, colnames(counts) == i]
  }
  
  # iterate through all potential interactions i and celltypes j
  for(i in 1:nrow(interactions)){
    for(j in cts_emitter){
      
      # if ligand of interaction i is found in the expressed genes of ct j
      if(grepl(interactions$Ligand_ENSEMBL[i], c(temp_expr[[j]]))){
        # add to the dataframe
        interactionmatPlus[
          # only into row corresponding to interaction i
          # only into columns that are denoted as cells from celltype j
          i, which(datasheet$celltype == j)] <- counts[
            # only from the one row in counts
            # whose gene ID is in the receptor col of lrdb_part interaction i
            rownames(counts) == interactions[i, 2], 
            which(datasheet$celltype == j)]
        # fill other slots with NAs
      }else{
        interactionmatPlus[i, which(datasheet$celltype == j)] <- NA
      }
    }
    # repeat for receivers
    for(j in cts_receiver){
      if(grepl(interactions$Receptor_ENSEMBL[i], c(temp_expr[[j]]))){
        interactionmatPlus[i, which(datasheet$celltype == j)] <- counts[
          rownames(counts) == interactions[i,1], which(datasheet$celltype == j)]
      }else{
        interactionmatPlus[i, which(datasheet$celltype == j)] <- NA
      }
    }
    
    if(i == 1000){
      print("................ 1000 done") #track progress 
    }else if(i == 2000){
      print("....................... 2000 done")
    }
  }
  # remove potential negative values
  interactionmatPlus[interactionmatPlus < 0] <- NA
  return(interactionmatPlus)
}

# returns an interaction matrix with cols = cells and rows = interactions
# extracts counts for ligand genes for emitters and receptor genes for receivers
# extracts counts only for expressed genes that remained after cutoff


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------

# MAIN FUNCTIONS (as seen in cpi_calculation.R)

#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# this function requires other functions:
# - add_annotation()
# - make_datasheet()
# - expr_cells_perc()
# - expr_cells_cutoff()
# - extract_matrix()

interaction_matrix <- function(sce, cutoff, annotation, ct_remove = NULL,
                               assay = NULL){
  # pre-processed single-cell experiment object with annotations
  # a cutoff between 0 and 1 to remove cells expressed in few cells/cell type
  # a vector of cell types that should be removed from the data
  
  # set logcounts as default if no assay type is given
  if(is.null(assay) == TRUE){
    assay <- "logcounts"
  }
  print(assay)

  # subset SCE objects 
  print(".... Preparing objects.") #############################################
  
  # extract and check cell types 
  cts1 <- unique(unfactor(colData(sce)$Identity))
  cts2 <- levels(colData(sce)$Identity)
  cts3 <- names(annotation)
  
  if(!setequal(cts1, cts2)){
    stop("Cell type factors do not correspond to cell types")
  }
  if(length(grep("FALSE", cts1 %in% cts3)) > 1){
    stop("Cell type annotation is missing")
  }
  
  cts <- cts2
    
  #subset sce to remove cell types
  if(is.null(ct_remove) == FALSE){
    
    ct_remove_pos <- list()
    for(i in 1:length(ct_remove)){
      ct_remove_pos[[i]] <- which(sce$Identity == ct_remove[i])
    }
    if(length(ct_remove_pos[[i]]) > 1){
      sce <- sce[,-unlist(ct_remove_pos)]
      # remove empty elements from cts not found in Identity afterward subsetting 
      empty_pos <- which(is.na(match(levels(sce$Identity), unique(sce$Identity))))
      cts <- cts[-empty_pos]
    }
  }
  
  # subset sce by genes from lrdb_part to reduce computing time
  lrdb_ids <- unique(c(unique(lrdb_part$Receptor_ENSEMBL), 
                       unique(lrdb_part$Ligand_ENSEMBL)))
  lrdb_ids <- lrdb_ids[!is.na(match(lrdb_ids, rowData(sce)$ENSMUS_ID))]
  sce <- sce[!is.na(match(rowData(sce)$ENSMUS_ID, lrdb_ids)),]
  
  # create countsmatrix and datasheet for extract_matrix() 
  print(".... Extracting data.") ###############################################

  countsmatrix <- assays(sce)[[grep(assay, names(assays(sce)))]]
  rownames(countsmatrix) <- rowData(sce)$ENSMUS_ID
  
  sce <- add_annotation(sce, annotation = annotation)
  datasheet <- make_datasheet(sce) 
  
  # extract percentage of expressing cells and implement cutoff ################
  print(".... Implementing gene cut off by percentage of expressing cells.")
  
  cts_list <- as.list(cts)
  pct_expr_cells <- expr_cells_perc(sce = sce, 
                                    cts = cts,
                                    cts_list = cts_list,
                                    assay = assay) 
  # get positions of filtered genes per ct
  expr_genes_pos <- lapply(cts_list, 
                           pct_expr_cells = pct_expr_cells, 
                           cutoff = cutoff, 
                           FUN = expr_cells_cutoff) 
  names(expr_genes_pos) <- cts
  
  # create interaction matrix 
  print(".... Generating interaction matrix.") #################################
  
  sce_im <- extract_matrix(counts = countsmatrix, 
                           datasheet = datasheet,
                           expr_genes_pos = expr_genes_pos, 
                           interactions = lrdb_part) 
  
  # return the relevant data
  interactionmatrix_list <- list(
    interactionmatrix = sce_im,
    datasheet = datasheet # keep datasheet for later functions
  )
  return(interactionmatrix_list)
}

# performs subsetting and filtering steps of a pre-processed sce object
# returns an interaction matrix with cols = cells and rows = interactions


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


# this function was adapted from Adrien's interactionranking() function

interaction_ranking <- function(interactionmatrix_list, 
                                cutoff = TRUE,
                                level = TRUE, 
                                top_level){
  # interactionmatrix_list from interaction_matrix
  # cutoff as indicator whether cutoff was implemented
  # level to indicate if ranks should be leveled to a top_level 
  # top_level = a value for leveling that must be higher than the highest rank

  # Preparation ################################################################
  
  # extract data from list
  interactionmatrix <- interactionmatrix_list[[1]]
  datasheet <- interactionmatrix_list[[2]]
  
  # reformat datasheet into two
  datasheetemi <- datasheet[datasheet[,3] %in% "emitter",]
  datasheetrec <- datasheet[datasheet[,3] %in% "receiver",]
  
  # empty objects 
  name <- vector()
  mat <- matrix(
    ncol = length(unique(datasheetemi[,2])) * length(unique(datasheetrec[,2])),
    nrow = nrow(interactionmatrix))
  Score <- mat
  rankLigs <- mat
  rankRecs <- mat
  
  # 1 will be added during each iteration 
  l <- 1
  
  # rank the ligand and receptor genes #########################################
  
  # for each ligand i and for each receptor j
  for(i in 1:length(unique(datasheetemi[,2]))){
    for (j in 1:length(unique(datasheetrec[,2]))){
      # get all positions of the current ligand i 
      pos_lig <- which(datasheet[,2] %in% unique(datasheetemi[,2])[i])
      # rank the mean expression values for those positions
      if(length(pos_lig) > 1){
        rankLig <- dense_rank(rowMeans(
          interactionmatrix[,pos_lig]))
      }else if(length(pos_lig) == 1){
        # if there is only one cell = one column, rowMeans does not work
        rankLig <- dense_rank(interactionmatrix[,pos_lig])
      }
      
      # same for receptors
      pos_rec <- which(datasheet[,2] %in% unique(datasheetrec[,2])[j])
      if(length(pos_rec) > 1){
        rankRec <- dense_rank(rowMeans(
          interactionmatrix[,pos_rec]))
      }else if(length(pos_rec) == 1){
        rankRec <- dense_rank(interactionmatrix[,pos_rec])
      }
      
      # level the ranks ########################################################
      
      # level ranks by adding up to top_level to reduce influences on score
      if(level == TRUE && cutoff == FALSE){
        # a rank of 1 indicates 0 expression (no cut off)

        # determine the highest rank
        max_ints_lig <- max(rankLig, na.rm = TRUE)
        max_ints_rec <- max(rankRec, na.rm = TRUE)
        # determine the resulting lowest rank
        min_rank_lig <- top_level-max_ints_lig
        min_rank_rec <- top_level-max_ints_rec
        # stop the function if the resulting lowest rank is <0
        if(min_rank_lig < 0 | min_rank_rec < 0){
          stop("top_level is too low, set top_level higher")
        }
        
        # add resulting lowest rank to remaining ranks to match the highest rank 
        # but only if ranks are not 1, which indicates expression of 0
        rankLig[rankLig != 1] <- rankLig[rankLig != 1]+min_rank_lig
        rankRec[rankRec != 1] <- rankRec[rankRec != 1]+min_rank_rec
        # remove 1s to avoid tiny scores after score normalization later
        rankLig[rankLig == 1] <- 0
        rankRec[rankRec == 1] <- 0      
        
      }else if(level == TRUE && cutoff == TRUE){
        # a rank of 1 indicates lowest expressed rank (no non-expressing cells)

        max_ints_lig <- max(rankLig, na.rm = TRUE)
        max_ints_rec <- max(rankRec, na.rm = TRUE)
        min_rank_lig <- top_level-max_ints_lig
        min_rank_rec <- top_level-max_ints_rec
        
        if(min_rank_lig < 0 | min_rank_rec < 0){
          stop("max rank level is too low, set max rank level higher")
        }
        
        # add lowest rank only if ranks are not NA, indicating expression of 0
        rankLig[!is.na(rankLig)] <- rankLig[!is.na(rankLig)]+min_rank_lig
        rankRec[!is.na(rankRec)] <- rankRec[!is.na(rankRec)]+min_rank_rec
      }
    
      # put in empty mat, with col corresponding to the number of iteration
      rankLigs[,l] <- rankLig
      rankRecs[,l] <- rankRec

      
      # calculate Rank Sums, Score and match cell types ########################
      
      Score[,l] <- rankLig + rankRec
      
      # store the colnames for Score for later
      name[l] <- paste0(unique(datasheetemi[,2])[i], "&",
                        unique(datasheetrec[,2])[j])
      l=l+1
    }
  }
  
  # Score normalization for each pair of cells 
  for (i in 1:ncol(Score)){
    Score[,i] <- Score[,i]/max(Score[,i], na.rm = TRUE)
  }
  colnames(Score) <- name
  
  interactionranking_list <- list()
  interactionranking_list[[1]] <- Score
  interactionranking_list[[2]] <- rankRecs
  interactionranking_list[[3]] <- rankLigs
  interactionranking_list[[4]] <- datasheet  # keep datasheet
  return(interactionranking_list)
}

# interaction_ranking
# returns a list of three matrices containing for all ctp and interactions
# - interaction score 
# - receptor ranks 
# - ligand ranks 
# - and the datasheet


#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------


interaction_list <- function(interactionranking_list, lrdb_comp){
  # interactionranking_list from interaction_ranking
  # lrdb_comp containing more info on all interactions
  
  # extract info 
  ir_matrix <- list(interactionranking_list[[1]], # Score
                    interactionranking_list[[2]], # Receptorranks
                    interactionranking_list[[3]]) # Ligandranks
  datasheet <- interactionranking_list[[4]]
  lrdb <- lrdb_comp
  
  # add names 
  names(ir_matrix) <- c("Score", "Receptorrank", "Ligandrank")
  rownames(ir_matrix$Score) <- lrdb$interactions$interaction_pair
  rownames(ir_matrix$Receptorrank) <- lrdb$interactions$interaction_pair 
  rownames(ir_matrix$Ligandrank) <-lrdb$interactions$interaction_pair
  colnames(ir_matrix$Receptorrank) <- colnames(ir_matrix$Score)
  colnames(ir_matrix$Ligandrank) <- colnames(ir_matrix$Score)
  
  # coerce into dataframes
  ir_matrix$Score <- as.data.frame(ir_matrix$Score)
  ir_matrix$Receptorrank <- as.data.frame(ir_matrix$Receptorrank)
  ir_matrix$Ligandrank <- as.data.frame(ir_matrix$Ligandrank)
  
  # add colData-like slot = Celltypes ##########################################
  ir_matrix$Celltypes <- data.frame(
    ct_pair = colnames(ir_matrix$Score),
    emitter = vector(mode = "character", length = ncol(ir_matrix$Score)),
    receiver = vector(mode = "character", length = ncol(ir_matrix$Score))
  )
  
  # add annotation as emitter or receiver for better overview
  ct <- unique(datasheet$celltype)
  
  for(i in ct){
    if(datasheet$annotation[datasheet$celltype == i][1] == "emitter"){
      ir_matrix$Celltypes$emitter[grep(i, ir_matrix$Celltypes$ct_pair)] <- i
    }else if(datasheet$annotation[datasheet$celltype == i][1] == "receiver"){
      ir_matrix$Celltypes$receiver[grep(i, ir_matrix$Celltypes$ct_pair)] <- i
    }
  }
  
  # add rowData-like slot = Interactions #######################################
  ir_matrix$Interactions <- data.frame(
    interaction_pair = lrdb$interactions$interaction_pair,
    ligand_symbol = lrdb$interactions$ligand_symbol,
    ligand_ensembl_id = lrdb$interactions$ligand_ensembl_id,
    receptor_symbol = lrdb$interactions$receptor_symbol,
    receptor_ensembl_id = lrdb$interactions$receptor_ensembl_id,
    interaction_type = lrdb$interactions$interaction_type
  )
  return(ir_matrix)
}

# transforms a list from interaction_ranking() into into another list
# with more information 
