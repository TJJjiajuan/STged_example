
data_process <- function(sc_exp, sc_label, spot_exp, spot_loc,
                         gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                         nUMI = 100, verbose = FALSE, depthscale = 1e5, clean.only = TRUE) {
  # Ensure correct dimensions of inputs
  if(ncol(sc_exp) != length(sc_label)) {
    stop("Mismatch between number of cells in sc_exp and length of sc_label.")
  }
  if(ncol(spot_exp) != nrow(spot_loc)) {
    stop("Mismatch between number of spots in spot_exp and number of rows in spot_loc.")
  }
  
  # Process scRNA-seq data
  sc_matrix <- cleanCounts(sc_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                           expression_threshold = expression_threshold, nUMI = nUMI,
                           verbose = verbose, depthscale = depthscale, clean.only = clean.only)
  sc_matrix <- as.matrix(sc_matrix)
  ind <- match(colnames(sc_matrix), colnames(sc_exp))
  sc_label <- sc_label[ind]
  
  # Process SRT data
  st_matrix <- cleanCounts(spot_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                           expression_threshold = expression_threshold, nUMI = nUMI,
                           verbose = verbose, depthscale = depthscale, clean.only = clean.only)
  st_matrix <- as.matrix(st_matrix)
  ind_sp <- match(colnames(st_matrix), colnames(spot_exp))
  spot_loc <- spot_loc[ind_sp, ]
  
  # Find common genes between datasets
  com_gene <- intersect(rownames(sc_matrix), rownames(st_matrix))
  sc_exp <- sc_matrix[com_gene,]
  st_exp <- st_matrix[com_gene,]
  
  # Re-check and filter by nUMI counts
  index_sc <- colSums(sc_exp) >= nUMI
  sc_exp_filter <- sc_exp[, index_sc]
  sc_label_filter <- sc_label[index_sc]
  
  index_st <- colSums(st_exp) >= nUMI
  st_exp_filter <- st_exp[, index_st]
  spot_loc_filter <- spot_loc[index_st,]
  
  # Compile and return the database of cleaned data
  database <- list(sc_exp = sc_exp_filter, sc_label = sc_label_filter,
                   spot_exp = st_exp_filter, spot_loc = spot_loc_filter)
  return(database)
}


cleanCounts <- function(counts, gene_det_in_min_cells_per = 0.01, expression_threshold = 0,
                        nUMI = 100, verbose = FALSE, depthscale = 1, clean.only = FALSE) {
  
  if (!inherits(counts, "matrix")) stop("The 'counts' must be a matrix.")
  n = ncol(counts)
  if (verbose) message("Starting with ", n, " cells.")
  
  # Selecting genes detected in sufficient number of cells
  sufficient_cells = rowSums(counts > expression_threshold) > gene_det_in_min_cells_per * n
  counts = counts[sufficient_cells, ]
  
  # Filtering cells based on total UMI counts after gene filtering
  sufficient_umi = colSums(counts) > nUMI
  counts = counts[, sufficient_umi]
  
  if (clean.only) {
    if (verbose) message("Matrix filtered without scaling. Dimension: ", dim(counts))
    return(counts)
  } else {
    sf = colSums(counts) 
    scaled_counts = log(sweep(counts, 2, sf, FUN = "/") * depthscale + 1)
    if (verbose) message("Matrix scaled and log-transformed. Dimension: ", dim(scaled_counts))
    return(scaled_counts)
  }
}



source("runindi.R")

library(CARD)
library("spacexr")
library(Seurat)
library(SPOTlight)
library(Giotto)
library(SONAR)

library(here)
library(Matrix)
library(data.table)
library(Seurat)
library(matlabr)
library(R.matlab)


#python_env <- 'C:/Users/visitor01/.conda/envs/ccnet/python.exe'
#reticulate::use_python(python = python_env, required = TRUE)
library(reticulate)
use_condaenv("ccnet", required = TRUE)


## load the 100 um2
#sim_data = readRDS( "C:/Users/visitor01/Desktop/TJJtopic/STged/STged_reservised/STgedNAR/stgedV1/realdata/MAOP/FN7_hash_mpoa_sim-29_ref01.Rds")

sim_data = readRDS("D:/TJJtopic/STged/STged_reservised/STgedNAR/stgedV1/realdata/MAOP/FN7_hash_mpoa_sim-29_ref01.Rds")

# load the referene from scRNA
#ref_data = readRDS("C:/Users/visitor01/Desktop/TJJtopic/STged/STged_reservised/STgedNAR/simulation_V2/GSE113576/scRNA_MPOA.RDS")

ref_data = readRDS("D:/TJJtopic/STged/STged_reservised/STgedNAR/simulation_V2/GSE113576/scRNA_MPOA.RDS")


sc_exp = ref_data$scexp
sc_label = ref_data$sclabel
sc_label <- as.factor(sc_label)

#genes =rownames(sc_exp)
#saveRDS(genes,"merfish_genes.RDS")


## spot data with different resolution
i = 2
patch_size = c(100,50)[i]
mpoa = sim_data$mpoa_list[[i]]$`-0.29`


spot_raw = rownames(mpoa$patchGexp)
spots = paste0("spot",1:length(spot_raw))

split_strings <- strsplit(spot_raw, "_")
spot_loc <- do.call(rbind, lapply(split_strings, function(x) as.numeric(x)))
spot_exp <- t(mpoa$patchGexp)[,spot_raw ]


rownames(spot_loc)  = spots
colnames(spot_exp) = spots


database = data_process(sc_exp = sc_exp, 
                        sc_label = sc_label , 
                        spot_exp = spot_exp,
                        spot_loc = as.data.frame(spot_loc),
                        gene_det_in_min_cells_per = 0, expression_threshold = 0,
                        nUMI = 10, verbose = FALSE, depthscale = 1e6, 
                        clean.only = TRUE)



## use for steroscope and cell2location

spot_exp_clear <- database$spot_exp
spot_loc_clear <- database$spot_loc
ref_exp_clear <- database$sc_exp
ref_label_clear  <- database$sc_label 



if(i==1){
  
  write.csv(spot_exp_clear,"spot_exp_clear_100um.csv")
  write.csv(spot_loc_clear,"spot_loc_clear_100um.csv")
  write.csv(ref_exp_clear,"ref_exp_clear_100um.csv")
  write.csv(ref_label_clear,"ref_label_clear_100um.csv")
  
}else{
  
  write.csv(spot_exp_clear,"spot_exp_clear_50um.csv")
  write.csv(spot_loc_clear,"spot_loc_clear_50um.csv")
  write.csv(ref_exp_clear,"ref_exp_clear_50um.csv")
  write.csv(ref_label_clear,"ref_label_clear_50um.csv")
  
}


##run individual models

SpatialDWLS1 <- suppressWarnings(spatialDWLS_run(database, 
                                                 my_python_path = python_env))

#  saveRDS(SpatialDWLS1,"merfish_SpatialDWLS_100um.RDS")   
saveRDS(SpatialDWLS1,"merfish_SpatialDWLS_50um.RDS")   

RCTD1 <- RCTD_run(database)

#saveRDS(RCTD1,"merfish_RCTD_100um.RDS")    
saveRDS(RCTD1,"merfish_RCTD_50um.RDS") 


CARD1 <- CARD_run(database)

#saveRDS(CARD1,"merfish_CARD_100um.RDS") 
saveRDS(CARD1,"merfish_CARD_50um.RDS")


SPOTlight1 <- suppressWarnings(SPOTlight_run(database))

# For R code use 
#saveRDS(SPOTlight1,"merfish_SPOTlight_100um.RDS")    
saveRDS(SPOTlight1,"merfish_SPOTlight_50um.RDS")


dwls_res<- suppressWarnings(DWLS_run(database, parallel = FALSE, 
                                     is_select_DEGs = FALSE, python_env))

#saveRDS(dwls_res,"merfish_DWLS_100um.RDS")    
saveRDS(dwls_res,"merfish_DWLS_50um.RDS")


### 
#read all data results and 
if(i==1){
  
  Results.Deconv = list()
  
  cell_type = colnames(readRDS("merfish_CARD_100um.RDS"))
  
  Results.Deconv$CARD = readRDS("merfish_CARD_100um.RDS")
  Results.Deconv$RCTD  = readRDS("merfish_RCTD_100um.RDS")
  Results.Deconv$SpatialDWLS  = readRDS("merfish_SpatialDWLS_100um.RDS")
  Results.Deconv$SPOTlight  =readRDS("merfish_SPOTlight_100um.RDS")
  
  ## stereo
  Stere_ct_temp = read.csv("stereoscope_stereo_st100.csv",header = TRUE)
  Stere_ct = Stere_ct_temp[,-1]
  cell_label = Stere_ct_temp[,1] 
  rownames(Stere_ct) =  cell_label
  Results.Deconv$Stereoscope = Stere_ct
  
  
  ## cell2loco
  cell2loc_ct_temp = read.csv("cell2location_prop100.csv",header = TRUE)
  cell2loc_ct = cell2loc_ct_temp[,-1]
  cell_label = cell2loc_ct_temp[,1] 
  rownames(cell2loc_ct) =  cell_label
  Results.Deconv$cell2location = cell2loc_ct
  
  ## get the ensemble results
  fish_ensemble = solve_ensemble (Results.Deconv)
  
  # For R code use 
  saveRDS(fish_ensemble,"merfish_ensemble_100um.RDS")  
  
  Results.Deconv$EnDecon = fish_ensemble$H_norm
  
  ## DWLS
  Results.Deconv$DWLS = readRDS("merfish_DWLS_100um.RDS") 
  
  ## SONAR
  SONAR_ct_temp = read.csv("SONAR.results_100um.csv",header = TRUE)
  SONAR_ct = SONAR_ct_temp[,-1]
  cell_label = SONAR_ct_temp[,1] 
  rownames(SONAR_ct) =  cell_label
  Results.Deconv$SONAR = SONAR_ct
  
  saveRDS(Results.Deconv,"merfish_Results.Deconv_100um.RDS")  
}


if(i==2){
  
  Results.Deconv = list()
  
  Results.Deconv$CARD = readRDS("merfish_CARD_50um.RDS")  
  Results.Deconv$RCTD  = readRDS("merfish_RCTD_50um.RDS")
  Results.Deconv$SpatialDWLS  = readRDS("merfish_SpatialDWLS_50um.RDS")
  Results.Deconv$SPOTlight  =readRDS("merfish_SPOTlight_50um.RDS") 
  
  ## stereo
  Stere_ct_temp = read.csv("stereoscope_stereo_st50.csv",header = TRUE)
  Stere_ct = Stere_ct_temp[,-1]
  cell_label = Stere_ct_temp[,1] 
  rownames(Stere_ct) =  cell_label
  Results.Deconv$Stereoscope = Stere_ct
  
  
  ## cell2loco
  cell2loc_ct_temp = read.csv("cell2location_prop50.csv",header = TRUE)
  cell2loc_ct = cell2loc_ct_temp[,-1]
  cell_label = cell2loc_ct_temp[,1] 
  rownames(cell2loc_ct) =  cell_label
  Results.Deconv$cell2location = cell2loc_ct
  
  ## get the ensemble results
  fish_ensemble = solve_ensemble (Results.Deconv)
  
  # For R code use 
  saveRDS(fish_ensemble,"merfish_ensemble_50um.RDS")  
  
  Results.Deconv$EnDecon = fish_ensemble$H_norm
  
  
  ## DWLS
  Results.Deconv$DWLS = readRDS("merfish_DWLS_50um.RDS") 
  
  ## SONAR
  SONAR_ct_temp = read.csv("SONAR.results_50um.csv",header = TRUE)
  SONAR_ct = SONAR_ct_temp[,-1]
  cell_label = SONAR_ct_temp[,1] 
  rownames(SONAR_ct) =  cell_label
  Results.Deconv$SONAR = SONAR_ct
  
  names(Results.Deconv)
  
  saveRDS(Results.Deconv,"merfish_Results.Deconv_50um.RDS")  
}
