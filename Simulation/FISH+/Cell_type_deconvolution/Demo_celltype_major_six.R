
remove(list=ls())
#data process
data_process <- function(sc_exp, sc_label, spot_exp, spot_loc,
                         gene_det_in_min_cells_per = 0.01, expression_threshold = 1,
                         nUMI = 100, verbose = FALSE,depthscale = depthscale,
                         clean.only = TRUE){
  
  #gene by cell matrix
  if(ncol(sc_exp) != length(sc_label))
    stop("Require cell labels!")
  
  if(ncol(spot_exp) != nrow(spot_loc))
    stop("Require x , y coordinations")
  
  #### scRNA-seq data process
  sc_matrix = cleanCounts(sc_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                          expression_threshold = expression_threshold, nUMI = nUMI,
                          verbose = verbose,depthscale = depthscale,
                          clean.only = clean.only)
  
  
  
  sc_matrix= as.matrix(sc_matrix)
  ind = match(colnames(sc_matrix), colnames(sc_exp))
  sc_label = sc_label[ind]
  
  # cell_type = sort(unique(sc_label))
  
  #### SRT data process
  st_matrix = cleanCounts(spot_exp, gene_det_in_min_cells_per = gene_det_in_min_cells_per,
                          expression_threshold = expression_threshold, nUMI = nUMI,
                          verbose = verbose,depthscale = depthscale,
                          clean.only = clean.only)
  
  st_matrix= as.matrix(st_matrix)
  ind_sp = match(colnames(st_matrix), colnames(spot_exp))
  spot_loc = spot_loc[ind_sp, ]
  
  #### find common genes
  com_gene = intersect(rownames(sc_matrix),rownames(st_matrix))
  sc_exp = sc_matrix[com_gene,]
  st_exp = st_matrix[com_gene,]
  
  ### rechecking nUMI
  index_sc <- colSums(sc_exp) >= nUMI
  sc_exp_filter <- sc_exp[,index_sc]
  sc_label_filter <- sc_label[index_sc]
  
  index_st <- colSums(st_exp) >= nUMI
  st_exp_filter = st_exp[,index_st]
  spot_loc_filter <- spot_loc[index_st,]
  
  database <- list(sc_exp = sc_exp_filter, sc_label = sc_label_filter,
                   spot_exp = st_exp_filter, spot_loc = spot_loc_filter)
  return(database)
}


cleanCounts <- function (counts, gene_det_in_min_cells_per = 0.01,
                         expression_threshold = 1 ,
                         nUMI = 100, verbose = FALSE, depthscale = 1,
                         clean.only = FALSE) {
  
  #counts: gene by cells matrix
  
  if (clean.only){
    
    n = ncol(counts)
    
    ##### select of the genes
    filter_index_genes = Matrix::rowSums(counts > expression_threshold) >
      gene_det_in_min_cells_per*n
    
    #### filter the cell
    filter_index_cells = Matrix::colSums(counts[filter_index_genes,] >
                                           expression_threshold) > nUMI
    
    x = counts[ filter_index_genes,filter_index_cells]
    
    if (verbose) {
      message("Resulting matrix has ", ncol(counts), " cells and ", nrow(counts), " genes")
    }
    
    return(x) 
    
  } else {
    
    n = ncol(counts)
    
    ##### select of the genes
    filter_index_genes = Matrix::rowSums(counts > expression_threshold) >
      gene_det_in_min_cells_per*n
    
    #### filter the cell
    filter_index_cells = Matrix::colSums(counts[filter_index_genes,] >
                                           expression_threshold) >nUMI
    
    x = counts[ filter_index_genes,filter_index_cells]
    
    
    sf <- colSums(x)/median(colSums(x))
    
    return(log(sweep(x, 2, sf, '/')*depthscale+1))
    
  }
  
}




source("runindi.R")

library(CARD)
library("spacexr")
library(Seurat)
library(SPOTlight)
library(Giotto)

python_env <- 'C:/Users/visitor01/.conda/envs/ccnet/python.exe'
reticulate::use_python(python = python_env, required = TRUE)

#load data

database_temp = readRDS("SS_fish_plus_spot.Rds")

database = data_process(sc_exp = database_temp$ref_data$sc_exp, 
                        sc_label = as.factor(database_temp$ref_data$sc_lable), 
                        spot_exp = database_temp$spot_data$spot_exp,
                        spot_loc = as.data.frame(database_temp$spot_data$spot_loc),
                        gene_det_in_min_cells_per = 0.01, expression_threshold = 1,
                        nUMI = 100, verbose = FALSE, depthscale = 1, 
                        clean.only = TRUE)



## use for  steroscope and cell2location
spot_exp_clear <- database$spot_exp
spot_loc_clear <- database$spot_loc
ref_exp_clear <- database$sc_exp
ref_label_clear  <- database$sc_label 

write.csv(spot_exp_clear,"spot_exp_clear.csv")
write.csv(spot_loc_clear,"spot_loc_clear.csv")
write.csv(ref_exp_clear,"ref_exp_clear.csv")
write.csv(ref_label_clear,"ref_label_clear.csv")


##run individuall models

SpatialDWLS1 <- suppressWarnings(spatialDWLS_run(database, 
                                         my_python_path = python_env))

saveRDS(SpatialDWLS1,"fish_SpatialDWLS.RDS")   


RCTD1 <- RCTD_run(database)

saveRDS(RCTD1,"fish_RCTD.RDS")   


CARD1 <- CARD_run(database)

saveRDS(CARD1,"fish_CARD.RDS")



SPOTlight1 <- suppressWarnings(SPOTlight_run(database))

# For R code use 
saveRDS(SPOTlight1,"fish_SPOTlight.RDS")    


### 
#read all data resultes and 
source("runindi.R")
Results.Deconv = list()

Results.Deconv$CARD = readRDS("fish_SPOTlight.RDS")  
Results.Deconv$RCTD  = readRDS("fish_CARD.RDS")
Results.Deconv$SpatialDWLS  = readRDS("fish_SpatialDWLS.RDS")   
Results.Deconv$SPOTlight  = readRDS("fish_SPOTlight.RDS")  

 
## stereo
Stere_ct_temp = read.csv("stereoscope_stereo_st.csv",header = TRUE)
Stere_ct = Stere_ct_temp[,-1]
cell_label = Stere_ct_temp[,1] 
rownames(Stere_ct) =  cell_label
Results.Deconv$Stereoscope = Stere_ct


## cell2loco
cell2loc_ct_temp = read.csv("cell2location_prop.csv",header = TRUE)
cell2loc_ct = cell2loc_ct_temp[,-1]
cell_label = cell2loc_ct_temp[,1] 
rownames(cell2loc_ct) =  cell_label
Results.Deconv$cell2location = cell2loc_ct

## get the ensemble results
fish_ensemble = solve_ensemble (Results.Deconv)

# For R code use 
saveRDS(fish_ensemble,"fish_ensemble.RDS")    

