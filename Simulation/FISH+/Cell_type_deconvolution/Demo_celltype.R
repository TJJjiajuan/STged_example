


source("runindi.R")

library(CARD)
library("spacexr")
library(Seurat)
library(SPOTlight)
library(Giotto)

## use the data for the deconvolution

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

Results.Deconv$CARD = readRDS("fish_CARD.RDS")  
Results.Deconv$RCTD  = readRDS("fish_RCTD.RDS")
Results.Deconv$SpatialDWLS  = readRDS("fish_SpatialDWLS.RDS")   
Results.Deconv$SPOTlight  = readRDS("SPOTlight_res.RDS")  

 
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




########## For all results
#read all data results and 
i =1 
if(i==1){
  
  Results.Deconv = list()
  
  celltype = sort(colnames(readRDS("fish_CARD.RDS")))
  
  Results.Deconv$CARD = readRDS("fish_CARD.RDS")[,celltype]  
  Results.Deconv$RCTD  = readRDS("fish_RCTD.RDS")[,celltype]  
  Results.Deconv$SpatialDWLS  = readRDS("fish_SpatialDWLS.RDS") [,celltype]  
  Results.Deconv$SPOTlight  = readRDS("SPOTlight_res.RDS")#[,celltype]  
  
  ## stereo
  Stere_ct_temp = read.csv("stereoscope_stereo_st.csv",header = TRUE)
  Stere_ct = Stere_ct_temp[,-1]
  cell_label = Stere_ct_temp[,1] 
  rownames(Stere_ct) =  cell_label
  Results.Deconv$Stereoscope = Stere_ct[,celltype]  
  
  
  ## cell2loco
  cell2loc_ct_temp = read.csv("cell2location_prop.csv",header = TRUE)
  cell2loc_ct = cell2loc_ct_temp[,-1]
  cell_label = cell2loc_ct_temp[,1] 
  rownames(cell2loc_ct) =  cell_label
  Results.Deconv$cell2location = cell2loc_ct[,celltype] 
  
  ## get the ensemble results
  fish_ensemble = solve_ensemble (Results.Deconv)
  
  Results.Deconv$EnDecon = fish_ensemble$H_norm
  
  ## DWLS
  Results.Deconv$DWLS = readRDS("FISH_DWLS.RDS") 
  
  ## SONAR
  SONAR_ct_temp = read.csv("SONAR.results_FISH.csv",header = TRUE)
  SONAR_ct = SONAR_ct_temp[,-1]
  cell_label = SONAR_ct_temp[,1] 
  rownames(SONAR_ct) =  cell_label
  Results.Deconv$SONAR = SONAR_ct
  
  saveRDS(Results.Deconv,"FISH_Results.Deconv.RDS")  
}



