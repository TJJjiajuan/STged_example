
############## apply RCTD #####################
RCTD_run = function(database, CELL_MIN_INSTANCE = 10,min_UMI =10){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <-  droplevels(database$sc_label)
  spot_exp <- database$spot_exp
  cell_type <-  droplevels(sort(unique(sc_label)))
  spot_loc <- database$spot_loc
  
  sparse_sc_exp <- as(sc_exp, "sparseMatrix")
  sparse_spot_exp <- as(spot_exp, "sparseMatrix")
  
  ## the reference scRNA-seq data
  cellnames <- colnames(sc_exp)
  cell_types <- as.factor(sc_label)
  names(cell_types) <- cellnames
  sc_nUMI <- as.numeric(colSums(sc_exp))
  names(sc_nUMI) <- cellnames
  reference <- spacexr::Reference(sparse_sc_exp, 
                                  cell_types, 
                                  nUMI = sc_nUMI, min_UMI = 10)
  
  ### Create SpatialRNA object
  coords <- as.data.frame(spot_loc)
  # coords <- as.data.frame(matrix(1,dim(spot_exp)[2],2))
  rownames(coords) <- as.character(colnames(spot_exp))
  nUMI <- colSums(spot_exp)
  puck <- spacexr::SpatialRNA(coords, counts=sparse_spot_exp, nUMI=nUMI)
  
  myRCTD <- suppressMessages(spacexr::create.RCTD(puck, reference, max_cores = 1,CELL_MIN_INSTANCE = CELL_MIN_INSTANCE))
  myRCTD <- suppressMessages(spacexr::run.RCTD(myRCTD, doublet_mode = 'full'))
  results <- myRCTD@results
  
  temp <- as.matrix(results$weights)
  norm_weights_temp <- sweep(temp, 1, rowSums(temp), '/')
  RCTD_results <- norm_weights_temp[,cell_type]
  
  return(RCTD_results)
}

############## apply SPOTlight #####################
SPOTlight_run = function(database, cl_n = 100,
                         hvg = 100,
                         min_cont = 0.001){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <-  droplevels(sort(unique(sc_label)))
  
  sc_ref <- Seurat::CreateSeuratObject(counts = sc_exp)
  sc_ref@meta.data$subclass <- as.factor(sc_label)
  
  sc_ref <- Seurat::SCTransform(sc_ref, verbose = FALSE)
  
  Seurat::Idents(object = sc_ref) <- sc_ref@meta.data$subclass
  cluster_markers_all <- Seurat::FindAllMarkers(object = sc_ref,
                                                assay = "SCT",
                                                slot = "data",
                                                verbose = FALSE,
                                                only.pos = TRUE)
  
  # Downsample scRNAseq to select gene set and number of cells to train the model
  # with default parameters
  se_sc_down <- SPOTlight::downsample_se_obj(se_obj = sc_ref,
                                             clust_vr = "subclass",
                                             cluster_markers = cluster_markers_all,
                                             cl_n = cl_n,
                                             hvg = hvg)
  
  nmf_mod_ls <- SPOTlight::train_nmf(cluster_markers = cluster_markers_all,
                                     se_sc = se_sc_down,
                                     mtrx_spatial = spot_exp,
                                     clust_vr = "subclass",
                                     ntop = NULL,
                                     hvg = hvg,
                                     transf = "uv",
                                     method = "nsNMF")
  
  nmf_mod <- nmf_mod_ls[[1]]
  
  # get basis matrix W
  w <- NMF::basis(nmf_mod)
  
  # get coefficient matrix H
  h <- NMF::coef(nmf_mod)
  
  # Subset to genes used to train the model
  temp_index <- ! is.na(pmatch(rownames(spot_exp), rownames(w)))
  spot_exp_train <- spot_exp[temp_index, ]
  
  ct_topic_profiles <- SPOTlight::topic_profile_per_cluster_nmf(h = h,
                                                                train_cell_clust = nmf_mod_ls[[2]])
  
  decon_mtrx <- SPOTlight::mixture_deconvolution_nmf(nmf_mod = nmf_mod,
                                                     mixture_transcriptome = spot_exp_train,
                                                     transf = "uv",
                                                     reference_profiles = ct_topic_profiles,
                                                     min_cont = min_cont)
  
  SPOTlight_results <- as.matrix(decon_mtrx[,- ncol(decon_mtrx)])
  rownames(SPOTlight_results) <- colnames(spot_exp)
  colnames(SPOTlight_results) <- cell_type
  
  return(SPOTlight_results)
}



############## apply SpatialDWLS #####################
spatialDWLS_run = function(database, my_python_path, is_select_DEGs = FALSE,
                           findmarker_method = "gini",
                           ncp_spa = 100, dimensions_to_use = 10, k = 10,
                           resolution = 0.4, n_iterations = 1000, n_cell = 50){
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  cell_type <-  droplevels(sort(unique(sc_label)))
  
  instrs = Giotto::createGiottoInstructions(python_path = my_python_path)
  
  sc_cortex <- Giotto::createGiottoObject(raw_exprs = sc_exp,
                                          instructions = instrs)
  
  sc_cortex <- Giotto::normalizeGiotto(gobject = sc_cortex)
  
  sc_cortex@cell_metadata$leiden_clus <- as.character(sc_label)
  
  if(is_select_DEGs){
    gini_markers_subclusters = Giotto::findMarkers_one_vs_all(gobject = sc_cortex,
                                                              method = findmarker_method,
                                                              expression_values = 'normalized',
                                                              cluster_column = 'leiden_clus',
                                                              verbose = FALSE)
    topgenes_gini = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']
    sc_norm_exp <- 2^(sc_cortex@norm_expr)-1
    ExprSubset <- sc_norm_exp[as.character(topgenes_gini$genes),]
  }else{
    sc_norm_exp <- 2^(sc_cortex@norm_expr)-1
    ExprSubset <- sc_norm_exp
  }
  
  
  Sig<-NULL
  for (i in as.character(unique(sc_label))){
    Sig<-cbind(Sig,(apply(ExprSubset,1,function(y) mean(y[which(sc_label==i)]))))
  }
  colnames(Sig)<-as.character(unique(sc_label))
  
  
  
  grid_seqFish <- Giotto::createGiottoObject(raw_exprs = spot_exp,instructions = instrs)
  grid_seqFish <- Giotto::normalizeGiotto(gobject = grid_seqFish)
  grid_seqFish <- Giotto::calculateHVG(gobject = grid_seqFish,show_plot = FALSE,
                                       return_plot =  FALSE)
  gene_metadata = Giotto::fDataDT(grid_seqFish)
  featgenes = gene_metadata[hvg == 'yes']$gene_ID
  grid_seqFish <- Giotto::runPCA(gobject = grid_seqFish, genes_to_use = featgenes, scale_unit = F, ncp = ncp_spa)
  grid_seqFish <- Giotto::createNearestNetwork(gobject = grid_seqFish, dimensions_to_use = 1:dimensions_to_use, k = k)
  grid_seqFish <- Giotto::doLeidenCluster(gobject = grid_seqFish, resolution = resolution, n_iterations = n_iterations)
  
  grid_seqFish <- Giotto::runDWLSDeconv(gobject = grid_seqFish, sign_matrix = Sig, n_cell = n_cell)
  
  spatialDWLS_results <- as.matrix(grid_seqFish@spatial_enrichment$DWLS[,-1])
  spatialDWLS_results <- spatialDWLS_results[,cell_type]
  rownames(spatialDWLS_results) <- colnames(spot_exp)
  
  return(spatialDWLS_results)
}

################ run CARD ################
CARD_run <- function(database, minCountGene = 10, minCountSpot = 5){
  ### extract SC, CL, ST, ST.Loc from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  spot_loc <- database$spot_loc
  cell_type <-  droplevels(sort(unique(sc_label)))
  spot_names_ori <- colnames(spot_exp)
  ###
  spot_names <- paste0("spot", 1:nrow(spot_loc))
  colnames(spot_exp) <- spot_names
  rownames(spot_loc) <- spot_names
  colnames(spot_loc) <- c("x", "y")
  
  sampleInfo <- rep("sample1", ncol(sc_exp))
  names(sampleInfo) <- colnames(sc_exp)
  meta.data <- data.frame(cellID = colnames(sc_exp), cellType = sc_label, sampleInfo = sampleInfo)
  rownames(meta.data) <- colnames(sc_exp)
  
  CARD_obj <- CARD::createCARDObject(sc_count = sc_exp,
                                     sc_meta = meta.data,
                                     spatial_count = spot_exp,
                                     spatial_location = spot_loc,
                                     ct.varname = "cellType",
                                     ct.select = cell_type,
                                     sample.varname = "sampleInfo",
                                     minCountGene = minCountGene,
                                     minCountSpot = minCountSpot)
  CARD_obj <- CARD::CARD_deconvolution(CARD_object = CARD_obj)
  
  CARD_results <- as.matrix(CARD_obj@Proportion_CARD)
  # ccc <- setdiff(colnames(spot_exp),colnames(CARD_obj@spatial_countMat))
  #
  # cc <- spot_exp[,ccc]
  CARD_results <- CARD_results[colnames(spot_exp), cell_type]
  rownames(CARD_results) <- spot_names_ori
  return(CARD_results)
}

SONAT_run= function(database, minCountSpot = 5){
 
   ### extract SC, CL, ST, ST.Loc from database
  sc_exp <- database$sc_exp
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp
  spot_loc <- database$spot_loc
  cell_type <-  droplevels(sort(unique(sc_label)))
  spot_names_ori <- colnames(spot_exp)
  ###
  spot_names <- paste0("spot", 1:nrow(spot_loc))
  colnames(spot_exp) <- spot_names
  rownames(spot_loc) <- spot_names
  colnames(spot_loc) <- c("x", "y")
  
  #calculate the nUMI and nUMI_spot
  nUMI <- colSums(sc_exp)
  names(nUMI) <- colnames(sc_exp)
  nUMI_spot <- colSums(spot_exp)
  names(nUMI_spot) <- colnames(spot_exp)
  
  #preprocess the input data
  processed_input<-SONAR.preprocess(sc_count= sc_exp,
                                    sc_cell_type=sc_label,
                                    sc_nUMI=nUMI,
                                    sp_coords = spot_loc,
                                    sp_count = spot_exp,
                                    sp_nUMI = nUMI_spot,
                                    cores=1)
  
  #deliver the preprocessed data to SONAR
  SONAR.deliver(processed_data=processed_input,path=code_path)
  
  #define the bandwidth, default is 1.2 times minimal distance
  temp<-dist(coords)
  temp<-Matrix::Matrix(temp)
  temp[temp==0] <- NA
  mindist <- min(temp,na.rm = T)
  h <- 1.2*mindist
}
####### Utils of data process of DWLS ##########
create_group_exp <- function(sc_exp,sc_label) {
  
  #sc_exp single cell gene expression datasets
  #sc_label  cell annotation of the single cells of the reference
  
  ##group cells
  # reference matrix (C) + refProfiles.var from TRAINING dataset
  cell_type = sort(unique(sc_label))
  group = list()
  for(i in 1:length(cell_type)){
    temp_use <- which(sc_label == cell_type[i])
    names(temp_use) <- NULL
    group[[i]] <- temp_use
  }
  sc_group_exp = sapply(group,function(x) Matrix::rowMeans(sc_exp[,x]))
  #sapply
  sc_group_exp = as.matrix(sc_group_exp)
  colnames(sc_group_exp) = cell_type
  return(sc_group_exp)
}

########## run DWLS ##############
DWLS_run = function(database, parallel = TRUE,  is_select_DEGs = TRUE, python_env){
  
  ### extract SC, CL, ST from database
  sc_exp <- database$sc_exp * 1e-4
  sc_label <- database$sc_label
  spot_exp <- database$spot_exp * 1e-4
  cell_type <- sort(unique(sc_label))
  
  
  if(is_select_DEGs){
    instrs = Giotto::createGiottoInstructions(python_path = python_env)
    
    sc_cortex <- Giotto::createGiottoObject(raw_exprs = sc_exp,
                                            instructions = instrs)
    
    sc_cortex <- Giotto::normalizeGiotto(gobject = sc_cortex)
    
    sc_cortex@cell_metadata$leiden_clus <- as.character(sc_label)
    
    gini_markers_subclusters = Giotto::findMarkers_one_vs_all(gobject = sc_cortex,
                                                              method = 'gini',
                                                              expression_values = 'normalized',
                                                              cluster_column = 'leiden_clus',
                                                              verbose = FALSE)
    topgenes_mast = gini_markers_subclusters[, head(.SD, 100), by = 'cluster']
    
    sc_exp <- sc_exp[topgenes_mast$genes, ]
    spot_exp <- spot_exp[topgenes_mast$genes, ]
  }
  cellAnnots <- data.frame(CellID = colnames(sc_exp),
                           cellType = sc_label)
  
  cell_type_exp <- create_group_exp(sc_exp,sc_label)
   
  if(parallel){
    #### because the cores can't large than 2 in using R CMD check processing
    chk <- Sys.getenv("_R_CHECK_LIMIT_CORES_", "")
    if (nzchar(chk) && chk == "TRUE") {
      # use 2 cores in CRAN/Travis/AppVeyor
      num.cores <- 2L
    } else {
      # use all cores in devtools::test()
      num.cores <- parallel::detectCores() - 1
    }
    myCluster <- parallel::makeCluster(num.cores)
    doParallel::registerDoParallel(myCluster)
    
    pred.DWLS <- Sys.time()
    DampenedWLS <- foreach::foreach(i=1:ncol(spot_exp), .combine = 'rbind', .inorder = TRUE) %dopar% {
      DWLS::solveDampenedWLS(cell_type_exp, spot_exp[,i])
    }
    end.DWLS <- Sys.time()
    time.DWLS <- difftime(end.DWLS, pred.DWLS, units = "mins")
    
    parallel::stopCluster(myCluster)
    
  }else{
    
   DampenedWLS = matrix(NA, ncol(spot_exp), length(cell_type))
    pred.DWLS <- Sys.time()
    for (i in seq(ncol(spot_exp))) {
      DampenedWLS[i,] <- DWLS::solveDampenedWLS(cell_type_exp, spot_exp[,i])
    }
    end.DWLS <- Sys.time()
    time.DWLS <- difftime(end.DWLS, pred.DWLS, units = "mins")
    
  }
  
  rownames(DampenedWLS) <- colnames(spot_exp)
  
  colnames(DampenedWLS) <- cell_type
  
  return(DampenedWLS)
}




solve_ensemble <- function(Results.Deconv, lambda = NULL, prob.quantile = 0.5,
                           niter = 100, epsilon = 1e-5){
  # Results.Deconv <- Results.Deconv.all[[1]]
  num.methods <- length(Results.Deconv)
  num.spots <- nrow(Results.Deconv[[1]])
  num.cell.type <- ncol(Results.Deconv[[1]])
  
  ## initialization V by the mean of individual values
  w <- c(rep(1/num.methods, num.methods))
  H <-  Reduce("+", Map("*", Results.Deconv, w))
  
  if(is.null(lambda)){
    cat("We will adpote a value for lambda in our algorithm...", "\n")
  }
  
  k <- 1
  while (k <= niter) {
    if(k == 1){
      loss_all_temp <- 0
      temp2 <-  sapply(Results.Deconv, L1_norm, Y = H)
      lambda <- quantile(temp2, probs = prob.quantile)
    }else{
      loss_all_temp <- loss_all
    }
    ##### update w
    temp2 <-  sapply(Results.Deconv, L1_norm, Y = H)
    w <- exp(-temp2/lambda)/sum(exp(-temp2/lambda))
    ##### update V
    temp4 <- do.call(abind::abind, c(Results.Deconv, list(along = 3)))
    H <- apply(temp4, 1:2, median_weighted, w = w)
    
    loss_main <- sum(sapply(Results.Deconv, L1_norm, Y = H) * w)
    loss_entropy <- sum(w * log(w))
    loss_all <- loss_main + lambda * loss_entropy
    
    cat("iter: ", k, "loss_main: ", loss_main, "loss_entropy: ", loss_entropy,
        "loss_all: ", loss_all, "lambda: ", lambda, "\n")
    if(k == niter)
      cat("The method maybe not convergens, the algorithm need an larger max_epoches!", "\n")
    
    if(abs(loss_all - loss_all_temp) < epsilon | k >= niter)
      break
    k <- k + 1
  }
  
  colnames(H) <- colnames(Results.Deconv[[1]])
  H_norm <- sweep(H, 1, rowSums(H), "/")
  return(list(H_norm = H_norm, w = w))
}

L1_norm <- function(X, Y){
  return(sum(abs(X-Y)))
}

median_weighted <- function(x, w){
  index_non_zero <- which(x != 0, arr.ind = TRUE)
  x_use <- x[index_non_zero]
  w_use <- w[index_non_zero]
  results <- weighted.median(x,w)#spatstat.geom::weighted.median(x,w)
  return(results)
}