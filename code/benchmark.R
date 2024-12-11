res_common_gene = function(estres, commgene){
  
  K = length(estres)
  
  for(i in 1:K){
    
    temp  =  estres[[i]][commgene,]
    
    estres[[i]] = temp
  }
  
  
  return(estres)
}

res_library_size_normaization = function(estres, depthscale = 1e6, log = TRUE ){
  
  K = length(estres)
  if(log){
    for(i in 1:K){
      
      sf <- colSums(estres[[i]])  
      
      temp  = log(sweep(estres[[i]], 2, sf, '/')*depthscale+1)
      
      temp[is.na(temp)] <- 0
      
      estres[[i]] = temp
    }
  }else{
    for(i in 1:K){
      
      sf <- colSums(estres[[i]])  
      
      temp  = sweep(estres[[i]], 2, sf, '/')
      
      temp[is.na(temp)] <- 0
      
      estres[[i]] = temp*depthscale
    } 
  }
  
  return(estres)
}


RCTD_run = function(sc_exp  = sc_exp, 
                    sc_label=  sc_label, spot_exp=  spot_exp,
                    spot_loc  = spot_loc ,CELL_MIN_INSTANCE = 20){
  ### extract SC, CL, ST from database
  
  cell_type <- sort(unique(sc_label))
  
  sparse_sc_exp <- as(sc_exp, "sparseMatrix")
  sparse_spot_exp <- as(spot_exp, "sparseMatrix")
  
  ## the reference scRNA-seq data
  cellnames <- colnames(sc_exp)
  cell_types <- as.factor(sc_label)
  names(cell_types) <- cellnames
  sc_nUMI <- as.numeric(colSums(sc_exp))
  names(sc_nUMI) <-cellnames
  reference <- spacexr::Reference(sparse_sc_exp, cell_types, nUMI = sc_nUMI)
  
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

###### the result from RCTD
RCTDexp = function(beta,  srt_exp, ref_mu, cutoff = 0.05 ){
  
  # beta is a spot_by_cell type matrix
  # srt_exp is a gene_by_spot matrix
  # ref_mu is a gene_by_cell type matrix
  # the out put   exp.est is a list with length K, and each list with a gene_by_spot gene expression matrix
  beta[beta<cutoff] = 0
  epsilon = 1e-10
  # for the whole cell types 
  denom = beta %*% t(ref_mu)+epsilon
  
  K = ncol(ref_mu)
  
  exp.est = matrix(list(), K,1)
  names( exp.est) = colnames(ref_mu)
  
  for(k in 1: K) {
    
    posterior = (beta[,k] %*% t(ref_mu[,k])) / denom
    
    exp.est [[k]] = t( posterior*t(srt_exp))
    
  }
  
  return(  exp.est )
  
}

st_mu_est = function(srt_exp, beta, cutoff = 0.05){
  
  beta[beta<cutoff] = 0
  if(is.null(colnames(beta))){
    
    colnames(beta) = paste0("celltype",1:ncol(beta))
  }
  beta_ind = beta > 0
  
  x = beta
  y = t(srt_exp)
  
  fit1 =lsfit(x, y,int=F)
  
  st_mu = fit1$coef
  
  #reshape to the list with the all cell type
  
  K = ncol(beta)
  n = nrow(beta)
  p = ncol(y )
  F_list = matrix(list(), K ,1)
  names(F_list) = colnames(beta)
  genenames = rownames(srt_exp)
  spot =  rownames(beta)
  
  for (i in 1:K){
    
    temp = matrix( rep(st_mu[i,],n), nrow = p, byrow = FALSE)
    temp1 =  matrix( rep(beta_ind[,i],p), nrow = p, byrow = TRUE)
    
    temp_res = temp* temp1
    rownames(temp_res) =  genenames 
    colnames(temp_res) =   spot  
    F_list[[i]]  = temp_res
    
  }
  
  
  
  return(F_list) 
}


sc_mu_est = function( sc_mu, beta, cutoff = 0.05){
  
  
  beta[beta<cutoff] = 0
  if(is.null(colnames(beta))){
    
    colnames(beta) = paste0("celltype",1:ncol(beta))
  }
  
  beta_ind = beta > 0
  #reshape to the list with the all cell type
  
  K = ncol(beta)
  n = nrow(beta)
  p = nrow(sc_mu)
  F_list = matrix(list(), K ,1)
  names(F_list) = colnames(beta)
  genenames = rownames(sc_mu)
  spot =  rownames(beta)
  for (i in 1:K){
    temp =   matrix( rep(sc_mu[,i],n), ncol =  n, byrow = FALSE)
    
    temp1 =  matrix( rep(beta_ind[,i],p), nrow = p, byrow = TRUE)
    temp_res = temp* temp1
    rownames(temp_res) =  genenames 
    colnames(temp_res) =   spot  
    F_list[[i]]  = temp_res
  }
  
  return(F_list) 
}


bMIND_est = function( spot_exp, ref_exp, beta,ncore = 1){
  
  #beta[beta<cutoff] = 0
  cell_type = colnames(beta)
  exp_bmind = bMIND(bulk = spot_exp, profile =  ref_exp, frac = beta,ncore=ncore)$A
  
  dims <- dim(exp_bmind) 
  exp_mind_list <- vector("list", length = dims[2])
  names(exp_mind_list) <- cell_type
  
  for (i in 1:dims[2]) {
    exp_matrix <- exp_bmind[,  i,]
    exp_mind_list[[i]] <- exp_matrix
  }
  
  
  return(exp_mind_list) 
}


ENIGMA_est = function( spot_exp, ref_exp, beta){
  
  cell_type = colnames(beta)
  ENIGMA_trace.v <- cell_deconvolve_trace(O = spot_exp,
                                          theta=beta,
                                          R=ref_exp)
  
  dims <- dim(ENIGMA_trace.v) 
  exp_ENIGMA_list <- vector("list", length = dims[3])
  names(exp_ENIGMA_list) <- cell_type
  
  for (i in 1:dims[3]) {
    exp_matrix <- ENIGMA_trace.v[, , i]
    exp_ENIGMA_list[[i]] <- exp_matrix
  }
  
  
  return(exp_ENIGMA_list) 
}


spotdecon_est = function(srt_exp, beta, cutoff = 0.05){
  
  beta[beta<cutoff] = 0
  
  K = ncol(beta)
  
  p = nrow(srt_exp )
  F_list = matrix(list(), K ,1)
  names(F_list) = colnames(beta)
  genenames = rownames(srt_exp)
  spot =  rownames(beta)
  
  for (i in 1:K){
    
    temp =  matrix( rep(beta[,i],p), nrow = p, byrow = TRUE)
    temp_res = srt_exp * temp
    rownames(temp_res) =  genenames 
    colnames(temp_res) =   spot  
    F_list[[i]]  = temp_res
    
    
  }
  
  return(F_list) 
}

corr.cal.spot <- function(F.est, F.true, beta.ind, row_filter = NULL) {
  
  # Validate input dimensions
  if (length(F.est) != length(F.true) || length(F.est) != length(beta.ind)) {
    stop("All input lists must have the same length.")
  }
  
  # Apply row_filter if provided
  if (!is.null(row_filter)) {
    F.est.filtered <- Map(function(mat, filter) mat[filter, , drop = FALSE], F.est, asplit(row_filter, 1))
    F.true.filtered <- Map(function(mat, filter) mat[filter, , drop = FALSE], F.true, asplit(row_filter, 1))
    beta.ind.filtered <- Map(function(mat, filter) mat[filter, , drop = FALSE], beta.ind, asplit(row_filter, 1))
  } else {
    F.est.filtered <- F.est
    F.true.filtered <- F.true
    beta.ind.filtered <- beta.ind
  }
  
  # Combine all filtered matrices
  F.est.all <- do.call(rbind, F.est.filtered)
  F.true.all <- do.call(rbind, F.true.filtered)
  beta.ind.all <- do.call(rbind, beta.ind.filtered)
  
  # Safe log function to handle small values and avoid issues
  safe_log <- function(x, base = 2) {
    x[x < 1e-10] <- 1e-10  # Replace very small values with a threshold to avoid log(0)
    return(log(x, base = base))
  }
  
  # Calculate metrics for each column
  calc_metrics <- function(i) {
    
    P <- F.est.all[, i][beta.ind.all[, i]]
    Q <- F.true.all[, i][beta.ind.all[, i]]
    
    # Handle edge cases
    if (length(P) < 2 || length(Q) < 2 || any(is.na(P)) || any(is.na(Q)) || var(P) == 0 || var(Q) == 0) {
      return(rep(NA, 4))
    }
    
    rmse <- sqrt(mean((P - Q) ^ 2, na.rm = TRUE))
    
    # Normalize P and Q ONLY for Jensen-Shannon Divergence calculation
    epsilon <- 1e-10  # Small value to avoid division by zero
    P_norm <- (P + epsilon) / sum(P + epsilon)
    Q_norm <- (Q + epsilon) / sum(Q + epsilon)
    
    # Calculate metrics without normalization (RMSE, Pearson, Cosine)
    pearson_corr <- cor(P_norm, Q_norm, use = "complete.obs")
    cosine_sim <- 1 - abdiv::cosine_distance(P_norm, Q_norm)
    
    # Calculate Jensen-Shannon divergence
    jsd <- tryCatch({
      M <- 0.5 * (P_norm + Q_norm)
      jsd_value <- 0.5 * (sum(P_norm * safe_log(P_norm / M)) + sum(Q_norm * safe_log(Q_norm / M)))
      if (jsd_value < 0) {
        jsd_value <- 0
      }
      jsd_value
    }, error = function(e) NaN)
    
    return(c(pearson_corr, cosine_sim, rmse, jsd))
  }
  
  # Apply the metrics calculation across columns
  est.res <- t(vapply(1:ncol(F.est[[1]]), calc_metrics, numeric(4)))
  
  return(est.res)
}


save_spot_level_results <- function(cell_type_matrices, normalize = TRUE) {
  
  # Number of cell types (length of the list)
  k <- length(cell_type_matrices)
  
  # Number of genes (rows) and spots (columns) in each matrix
  p <- nrow(cell_type_matrices[[1]])
  n <- ncol(cell_type_matrices[[1]])
  
  # Create a list to store spot-level matrices
  spot_level_matrices <- vector("list", n)
  
  # Iterate over each spot
  for (spot_idx in 1:n) {
    # For each spot, extract the column from each cell type matrix
    spot_matrix <- sapply(cell_type_matrices, function(mat) mat[, spot_idx])
    
    # Optionally normalize the spot matrix
    if (normalize) {
      # Normalize each row (gene expression) by the sum of that row across all cell types
      spot_matrix <- spot_matrix/sum(spot_matrix)
        #sweep(spot_matrix, 1, rowSums(spot_matrix), FUN = "/")
      #spot_matrix[is.na(spot_matrix)] <- 0  # Handle any divisions by zero
    }
    
    # Store the result as a p by k matrix (genes x cell types)
    spot_level_matrices[[spot_idx]] <- spot_matrix
  }
  
  return(spot_level_matrices)
}

evaluate_spot_metrics <- function(spot_F_true, spot_F_est, spot_beta_ind) {
  
  # Function to compute various metrics for a single spot
  compute_metrics <- function(P, Q) {
    
    # Handle edge cases
    if (length(P) < 2 || length(Q) < 2 || any(is.na(P)) || any(is.na(Q)) || var(P) == 0 || var(Q) == 0) {
      return(c(pearson_corr = NA, cosine_sim = NA, rmse = NA, jsd = NA))
    }
    
    # Calculate Pearson correlation
    pearson_corr <- cor(P, Q, use = "complete.obs", method =  "spearman")

    
    # Calculate cosine similarity
    cosine_sim <- 1 - abdiv::cosine_distance(P, Q)
    
    # Calculate RMSE
    rmse <- sqrt(mean((P - Q) ^ 2, na.rm = TRUE))
    
    # Normalize for Jensen-Shannon Divergence (JSD)
    epsilon <- 1e-12
    P_norm <- (P + epsilon) / sum(P + epsilon)
    Q_norm <- (Q + epsilon) / sum(Q + epsilon)
    
    # Calculate Jensen-Shannon Divergence (JSD)
    jsd <- tryCatch({
      M <- 0.5 * (P_norm + Q_norm)
      jsd_value <- 0.5 * (sum(P_norm * log(P_norm / M)) + sum(Q_norm * log(Q_norm / M)))
      if (jsd_value < 0) {
        jsd_value <- 0
      }
      jsd_value
    }, error = function(e) NaN)
    
    return(c(pearson_corr = pearson_corr, cosine_sim = cosine_sim, rmse = rmse, jsd = jsd))
  }
  
  # Use lapply to apply compute_metrics over all spots
  results_list <- sapply(1:length(spot_F_true), function(i) {
    beta_ind <- spot_beta_ind[[i]]
    P <- spot_F_true[[i]][beta_ind]
    Q <- spot_F_est[[i]][beta_ind]
    
    # Compute metrics for this spot
    metrics <- compute_metrics(P, Q)
    
    # Optionally print or store results for each spot
   # print(paste("Spot", i, "Metrics:"))
   # print(metrics)
    
    # Return the metrics for this spot
    return(metrics)
  })
  
  # Return the list of results 
   results_list = t(results_list)
  return(results_list)
}





#################
cellcorr.plot = function(corr.list, beta.type, title = "Cell", ytitle ="Spearman correlation",   plot= TRUE){
  
  beta.ind.temp = beta.type>0.05
  nnum = apply(beta.ind.temp,2,sum)
  
  K = length(corr.list)
  methods = names(corr.list)
  cell_type = colnames(beta.type)
  lebel_method = rep(cell_type, c(nnum))
  
  cellcorrs <- data.frame(corr = unlist(corr.list),
                          lable = factor(rep(lebel_method,  K),levels = cell_type),
                          methods = factor(rep(methods, each = sum(nnum)),levels =  methods ))
  
  p1 <- ggplot(cellcorrs,  aes(x =lable, y = corr))+ 
    geom_boxplot()+ aes(fill = methods)+
    labs(title = title , x=" ", y= ytitle) + theme_classic()+
    theme(
      axis.text.x.bottom = element_text(size = 12,hjust = 1,angle =45), 
      axis.text.y.left = element_text(size = 12),
      axis.title.x = element_blank(),##element_text(size = 14,hjust = 0.5), 
      axis.title.y = element_text(size = 14),
      plot.title = element_text( size=14,hjust = 0.5),
      #legend.title = element_blank(), 
      #legend.position = "none",
      panel.grid.major = element_blank(), 
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      axis.line = element_line(colour = "black")
    )
  
  if(plot){
    
    print(p1)
  }
  
  
  return(p1)
  
}

