

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
      
      estres[[i]] = temp
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

TCA_est<- function(spot_exp, beta = beta,C1 = NULL, C2 = NULL,
                   parallel = TRUE, num_cores= 20) {
  
  tca.mdl <- tca(X = sqrt(spot_exp), W = beta, C1 = C1, C2 = C2,
                 parallel = parallel,  num_cores = num_cores)
  
  Z_hat <- tensor(X = sqrt(spot_exp), tca.mdl)
  

  return(Z_hat)
}


ENIGMA_est <- function(spot_exp, sc_exp, sc_label, beta = NULL) {
  
  cell_type <- sort(unique(sc_label))
  egm =  create_ENIGMA(bulk= spot_exp,ref= sc_exp,
                       ref_type = "sort",
                       meta_ref = as.data.frame(sc_label))

  
  # Estimate cell-type proportions if beta is not provided
  if (is.null(beta)) {

    cat("Estimation of cell type proportion\n")

    egm = get_cell_proportion(egm, method = "RLR")

  } else {
    egm@result_cell_proportion =beta
  }
  
  # Perform deconvolution using ENIGMA
  egm = ENIGMA_trace_norm(egm, model_tracker = TRUE,model_name = "log", preprocess = "log")
  ENIGMA_trace.v = sce2array(egm,  norm_output = FALSE,model_name = "log")

  # Extract dimensions from the result
  dims <- dim(ENIGMA_trace.v )
  # Initialize list to store results for each cell type
  exp_ENIGMA_list <- vector("list", length = dims[3])
  names(exp_ENIGMA_list) <- cell_type
  
  # Store each deconvolved matrix in the list
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
  epsilon = 1e-12
  # Safe log function to handle small values and avoid issues
  safe_log <- function(x, base = 2) {
    x[x < epsilon] <- epsilon # Replace very small values with a threshold to avoid log(0)
    return(log(x, base = base))
  }
  
  # Calculate metrics for each column
  calc_metrics <- function(i) {
    P <- F.est.all[, i][beta.ind.all[, i]]
    Q <- F.true.all[, i][beta.ind.all[, i]]
    
    P_norm <- (P + epsilon) / sum(P + epsilon)
    Q_norm <- (Q + epsilon) / sum(Q + epsilon)

    # Handle edge cases
    if (length(P) < 2 || length(Q) < 2 || any(is.na(P)) || any(is.na(Q)) || var(P) == 0 || var(Q) == 0) {
      return(rep(NA, 4))
    }
    
    rmse <- sqrt(mean((P_norm - Q_norm) ^ 2, na.rm = TRUE))
    
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


corr.cal.gene <- function(F.est, F.true, beta.ind) {
  
  # Validate input dimensions
  if (length(F.est) != length(F.true) || length(F.est) != length(beta.ind)) {
    stop("All input lists must have the same length.")
  }
  
  # Combine filtered matrices
  F.est.all <- do.call(cbind, F.est)
  F.true.all <- do.call(cbind, F.true)
  beta.ind.all <- do.call(cbind, beta.ind)
  
  epsilon <- 1e-10  
  safe_log <- function(x, base = 2) {
    x[x < epsilon] <- epsilon
    return(log(x, base = base))
  }
  
  calc_metrics <- function(i) {
    P <- F.est.all[i, beta.ind.all[i, ]]
    Q <- F.true.all[i, beta.ind.all[i, ]]
    
    # Handle edge cases
    if (length(P) < 2 || length(Q) < 2 || any(is.na(P)) || any(is.na(Q)) || var(P) == 0 || var(Q) == 0) {
      return(rep(NA, 4))
    }
    
    rmse <- sqrt(mean((P - Q) ^ 2, na.rm = TRUE))
    pearson_corr <- cor(P, Q, use = "complete.obs")
    cosine_sim <- 1 - abdiv::cosine_distance(P, Q)
    
  
    P_norm <- (P + epsilon) / sum(P + epsilon)
    Q_norm <- (Q + epsilon) / sum(Q + epsilon)

    

    jsd <- tryCatch({
      M <- 0.5 * (P_norm + Q_norm)
      jsd_value <- 0.5 * (sum(P_norm * safe_log(P_norm / M)) + sum(Q_norm * safe_log(Q_norm / M)))
      max(jsd_value, 0)  # Ensure non-negative value
    }, error = function(e) NaN)
    
    return(c(pearson_corr, cosine_sim, rmse, jsd))
  }
  
  # Apply the metrics calculation across rows (genes)
  est.res <- t(vapply(1:nrow(F.est[[1]]), calc_metrics, numeric(4)))
  
  return(est.res)
}
