mHAG_func <- function(X_scaled, Y_scaled, alpha = 0.5, num_cores = 50,
                      maxit = 50000, lambda_choice = "lambda.1se") {
  
  # Identify constant columns in Y_scaled
  is_constant_column <- apply(Y_scaled, 2, function(col) length(unique(col)) == 3)
  
  # Remove constant columns from Y_scaled
  Y_scaled <- Y_scaled[, !is_constant_column]
  
  # Initialize matrix to store the coefficients; rows correspond to features (columns of X), columns to genes (columns of Y)
  gene_results_matrix <- matrix(0, nrow = ncol(X_scaled), ncol = ncol(Y_scaled))
  rownames(gene_results_matrix) <- colnames(X_scaled)  # Set row names to match X's column names
  colnames(gene_results_matrix) <- colnames(Y_scaled)  # Set column names to match Y's column names
  
  # Setup parallel processing for cv.glmnet
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  on.exit(stopCluster(cl))  # Ensure the cluster is always stopped, even if an error occurs
  
  tryCatch({
    for (gene_idx in seq_len(ncol(Y_scaled))) {
      Y_gene <- Y_scaled[, gene_idx, drop = FALSE]
      
      # Check if Y_gene has constant values
      if (length(unique(Y_gene)) == 1) {
        # Skip constant Y_gene, and leave coefficients as zeros in the result matrix
        cat(paste("Skipping gene", gene_idx, "because Y is constant.\n"))
        next
      }
      
      # Perform cross-validation in parallel
      cv_model <- cv.glmnet(X_scaled, Y_gene, alpha = alpha, nfolds  = 5,
                            family = "gaussian", parallel = TRUE, type.measure = "mse")
      
      # Select the appropriate lambda based on the input parameter
      if (lambda_choice == "lambda.min") {
        best_lambda <- cv_model$lambda.min
      } else if (lambda_choice == "lambda.1se") {
        best_lambda <- cv_model$lambda.1se
      } else {
        stop("Invalid lambda_choice. Use 'lambda.min' or 'lambda.1se'.")
      }
      
      # Fit the final model using the selected lambda
      final_model <- glmnet(X_scaled, Y_gene, alpha = alpha, 
                            lambda = best_lambda, family = "gaussian", maxit = maxit)
      
      # Extract coefficients from the model (excluding intercept)
      coef_final <- as.matrix(coef(final_model, s = best_lambda))[-1, , drop = FALSE]  # Remove intercept
      
      # Store coefficients in the corresponding column of the result matrix
      gene_results_matrix[, gene_idx] <- coef_final
    }
    
    return(gene_results_matrix)  # Return the results matrix
  }, error = function(e) {
    print(paste("Error in model processing: ", e$message))
  })
}


get_non_zero_gene_names <- function(mHAG_est) {
  # Identify rows (genes) that have non-zero values
  non_zero_rows <- apply(mHAG_est, 2, function(col) any(col != 0))
  
  # Get the gene names for the rows that contain non-zero values
  non_zero_gene_names <- colnames(mHAG_est)[non_zero_rows]
  
  # Return the gene names
  return(non_zero_gene_names)
}

mHAG_module <- function(mHAG_est, threshold = 1e-5) {
  
  # Create a normalized copy of mHAG_est and set very small values to zero
  mHAG_est <- mHAG_est
  mHAG_est[abs(mHAG_est) < threshold] <- 0
  
  # Extract names of non-zero genes
  mHVG <- get_non_zero_gene_names(mHAG_est)
  
  # Transpose the normalized matrix and subset to only mHVG columns
  coeff <- t(mHAG_est[, mHVG])
  
  # Identify coefficients greater than zero
  ctHVG <- abs(coeff) > 0
  
  # Return relevant results
  return(list(mHVG = mHVG, coeff = coeff, ctHVG = ctHVG))
}

pval.srt.dot = function(spot_exp, spot_loc, spot_prop,
                        jointpcc = Joinscores[[i]], size = 3, seck = 4){
  
  genes <- rownames(jointpcc) 
  ordered_genes <- genes[order(jointpcc, decreasing = TRUE)] 
  
  h1gene = c(ordered_genes[1:seck],
             ordered_genes[(length(ordered_genes)-(seck-1)):length(ordered_genes)])
  
  plots = list()
  for (k in 1:length(h1gene)){
    
    spot_svg =  spot_exp[h1gene[k],]/max(spot_exp[h1gene[k],])
    data_svg = data.frame(x = spot_loc[,1], 
                          y= spot_loc[,2],
                          exp  = spot_svg,
                          prop = spot_prop)
    
    data_svg$color <- ifelse(data_svg$prop == 0, "gray", NA)  
    
    plots[[k]] = ggplot(data = data_svg, aes(x = x, y = y)) +
      geom_point(aes(color = exp), size = size, data = subset(data_svg, prop != 0)) +
      geom_point(color = "gray", size = size, data = subset(data_svg, prop == 0)) +
      scale_color_viridis() +
      labs(y = "", x = " ", title = as.character(h1gene[k]), color = "", size = "") +
      theme_classic() +
      theme(
        axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.border = element_blank(),
        axis.line = element_blank(),
        plot.title = element_text(size = 8, hjust = 0.5)
      )
  }
  
  return(plots)
}

plot_scatterpie <- function(spot_loc_raw, prop_Sec, use_color, 
                            pie_r = 0.5, circle_r = 0.7, point_size = 6.75) {
  
  # Ensure the rownames of prop_Sec exist in spot_loc_raw
  matching_spots <- intersect(rownames(prop_Sec), rownames(spot_loc_raw))
  
  if (length(matching_spots) == 0) {
    stop("No matching spots found between 'spot_loc_raw' and 'prop_Sec'.")
  }
  
  # Subset spot_loc_raw and prop_Sec to only include matching spots
  pie_data <- cbind(spot_loc_raw[matching_spots, ], prop_Sec[matching_spots, ])
  
  # Base plot with gray points for all spots, and larger circle size for the base
  base_plot <- ggplot(spot_loc_raw, aes(x = x, y = y)) +
    geom_point(shape = 21, fill = "grey90", color = "black", size = circle_r * point_size) + 
    coord_equal() +
    theme_minimal() +
    theme(axis.ticks = element_blank(), 
          legend.title = element_blank(), 
          axis.text = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_blank(),
          axis.title = element_blank())
  
  # Overlay pie charts for selected spots in prop_Sec, with pie charts smaller than circles
  pie_plot <- base_plot +
    geom_scatterpie(aes(x = x, y = y), data = pie_data, 
                    cols = colnames(prop_Sec), pie_scale = pie_r, color = NA) +
    scale_fill_manual(values = use_color$col_vector,
                      breaks = use_color$cell_types)
  
  return(pie_plot)
}


PCC.score = function(spot_exp = Y_input, coeff = coeff,
                     spot_prop, method = "pearson") {
  
  
  cell.type = colnames(coeff)
  mhvg = rownames(coeff)
  proall = as.data.frame(spot_prop[, cell.type])
  
  joint.score = matrix(list(),length(cell.type),1)
  names(joint.score) = cell.type
  
  for (i in cell.type) {
    
    genes = mhvg[coeff[, i] != 0]
    
    # Ensure that genes exist
    if (length(genes) > 0) {
      x = as.matrix(spot_exp[, genes])
      y = proall[, i]
      
      # Compute correlation
      corr = cor(x, y, method = method)
      rownames(corr) = genes
      
      joint.score[[i]] = corr
    } else {
      joint.score[[i]] = NULL  # No significant genes for this cell type
    }
  }
  
  return(joint.score)
}

Jointscore.barplot = function(Score = Score, y.title="Count", x.title="Joint score",
                              plot = TRUE, nrow = 1, fill = "#56B4E9", color = "black") {
  
  cell.type = names(Score)
  plots = list()
  
  for(i in cell.type) {
    
    # Check if the length of Score[[i]] is less than 10
    if (length(Score[[i]]) < 10) {
      next  # Skip to the next iteration if length is less than 10
    }
    
    type_two_gene = data.frame(cor = Score[[i]])
    
    plots[[i]] = ggplot(type_two_gene, aes(x= cor)) +
      geom_histogram(fill = fill ,bins = 40, color = color) + theme_bw() +
      labs(y = y.title, x = x.title, title = i) +
      theme(
        axis.text.x.bottom = element_text(size = 12, hjust = 1), 
        axis.text.y.left = element_text(size = 12),
        axis.title.x = element_text(size = 14, hjust = 0.5), 
        axis.title.y = element_text(size = 14),
        plot.title = element_text(size = 14, hjust = 0.5),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        axis.line = element_line(colour = "black"))
    
  }
  
  if(plot && length(plots) > 0) {
    pfig = plot_grid(plotlist = plots, nrow = nrow)
    print(pfig)
  }
  
  return(plots)
}



cellchat_ana = function(counts = counts, raw = T, min.cells = 20,cutoff = 0.01,
                        label = label,  CellChatDB.use =  CellChatDB.human){
  
  if(raw){
    
    data.input = normalizeData(counts)
    
  }else{
    counts[counts<cutoff] = 0
    data.input = counts
  }
  
  meta <- data.frame(group = label, row.names = names(label))
  
  cellchat <- createCellChat(object = data.input, meta = meta, group.by = "group")
  cellchat <- addMeta(cellchat, meta = meta, meta.name = "labels")
  cellchat <- setIdent(cellchat, ident.use = "labels")  #
  
  cellchat@DB <- CellChatDB.use
  cellchat <- subsetData(cellchat)
  
  cellchat <- identifyOverExpressedGenes(cellchat)
  cellchat <- identifyOverExpressedInteractions(cellchat)
  
  cellchat <- computeCommunProb(cellchat)
  cellchat <- filterCommunication(cellchat, min.cells =  min.cells)
  cellchat <- aggregateNet(cellchat)
  
  return(cellchat)
}



seurat_cluster <- function(counts, min_cells = 20, min_features = 0, 
                           nfeatures = 3000, resolution = 0.2, dims = 1:30,
                           k_param = 30) {
  # Create Seurat object
  pbmc <- CreateSeuratObject(counts = counts, min.cells = min_cells, min.features = min_features)
  
  # Normalize the data
  pbmc <- NormalizeData(pbmc, normalization.method = "LogNormalize")
  
  # Find variable features
  pbmc <- FindVariableFeatures(pbmc, nfeatures = nfeatures)
  
  # Scale the data based on variable features
  all.genes <- VariableFeatures(object = pbmc)
  pbmc <- ScaleData(pbmc, features = all.genes)
  
  # Perform PCA
  pbmc <- RunPCA(pbmc, features = VariableFeatures(object = pbmc))
  
  # Perform UMAP
  pbmc <- RunUMAP(pbmc, dims = dims)
  
  # Find neighbors
  pbmc <- FindNeighbors(pbmc, dims = dims, k.param = k_param)
  
  # Perform clustering
  pbmc <- FindClusters(pbmc, resolution = resolution)
  
  # Return the final Seurat object with clusters
  return(pbmc)
}
## module of wcor
calculateSpatialSimMatrix <- function(st_exp, w) {
  
  # Compute the similarity matrix using the CellTypeSpatialCrossCor function
  sim_mat <- CellTypeSpatialCrossCor(gexpA = t(st_exp), gexpB = t(st_exp), weight = w)
  WCor = sim_mat
  
  # Normalize the similarity matrix
  sim_mat <- sim_mat - min(sim_mat, na.rm = TRUE)
  sim_mat <- sim_mat / max(sim_mat, na.rm = TRUE)
  
  # Make the similarity matrix symmetric
  sim_mat <- (sim_mat + t(sim_mat)) / 2
  
  # Convert to a distance matrix
  sim_mat <- 1 - sim_mat
  
  return(list(WCor = WCor, sim_mat = sim_mat))
}

# Helper function for computing spatial cross-correlation
CellTypeSpatialCrossCor <- function(gexpA, gexpB, weight) {
  n1 <- nrow(gexpA)
  n2 <- nrow(gexpB)
  colmeanA <- colMeans(gexpA)
  colmeanB <- colMeans(gexpB)
  
  N <- n1 + n2
  W <- 2 * sum(weight)
  
  x <- gexpA - matrix(rep(colmeanA, n1), n1, byrow = TRUE)
  y <- gexpB - matrix(rep(colmeanB, n2), n2, byrow = TRUE)
  
  cv1 <- t(x) %*% weight %*% y
  cv2 <- t(x) %*% t(weight) %*% y
  cv <- (cv1 + cv2) / 2
  
  v <- sqrt(colSums(x^2)) %o% sqrt(colSums(y^2))
  
  iSCI <- (N / W) * (cv / v)
  
  return(iSCI)
}


