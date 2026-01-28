library(Seurat)
library(segmented)
library(ggplot2)
library(dplyr)
library(RANN)
library(igraph)
library(rpart)
library(future)
library(future.apply)
library(parallelly)




#' Faster UMAP implementation using subset-and-project approach
#'
#' @param seurat_obj The Seurat object
#' @param subset_size Number or proportion of cells to use for UMAP model training
#' @param n.components Number of UMAP dimensions
#' @param n.neighbors Number of neighbors for UMAP
#' @param random_seed Random seed for reproducibility
#' @param features Features to use (defaults to variable features)
#' @param metric Distance metric to use
#' @param min.dist Minimum distance parameter for UMAP
#' @param verbose Print progress messages
#' @return Updated Seurat object with UMAP embedding
FastRunUMAP <- function(seurat_obj,
                        subset_size = 50000,
                        n.components = 2,
                        n.neighbors = 15,
                        random_seed = 42,
                        features = NULL,
                        metric = "correlation",
                        min.dist = 0.1,
                        verbose = TRUE) {
  require(Seurat)
  require(uwot)
  require(Matrix)
  require(future)
  require(future.apply)
  
  # Set up parallel processing
  if (verbose) message("Preparing data for UMAP...")
  
  # Get features to use
  if (is.null(features)) {
    features <- VariableFeatures(seurat_obj)
  }
  
  # Get dimensionality reduction to use as input (PCA by default)
  if ("pca" %in% names(seurat_obj@reductions)) {
    use_dr <- "pca"
    input_data <- Embeddings(seurat_obj, reduction = "pca")
  } else {
    if (verbose) message("PCA not found, using scaled data...")
    input_data <- t(GetAssayData(seurat_obj, slot = "scale.data")[features, ])
  }
  
  # Determine subset size (if fraction, convert to number)
  cell_count <- ncol(seurat_obj)
  if (subset_size < 1) {
    subset_size <- ceiling(cell_count * subset_size)
  }
  subset_size <- min(subset_size, cell_count)
  
  if (verbose) {
    message(sprintf(
      "Using %d cells (out of %d) to train UMAP model...",
      subset_size, cell_count
    ))
  }
  
  # Randomly select subset of cells for UMAP model training
  set.seed(random_seed)
  subset_idx <- sample(cell_count, size = subset_size)
  
  # Train UMAP on subset
  if (verbose) message("Training UMAP model on subset...")
  
  subset_data <- input_data[subset_idx, ]
  
  # Train UMAP with uwot directly for more speed control
  umap_model <- uwot::umap(
    X = subset_data,
    n_components = n.components,
    n_neighbors = n.neighbors,
    metric = metric,
    min_dist = min.dist,
    ret_model = TRUE,
    n_epochs = 200, # Reduced for speed
    n_threads = future::availableCores() - 1,
    n_sgd_threads = "auto",
    batch = TRUE, # Use batched processing
    verbose = verbose
  )
  
  if (verbose) message("Projecting all cells onto UMAP model...")
  
  # Project all cells onto trained UMAP model
  if (cell_count <= 100000) {
    # For datasets that fit in memory, do it all at once
    all_umap <- uwot::umap_transform(
      X = input_data,
      model = umap_model,
      n_threads = future::availableCores() - 1,
      n_sgd_threads = "auto",
      verbose = verbose
    )
  } else {
    # For very large datasets, process in chunks
    chunk_size <- 100000
    n_chunks <- ceiling(cell_count / chunk_size)
    all_umap <- matrix(0, nrow = cell_count, ncol = n.components)
    
    for (i in 1:n_chunks) {
      if (verbose) message(sprintf("Processing chunk %d of %d...", i, n_chunks))
      start_idx <- (i - 1) * chunk_size + 1
      end_idx <- min(i * chunk_size, cell_count)
      chunk_data <- input_data[start_idx:end_idx, ]
      
      all_umap[start_idx:end_idx, ] <- uwot::umap_transform(
        X = chunk_data,
        model = umap_model,
        n_threads = future::availableCores() - 1,
        n_sgd_threads = "auto",
        verbose = TRUE
      )
    }
  }
  
  # Create new dimensional reduction object
  colnames(all_umap) <- paste0("UMAP_", 1:n.components)
  
  rownames(all_umap) <- colnames(seurat_obj)
  
  umap_dimred <- CreateDimReducObject(
    embeddings = all_umap,
    key = "UMAP_",
    assay = DefaultAssay(seurat_obj)
  )
  
  # Add UMAP to Seurat object
  seurat_obj[["umap"]] <- umap_dimred
  
  if (verbose) message("UMAP computation complete.")
  return(seurat_obj)
}





threshold_groups <- function(input_vector, thresholds) {
  # Initialize an empty vector to hold the group assignments
  groups <- numeric(length(input_vector))
  
  # Loop through each threshold in the list
  for (i in seq_along(thresholds)) {
    # Create a logical vector indicating which elements in the input_vector are above the threshold
    group <- input_vector >= thresholds[i]
    
    # Assign group i to all elements in input_vector that are above threshold i and below threshold (i+1)
    if (i < length(thresholds)) {
      next_group <- input_vector >= thresholds[i + 1]
      group <- group & !next_group
    }
    
    # Assign group i to all elements in input_vector that are above threshold i
    groups[group] <- i
  }
  
  # Return the vector of group assignments
  return(groups)
}




threshold_function <- function(k, data_anb, clusters) {
  ## Save the name of antibodies
  name <- colnames(data_anb[k])
  
  ## add louvain cluster to data intensity
  data_intensity_cluster <- data_frame(data_anb, lv_clusters = clusters)
  
  # Order louvain cluster by median
  group_ordered <- with(
    data_intensity_cluster,
    reorder(
      lv_clusters,
      get(name), median
    )
  )
  
  # Create data with reordered group levels
  data_ordered <- data_intensity_cluster # keep order of orginial data
  data_ordered$group <- factor(data_ordered$lv_clusters,
                               levels = levels(group_ordered)
  )
  
  # Create data with median value in each cluster
  data_median <- data.frame(median_value = attr(group_ordered, "scores"), lv_clusters = names(attr(group_ordered, "scores")))
  data_median <- data_median[order(data_median$median_value, decreasing = F), ]
  
  
  if (sum(data_median$median_value) == 0) {
    threshold_point_final <- max(data_ordered[, name]) # if all median louvain = 0, we do not have positive cells (threshold at max)
  } else {
    x <- 1:nrow(data_median)
    y <- data_median$median_value
    fit <- lm(y ~ x)
    
    # fit piecewise regression model to original model, estimating a breakpoint
    AIC_value <- NULL
    for (i in 1:3) {
      out <- tryCatch(
        {
          segmented.fit <- segmented(fit, seg.Z = ~x, npsi = i)
        },
        error = function(e) {
          e
        }
      )
      if (any(class(out) == "error") == T) {
        AIC_value[i] <- 100000
      } else {
        segmented.fit <- segmented(fit, seg.Z = ~x, npsi = i)
        AIC_value[i] <- AIC(segmented.fit)
      }
    }
    
    ## Choose the smallest AIC for better fit
    segmented.fit <- segmented(fit, seg.Z = ~x, npsi = which.min(AIC_value))
    
    ## get the breakpoint value
    breakpoint_point <- data.frame(x = segmented.fit$psi[, 2], y = segmented.fit$fitted.values[floor(segmented.fit$psi[, 2])])
    
    
    ## assign group base on multiple breakpoints
    Group_expression <- threshold_groups(thresholds = breakpoint_point$x, input_vector = 1:nrow(data_median))
    data_median_update <- data.frame(data_median, Group_expression)
    new_data_anb <- merge(data_intensity_cluster, data_median_update, "lv_clusters")
    
    
    
    
    
    
    ## assign lowest intensity group (0) and highest intensity group (4)
    if (which.min(AIC_value) == 3) {
      new_data_anb$Group_expression <- ifelse(new_data_anb$Group_expression == 3, 4, new_data_anb$Group_expression)
      slope <- as.numeric(segmented.fit$coefficients[c(2, 3, 4, 5)])
      if (slope[2] < 0.05) {
        new_data_anb$Group_expression <- ifelse(new_data_anb$Group_expression == 1, 0, new_data_anb$Group_expression)
      } else {
        new_data_anb$Group_expression <- new_data_anb$Group_expression
      }
    }
    
    if (which.min(AIC_value) == 2) {
      new_data_anb$Group_expression <- ifelse(new_data_anb$Group_expression == 2, 4, new_data_anb$Group_expression)
      slope <- as.numeric(segmented.fit$coefficients[c(2, 3, 4)])
      if (slope[2] < 0.05) {
        new_data_anb$Group_expression <- ifelse(new_data_anb$Group_expression == 1, 0, new_data_anb$Group_expression)
      } else {
        new_data_anb$Group_expression <- new_data_anb$Group_expression
      }
    }
    
    if (which.min(AIC_value) == 1) {
      new_data_anb$Group_expression <- ifelse(new_data_anb$Group_expression == 1, 4, new_data_anb$Group_expression)
      slope <- as.numeric(segmented.fit$coefficients[c(2, 3)])
      if (slope[2] < 0.05) {
        new_data_anb$Group_expression <- ifelse(new_data_anb$Group_expression == 1, 0, new_data_anb$Group_expression)
      } else {
        new_data_anb$Group_expression <- new_data_anb$Group_expression
      }
    }
    
    
    
    
    
    ## create a vector of threshold
    cut_seq <- seq(as.numeric(quantile(new_data_anb[, k + 1], 0.05)), as.numeric(quantile(new_data_anb[, k + 1], 0.99)), length.out = 1000)
    cut_seq <- cut_seq[-1] # sometime it 0, I would remove
    threshold_point_final <- NULL
    low_high <- high_low <- high_high <- low_low <- NULL
    for (g in 1:length(cut_seq)) {
      high_data <- new_data_anb[which(new_data_anb$Group_expression == 4), k + 1]
      low_data <- new_data_anb[which(new_data_anb$Group_expression == 0), k + 1]
      
      # low_high: low expr in high group
      # high_low: high expr in low group
      
      low_high[g] <- length(high_data[which(high_data < cut_seq[g] | high_data == cut_seq[g])]) / length(high_data)
      high_low[g] <- length(low_data[which(low_data > cut_seq[g] | low_data == cut_seq[g])]) / length(low_data)
      
      # high_high: high expr in high group
      # low_low: low expr in low group
      
      high_high[g] <- length(high_data[which(high_data > cut_seq[g] | high_data == cut_seq[g])]) / length(high_data)
      low_low[g] <- length(low_data[which(low_data < cut_seq[g] | low_data == cut_seq[g])]) / length(low_data)
    }
    threshold_point_final <- cut_seq[which.min(low_high + high_low)]
  }
  
  ## create matrix binarization
  Group_Threshold_new <- ifelse(as.numeric(t(data_intensity_cluster[, name])) > threshold_point_final, 1, 0)
  
  
  return(list(Group_Threshold_new, threshold_point_final))
}

TACIT <- function(data_expression, r, p, Signature) {
  

  data_anb=data_expression[,colnames(Signature)]

  ### Processing data
  #----
  ## Signature data (St)
  Signature[is.na(Signature) == T] <- 0 ## assign all NA as 0
  Signature <- as.data.frame(Signature)

  # Intensity data (A_ij)

  main_cell_types_with_markers <- find_main_cell_types_with_subsets_and_markers(Signature)
  ct_main <- names(main_cell_types_with_markers)
  marker_subset <- unique(unlist(main_cell_types_with_markers))



  orig_values <- as.matrix(data_anb)
  rownames(orig_values) <- 1:nrow(data_anb)
  orig_values <- t(orig_values)
  orig_values_metadata <- data.frame("CellID" = 1:nrow(data_anb))
  row.names(orig_values_metadata) <- orig_values_metadata$CellID
  scfp <- CreateSeuratObject(counts = orig_values, meta.data = orig_values_metadata)


  ### Create Seurat object
  #----
  print("-------------------REACHED HERE: NormalizeData---------------")
  # Set up parallel processing
  plan(sequential)
  scfp <- NormalizeData(scfp, normalization.method = "CLR", margin = 2) # do not normalize data
  print("-------------------REACHED HERE: FindVariableFeatures---------------")
  scfp <- FindVariableFeatures(scfp, selection.method = "disp", nfeatures = 500)
  print("-------------------REACHED HERE: ScaleData---------------")
  plan(multisession, workers = parallelly::availableCores() - 1)
  scfp <- ScaleData(scfp, features = rownames(scfp))





  #### Clustering
  #----
  print("-------------------REACHED HERE: RunPCA---------------")
  scfp <- RunPCA(scfp, features = VariableFeatures(object = scfp))
  print("-------------------REACHED HERE: RunUMAP---------------")


  ####Clustering
  #----
  scfp <- RunPCA(scfp, features = VariableFeatures(object = scfp))
  print("-------------------REACHED HERE: RunUMAP---------------")
  total_cells <- ncol(scfp)
  #### Clustering
  #----
  if (total_cells < 100000) {
    print("Less than 100000 cells, running UMAP on the entire dataset")
    scfp <- RunUMAP(scfp,
                    features = VariableFeatures(object = scfp), n.components = p, metric = "correlation", umap.method = "uwot",
                    n.neighbors = 15, learning.rate = 1, spread = 1, min.dist = 0.01, set.op.mix.ratio = 1, local.connectivity = 1
    )
  } else {
    print("More than 100000 cells, training UMAP using 30% of cells for training or max 100000 cells, whichever is smaller")
    # Use 30% of cells for training or max 100000 cells, whichever is smaller
    subset_size <- min(100000, ceiling(total_cells * 0.3))
    # For very small datasets, use all cells
    subset_size <- min(subset_size, total_cells)

    scfp <- FastRunUMAP(
      seurat_obj = scfp,
      subset_size = subset_size,
      n.components = p,
      n.neighbors = 15,
      metric = "correlation",
      min.dist = 0.01, # Same as original
      random_seed = 1, # For reproducibility
      verbose = TRUE
    )
  }

  print("-------------------REACHED HERE: FindNeighbors---------------")
  scfp <- FindNeighbors(scfp, dims = 1:p, reduction = "umap")




  # Run clusters
  print("-------------------REACHED HERE: FindClusters---------------")
  scfp <- FindClusters(scfp, resolution = r, random.seed = 1, algorithm = 4, n.iter = 2)
  print("---------------- now here ---------------")
  # Collected clusters ID for each cells
  clusters <- as.numeric(scfp@meta.data[["seurat_clusters"]])

  print(paste0("Number of clusters: ", length(unique(clusters)), " with ", round(mean(as.numeric(table(clusters))), 0), " cells per clusters"))

  ### Cell type score
  #----
  data_anb_matrix <- as.matrix((data_anb)) # scale data intensity and convert to matrix from data frame
  Signature_matrix <- as.matrix(Signature) # convert signature to matrix from data frame
  ct <- data_anb_matrix %*% t(Signature_matrix)
  colnames(ct) <- Signature


  #  data_anb=scale(data_anb)

  ### Threshold for cell type
  #----
  library(parallel)
  library(foreach)
  library(doParallel)
  library(caret)

  ct <- as.data.frame(ct)

  num_cores <- parallelly::availableCores() - 1 # Leave one core free for system processes
  registerDoParallel(cores = num_cores)

  # Parallelize first loop
  Group_threshold_data <- NULL
  final_threshold <- NULL

  cat("==== Phase 1: Processing Main Data Columns ====\n")
  cat("Processing", ncol(ct), "columns in parallel...\n")
  print('Entering first loop...')
  final_threshold <- foreach(
    k = 1:ncol(ct), .combine = c,
    .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
  ) %dopar% {
    # Return the second element of the threshold_function result
    vector <- threshold_function(k, ct, clusters)
    return(vector[[2]])
  }

  print('entering second loop...')
  # Process Group_threshold_data in parallel
  Group_threshold_data <- foreach(
    k = 1:ncol(ct), .combine = cbind,
    .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
  ) %dopar% {
    # Return the first element of the threshold_function result
    vector <- threshold_function(k, ct, clusters)
    return(vector[[1]])
  }


  colnames(Group_threshold_data) <- colnames(ct)
  final_threshold <- data.frame(value = final_threshold, Name = colnames(ct)) ########### . Output 1
  # 
  # # Conditionally handle marker subset processing
  # if (length(marker_subset) > 1) {
  #   cat("==== Phase 2: Processing Marker Subset Data ====\n")
  #   cat("Detected", length(marker_subset), "marker subsets to process\n")
  #   
  #   ct2 <- data_anb[, marker_subset]
  #   ct2 <- as.data.frame(ct2)
  #   colnames(ct2) <- marker_subset
  #   
  #   cat("Processing marker subset columns in parallel...\n")
  #   
  #   # Parallelize second loop
  #   final_threshold2 <- foreach(
  #     k = 1:ncol(ct2), .combine = c,
  #     .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
  #   ) %dopar% {
  #     vector <- threshold_function(k, data_anb = ct2, clusters = clusters)
  #     return(vector[[2]])
  #   }
  #   
  #   Group_threshold_data2 <- foreach(
  #     k = 1:ncol(ct2), .combine = cbind,
  #     .packages = c("stats", "segmented", "dplyr"), .export = c("threshold_function", "threshold_groups")
  #   ) %dopar% {
  #     vector <- threshold_function(k, data_anb = ct2, clusters = clusters)
  #     return(vector[[1]])
  #   }
  #   
  #   final_threshold2 <- data.frame(value = final_threshold2, Name = marker_subset) ########### . Output 1
  #   colnames(Group_threshold_data2) <- marker_subset
  #   
  #   cat("==== Phase 3: Processing Cell Types and Markers ====\n")
  #   cat("Processing", length(ct_main), "main cell types...\n")
  #   
  #   # Vectorize operations where possible
  #   for (p in 1:length(ct_main)) {
  #     cat("  - Processing cell type", p, "of", length(ct_main), ":", ct_main[p], "\n")
  #     ct_score_current <- Group_threshold_data[, ct_main[p]]
  #     
  #     if (length(unlist(main_cell_types_with_markers[ct_main[p]])) == 1) {
  #       rSum <- (Group_threshold_data2[, unlist(main_cell_types_with_markers[ct_main[p]])])
  #     } else {
  #       rSum <- rowSums(Group_threshold_data2[, unlist(main_cell_types_with_markers[ct_main[p]])])
  #     }
  #     Group_threshold_data[, ct_main[p]] <- ifelse(rSum > 0, 0, ct_score_current)
  #   }
  #   
  #   unique_markers_for_subsets <- find_unique_marker_for_subsets(Signature)
  #   if (length(unique_markers_for_subsets) > 0) {
  #     cat("\nProcessing", length(unique_markers_for_subsets), "unique markers for subsets...\n")
  #     # Process unique markers more efficiently
  #     for (p in 1:length(unique_markers_for_subsets)) {
  #       ctype <- names(unique_markers_for_subsets[p])
  #       marker_id <- unique_markers_for_subsets[ctype][[1]]
  #       Group_threshold_data[, ctype] <- ifelse(Group_threshold_data2[, marker_id] == 0, 0, Group_threshold_data[, ctype])
  #     }
  #   }
  # }
  # 
  # cat("==== Phase 4: Post-Processing and Analysis ====\n")
  # cat("Generating potential cell types...\n")
  # potential_ct <- replace_zero_with_column_names_opt(Group_threshold_data)
  # 
  # cat("Creating test data...\n")
  # potential_ct_test <- replace_zero_with_column_names_opt(Group_threshold_data)
  # combined_vector <- apply(potential_ct_test, 1, function(x) paste(x[which(x != "0")], collapse = "::"))
  # 
  # cat("Calculating group mixtures...\n")
  # group_mix <- rowSums(Group_threshold_data)
  # Group_threshold_data <- data.frame(ID = 1:nrow(Group_threshold_data), clusters, Group_threshold_data)
  # data_cluster <- data.frame(clusters, ct)
  # colnames(data_cluster)[-1] <- Signature$cell_type
  # 
  # cat("Calculating median values per cluster...\n")
  # medians_per_cluster <- aggregate(. ~ clusters, data = data_cluster, FUN = median)
  # medians_per_cluster_filter <- medians_per_cluster[, -1]
  # 
  # cat("Applying thresholds and converting to binary values...\n")
  # # Convert final_threshold to a named vector for easier lookup
  # thresholds <- setNames(final_threshold$value, final_threshold$Name)
  # # Apply thresholding to convert values to binary
  # binary_df <- medians_per_cluster_filter
  # binary_df[] <- sapply(names(medians_per_cluster_filter), function(column_name) {
  #   if (column_name %in% names(thresholds)) {
  #     as.integer(medians_per_cluster_filter[[column_name]] > thresholds[[column_name]])
  #   } else {
  #     medians_per_cluster_filter[[column_name]]
  #   }
  # })
  # 
  # cat("Generating potential cell type clusters...\n")
  # potential_ct_cluster <- replace_zero_with_column_names_opt(binary_df)
  # combined_vector <- apply(potential_ct_cluster, 1, function(x) paste(x[which(x != "0")], collapse = "::"))
  # 
  # cat("Identifying clean and mixed clusters...\n")
  # cluster_clean <- medians_per_cluster$clusters[rowSums(binary_df) == 1]
  # cluster_mix <- medians_per_cluster$clusters[rowSums(binary_df) != 1]
  # 
  # data_cluster_group <- data.frame(clusters = medians_per_cluster$clusters, combined_vector)
  # data_merge_cluster <- merge(data_cluster_group, Group_threshold_data, "clusters")
  # data_merge_cluster <- data_merge_cluster[order(data_merge_cluster$ID), ]
  # 
  # 
  # cat("Separating clean and mixed groups...\n")
  # data_clean_group <- data_merge_cluster[which(data_merge_cluster$clusters %in% cluster_clean), ]
  # data_mix_group <- data_merge_cluster[which(data_merge_cluster$clusters %in% cluster_mix), ]
  # 
  # cat("\nCluster Summary:\n")
  # print(paste0("Number of clean MicroClusters: ", length(cluster_clean)))
  # print(paste0("Number of mixed MicroClusters: ", length(cluster_mix)))
  # cat("- Total MicroClusters:", length(cluster_clean) + length(cluster_mix), "\n\n")
  # 
  # colnames(data_mix_group)[-c(1, 2, 3)] <- Signature$cell_type
  # 
  # potential_ct <- replace_zero_with_column_names_opt(data_mix_group[, -c(1, 2, 3)])
  # 
  # ### Identify clean and mix data
  # 
  # cat("Processing cell type combinations...\n")
  # ## vector of potential cell type
  # combined_vector <- apply(potential_ct, 1, function(x) paste(x[which(x != "0")], collapse = "::"))
  # ## adding the vector combination into data
  # 
  # cat("Preparing final data matrices...\n")
  # data_anb_matrix <- data_anb_matrix[, colnames(Signature)[-1]]
  # data_check <- data.frame(combined_vector, data_anb_matrix[data_mix_group$ID, ], data_mix_group[, -c(1, 2, 3)], ID = data_mix_group$ID)
  # 
  # 
  # colnames(data_check)[2:sum(ncol(data_anb_matrix), 1)] <- colnames(data_anb_matrix)
  # colnames(data_check)[sum(ncol(data_anb_matrix), 2):sum(ncol(data_anb_matrix), 1, length(Signature[, 1]))] <- Signature[, 1]
  # 
  # ## Count the number of instances where a cell's value exceeds the corresponding threshold.
  # row_sums <- rowSums(data_check[, sum(ncol(data_anb_matrix), 2):sum(ncol(data_anb_matrix), 1, length(Signature[, 1]))])
  # 
  # ## Selected clean cell types
  # data_out <- data_check[row_sums %in% c(1, 0), ]
  # data_clean <- data.frame(mem = c(data_out$combined_vector), ID = c(data_out$ID))
  # 
  # ## Selected clean cell types
  # data_mixed <- data_check[!(row_sums %in% c(1, 0)), ]
  # data_mixed <- data_mixed[which(data_mixed$combined_vector %in% names(table(data_mixed$combined_vector))[as.numeric(table(data_mixed$combined_vector)) > 0]), ]
  # data_mixed <- data_mixed[order(rowSums(data_mixed[sum(ncol(data_anb_matrix), 2):sum(ncol(data_anb_matrix), 1, length(Signature[, 1]))])), ]
  # 
  # 
  # cat("==== Phase 5: Final Processing and KNN Analysis ====\n")
  # ######## Solving for mixed group
  # ## KNN
  # #----
  # library(parallel)
  # library(foreach)
  # library(doParallel)
  # library(caret)
  # 
  # data_out_test <- data_out
  # data_clean_final <- data_clean
  # data_clean_cluster <- data.frame(mem = data_clean_group$combined_vector, ID = data_clean_group$ID)
  # data_clean_final <- rbind(data_clean_final, data_clean_cluster)
  # 
  # 
  # print(paste0("Number of clean cells: ", nrow(data_clean_final), " with ", length(unique(data_clean_final$mem)), " cell types"))
  # print(paste0("Number of mixed cells: ", nrow(data_mixed), " with ", length(unique(data_mixed$combined_vector)), " mixed cell types"))
  # 
  # # Get unique mixed cell types
  # mix_group <- unique(data_mixed$combined_vector)
  # 
  # clean_cell_types <- unique(data_clean_final$mem)
  # mixed_by_type <- split(data_mixed, data_mixed$combined_vector)
  # 
  # # Get signature markers for each clean cell type
  # signature_markers <- list()
  # for (ct in clean_cell_types) {
  #   signature_markers[[ct]] <- names(which(colSums(Signature[Signature[, 1] == ct, ][-1]) > 0))
  # }
  # 
  # num_cores <- min(parallelly::availableCores() - 1, length(mix_group))
  # registerDoParallel(cores = num_cores)
  # 
  # results <- foreach(i = 1:length(mix_group), .combine = rbind, .packages = c("caret")) %dopar% {
  #   current_mix_type <- mix_group[i]
  #   
  #   # Get all cells with this specific mixed type
  #   data_choose_second <- mixed_by_type[[current_mix_type]]
  #   
  #   # Get potential component cell types
  #   mix_select <- strsplit(current_mix_type, "::")[[1]]
  #   mix_select <- intersect(mix_select, clean_cell_types)
  #   
  #   if (length(mix_select) == 0) {
  #     return(data.frame(ID = data_choose_second$ID, mem = "Others"))
  #   }
  #   
  #   # Selected markers that signature for component cell type
  #   anb_select <- names(which(colSums(Signature[Signature[, 1] %in% mix_select, ][-1]) > 0))
  #   
  #   if (length(anb_select) == 0) {
  #     return(data.frame(ID = data_choose_second$ID, mem = "Others"))
  #   }
  #   
  #   # Get the data intensity from the clean cells (all at once for these cell types)
  #   data_clean_work <- data_out_test[
  #     which(data_out_test$combined_vector %in% mix_select),
  #     c(anb_select, "combined_vector", "ID")
  #   ]
  #   
  #   if (nrow(data_clean_work) == 0) {
  #     return(data.frame(ID = data_choose_second$ID, mem = "Others"))
  #   } else {
  #     # Extract relevant data columns for mixed cells
  #     data_test <- data_choose_second[, c(anb_select, "ID")]
  #     
  #     # Add small random noise (vectorized operation)
  #     train_matrix <- as.matrix(data_clean_work[, anb_select]) +
  #       matrix(runif(nrow(data_clean_work) * length(anb_select), 0, 0.001),
  #              nrow = nrow(data_clean_work)
  #       )
  #     
  #     test_matrix <- as.matrix(data_test[, anb_select]) +
  #       matrix(runif(nrow(data_test) * length(anb_select), 0, 0.001),
  #              nrow = nrow(data_test)
  #       )
  #     
  #     knn_model <- caret::knn3(train_matrix,
  #                              as.factor(data_clean_work$combined_vector),
  #                              k = min(10, nrow(data_clean_work))
  #     )
  #     
  #     # Get predictions for all cells at once
  #     predictions <- predict(knn_model, test_matrix, type = "prob")
  #     
  #     max_indices <- max.col(predictions, ties.method = "first")
  #     predicted_types <- colnames(predictions)[max_indices]
  #     
  #     # Return batch results
  #     return(data.frame(ID = data_test$ID, mem = predicted_types))
  #   }
  # }
  # 
  # # Combine results with clean cells
  # data_clean_final <- rbind(data_clean_final, results)
  # 
  # # Final processing
  # data_clean_final$mem <- ifelse(data_clean_final$mem == "", "Others", data_clean_final$mem)
  # data_clean_final <- data_clean_final[order(as.numeric(data_clean_final$ID), decreasing = F), ]
  # 
  # data_final=data.frame(TACIT=data_clean_final$mem,UMAP1=as.numeric(scfp@reductions[["umap"]]@cell.embeddings[,1]),
  #                       UMAP2=as.numeric(scfp@reductions[["umap"]]@cell.embeddings[,2]),data_expression)
  # 
  # return(data_final)
  data_final=data.frame(UMAP1=as.numeric(scfp@reductions[["umap"]]@cell.embeddings[,1]),
                        UMAP2=as.numeric(scfp@reductions[["umap"]]@cell.embeddings[,2]),
                        data_expression)
  return(list(
    df = data_final, 
    threshold = final_threshold)
  )
}




find_main_cell_types_with_subsets_and_markers <- function(signature) {
  # Extract cell types and markers
  cell_types <- signature[, 1]
  markers <- colnames(signature)[-1] # Exclude the first column which is Cell_type
  
  # Initialize a list to store the main cell types and their subset markers
  main_cell_types <- list()
  
  # Loop over each cell type
  for (i in 1:length(cell_types)) {
    main_type <- cell_types[i]
    main_markers <- colnames(signature)[signature[i, ] > 0]
    
    # Store additional markers for subsets
    additional_subset_markers <- character()
    
    # Check against all other cell types for a subset
    for (j in 1:length(cell_types)) {
      if (i != j) {
        subset_type <- cell_types[j]
        subset_markers <- colnames(signature)[signature[j, ] > 0]
        
        # Check if subset has all markers of main and additional ones
        if (all(main_markers %in% subset_markers) && length(subset_markers) > length(main_markers)) {
          additional_subset_markers <- unique(c(additional_subset_markers, setdiff(subset_markers, main_markers)))
        }
      }
    }
    
    # Add to main cell types with their additional subset markers
    if (length(additional_subset_markers) > 0) {
      main_cell_types[[main_type]] <- additional_subset_markers
    }
  }
  
  return(main_cell_types)
}





find_unique_marker_for_subsets <- function(signature) {
  cell_types <- signature[, 1]
  markers <- colnames(signature)[-1] # Exclude the first column which is Cell_type
  
  # Initialize a list to store the unique marker for each subset cell type
  unique_markers <- list()
  
  # Loop over each cell type
  for (i in 1:length(cell_types)) {
    subset_type <- cell_types[i]
    subset_markers <- colnames(signature)[signature[i, ] > 0]
    
    # Find the immediate parent cell type
    parent_markers <- NULL
    for (j in 1:length(cell_types)) {
      if (i != j) {
        potential_parent_markers <- colnames(signature)[signature[j, ] > 0]
        if (all(potential_parent_markers %in% subset_markers) &&
            length(potential_parent_markers) < length(subset_markers)) {
          parent_markers <- potential_parent_markers
        }
      }
    }
    
    # Determine the unique marker
    if (!is.null(parent_markers)) {
      additional_markers <- setdiff(subset_markers, parent_markers)
      if (length(additional_markers) == 1) {
        unique_markers[[subset_type]] <- additional_markers
      }
    }
  }
  
  return(unique_markers)
}




















replace_zero_with_column_names <- function(binary_data) {
  result_matrix <- matrix("", nrow = nrow(binary_data), ncol = ncol(binary_data))
  
  for (i in 1:nrow(binary_data)) {
    for (j in 1:ncol(binary_data)) {
      if (binary_data[i, j] == 1) {
        result_matrix[i, j] <- colnames(binary_data)[j]
      } else {
        result_matrix[i, j] <- 0
      }
    }
  }
  
  return(result_matrix)
}

replace_zero_with_column_names_opt <- function(binary_data) {
  # Create result matrix
  result_matrix <- matrix("0", nrow = nrow(binary_data), ncol = ncol(binary_data))
  col_names <- colnames(binary_data)
  
  # Use matrix operations instead of loops
  for (j in 1:ncol(binary_data)) {
    # Find which rows have 1 in this column
    idx <- which(binary_data[, j] == 1)
    if (length(idx) > 0) {
      result_matrix[idx, j] <- col_names[j]
    }
  }
  
  return(result_matrix)
}



class_accuracy <- function(df) {
  unique_classes <- unique(df$true)
  accuracies <- numeric(length(unique_classes))
  for (i in 1:length(unique_classes)) {
    class <- unique_classes[i]
    correct <- sum(df$true == class & df$predict == class)
    total <- sum(df$true == class)
    accuracies[i] <- correct / total
  }
  data_return <- data.frame(cell_type = unique_classes, accuracies)
  return(data_return)
}











# Custom function to remove connections that exceed the threshold distance
removeFarConnections <- function(tri, cutoff) {
  # Extract the coordinates of the points from the triangulation
  point_coords <- tri$delsgs[, 1:4]
  
  # Compute the Euclidean distance between each pair of points
  distances <- sqrt((point_coords[, 1] - point_coords[, 3])^2 + (point_coords[, 2] - point_coords[, 4])^2)
  threshold_distance <- as.numeric(quantile(as.numeric(distances), cutoff))
  
  # Get the indices of edges that are within the threshold distance
  valid_edges <- distances <= threshold_distance
  
  # Get the indices of points involved in the valid edges
  valid_points <- unique(c(tri$delsgs[valid_edges, 1], tri$delsgs[valid_edges, 2]))
  
  # Create a new triangulation object with the valid connections
  new_tri <- tri
  new_tri$delsgs <- new_tri$delsgs[valid_edges, ]
  new_tri$dirsgs <- new_tri$dirsgs[valid_edges, ]
  new_tri$ptlist <- new_tri$ptlist[valid_points, ]
  
  return(new_tri)
}

library(mclust) # For Rand Index

calculate_metrics <- function(predictions, true_labels) {
  # Ensure both are factors and have the same levels
  combined_levels <- union(levels(factor(predictions)), levels(factor(true_labels)))
  predictions <- factor(predictions, levels = combined_levels)
  true_labels <- factor(true_labels, levels = combined_levels)
  
  # Overall accuracy
  overall_accuracy <- sum(predictions == true_labels) / length(true_labels)
  
  # Weighted accuracy (as provided)
  table_true <- table(true_labels)
  weighted_accuracy <- sum((table(predictions[true_labels == predictions]) / table_true) * (table_true / sum(table_true)), na.rm = TRUE)
  
  # Sensitivity (Recall) for each cell type
  sensitivity <- sapply(levels(true_labels), function(x) {
    true_positives <- sum(predictions[true_labels == x] == x)
    actual_positives <- sum(true_labels == x)
    if (actual_positives == 0) {
      return(NA)
    } # Avoid division by zero
    true_positives / actual_positives
  })
  names(sensitivity) <- levels(true_labels)
  sensitivity[is.na(sensitivity) == T] <- 0
  
  # Precision (Positive Predictive Value) for each cell type
  precision <- sapply(levels(predictions), function(x) {
    true_positives <- sum(predictions[true_labels == x] == x)
    predicted_positives <- sum(predictions == x)
    if (predicted_positives == 0) {
      return(NA)
    } # Avoid division by zero
    true_positives / predicted_positives
  })
  names(precision) <- levels(predictions)
  precision[is.na(precision) == T] <- 0
  
  
  # Calculate weighted recall and precision
  weights <- table_true / sum(table_true)
  weighted_recall <- sum(sensitivity * weights, na.rm = TRUE)
  weighted_precision <- sum(precision * weights, na.rm = TRUE)
  
  # F1 Score for each cell type
  f1_scores <- 2 * (precision * sensitivity) / (precision + sensitivity)
  names(f1_scores) <- levels(predictions)
  f1_scores[is.na(f1_scores) == T] <- 0
  
  weighted_f1 <- sum(f1_scores * weights, na.rm = TRUE)
  
  return(list(
    OverallAccuracy = overall_accuracy,
    WeightedAccuracy = weighted_accuracy,
    Sensitivity = sensitivity,
    Precision = precision,
    WeightedRecall = weighted_recall,
    WeightedPrecision = weighted_precision,
    F1Scores = f1_scores,
    WeightedF1 = weighted_f1
  ))
}





reorder_columns_based_on_signature <- function(Signature) {
  # Initialize the new order of columns, starting with cell_type
  new_col_order <- c("cell_type")
  
  # Loop through each cell type
  for (cell_type in Signature$cell_type) {
    # Find signature markers for this cell type (values == 1)
    signature_markers <- names(Signature)[which(Signature[Signature$cell_type == cell_type, ] > 0)]
    
    # Exclude 'cell_type' column and any already included columns
    signature_markers <- setdiff(signature_markers, c("cell_type", new_col_order))
    
    # Add these markers to the new column order
    new_col_order <- c(new_col_order, signature_markers)
  }
  
  # Add any remaining columns that have not been added yet
  remaining_cols <- setdiff(names(Signature), new_col_order)
  new_col_order <- c(new_col_order, remaining_cols)
  
  # Reorder the dataframe according to the new column order
  Signature_reordered <- Signature[, new_col_order]
  
  return(Signature_reordered)
}
