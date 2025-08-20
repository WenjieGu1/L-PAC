# load dataset
```{r}
# Load the dataset
data <- read.table(file.choose(), header = TRUE, sep = "\t", comment.char = "#")

data_ref <- read.table(file.choose(), header = TRUE, sep = "\t")
# get ids and sample names
data_id <- data_ref[,c('samplenames','patient_id')]

# Extract only the numeric CNA values (exclude chromosome, start, end, feature columns)
bin_val <- as.matrix(data[, 5:ncol(data)])

# transform value
bin_val <- 2^bin_val

# Normalize each column (sample) to have mean ~1
bin_val <- sweep(bin_val, 2, colMeans(bin_val), "/")

# Ensure at least two samples exist
if (ncol(bin_val) < 2) {
  stop("Input should have at least two samples (columns).")
}

```

```{r}
# Load new dataset
data <- read.table(file.choose(), header = TRUE, sep = "\t", comment.char = "#")

data_ref <- read.csv(file.choose(), header = TRUE)
# get ids and sample names
data_id <- data_ref[,c('samplenames','patient_id')]

# Extract only the numeric CNA values (exclude chromosome, start, end, feature columns)
bin_val <- as.matrix(data[, 5:ncol(data)])

# transform value
bin_val <- 2^bin_val

# Normalize each column (sample) to have mean ~1
bin_val <- sweep(bin_val, 2, colMeans(bin_val,na.rm=TRUE), "/")

# Ensure at least two samples exist
if (ncol(bin_val) < 2) {
  stop("Input should have at least two samples (columns).")
}

```

```{r}
# Load dataset 3
data <- read.csv(file.choose(), header = TRUE)

data_ref <- read.csv(file.choose(), header = TRUE)

# get ids and sample names
data_id <- data_ref[,c('samplenames','patient_id')]

```


# functions
```{r}
library(ggplot2)
library(tidyr)

# A main function to run cluster-LPAC and plot profile
## INPUT
#  bin_val     : Nb x Ns size array containing relative CNA values of bins(not normalized)
#  of equal size. With Nb the number of bins, and Ns the number of multi-region samples.
#  seg_val     : Nb x Ns size array containing values of segment length.
#  Cval        : define the cut-off consine similarity value
#  Dmin        : define initial minimum absCNA distance(CNH)
#  plot        : TRUE means plotting absolute copy number profile
#  seg_info    : dataframe that contains chromosome, start, end and feature columns
cluster_lpac_main <- function(bin_val,seg_len=NULL, Cval=0.96, plot=FALSE, seg_info=NULL, img_path){
  result <- cluster_LPAC(bin_val = bin_val, seg_len = seg_len, Cval = Cval)
  if (plot){
    if(is.null(seg_info)){stop(paste("Plot requires seg_info"))}
    plot_profile(bin_val,seg_info,result)
    ggsave(paste0(img_path,".png"), plot = plot_profile(bin_val,seg_info,result),width = 16, height = 5, dpi = 300)
  }
  return(result)
}


# Define a function to calculate weighted cosine similarity
cosine_similarity <- function(a, b, w) {
  sum(a * b * w) / (sqrt(sum(a * a * w)) * sqrt(sum(b * b * w)))
}

# Define a function to calculate weighted standard deviation
sdr <- function(val,len){
  # Calculate weighted mean
  weighted_mean <- sum(val * len) / sum(len)
  
  # Calculate weighted variance
  weighted_var <- sum(len * (val - weighted_mean)^2) / sum(len)
  
  # Return weighted standard deviation
  return(sqrt(weighted_var))
}

# Convert relative copy number to absolute copy number
rcn_to_acn <- function(value, purity, ploidy){
  a1 <- (purity*ploidy+2*(1-purity))/purity
  a2 <- 2*(1-purity)/purity
  q <- a1*value - a2
  return(q)
}

  
# Plot copy number profile
#  bin_val     : same to the cluster-LPAC input
#  seg_info    : a dataframe contains chromosome, start, end and feature columns
#  result      : cluster-LPAC output result
plot_profile <- function(bin_val,seg_info,result){
  acn_val <- bin_val
  for (i in ncol(bin_val)){
    acn_val[,i] <- rcn_to_acn(bin_val[,i],result$purities[i],result$ploidies[i])
  }
  profile_df <- cbind(seg_info, acn_val)
  return(plot_genome_cn_multi(profile_df,sample_cols = colnames(bin_val)))
}


# copy number profile for multi-samples
plot_genome_cn_multi <- function(data,
                                 chrom_col = "chromosome",
                                 start_col = "start",
                                 end_col = "end",
                                 sample_cols = c("sample1", "sample2", "sample3"), # multiple samples
                                 sample_colors = NULL,
                                 min_cn = 0,
                                 max_cn = NULL,
                                 line_size = 0.9) {
  
  # rename columns
  data <- data %>%
    rename(
      Chromosome = !!sym(chrom_col),
      Start = !!sym(start_col),
      End = !!sym(end_col)
    )
  
  # reshape to long format
  data_long <- data %>%
    pivot_longer(
      cols = all_of(sample_cols),
      names_to = "Sample",
      values_to = "CopyNumber"
    ) %>% 
    mutate(Sample = factor(Sample, levels = sample_cols))
  
  # sort chromosomes and calculate offsets
  chr_lengths <- data_long %>%
    group_by(Chromosome) %>%
    summarize(chr_length = max(End), .groups = "drop") %>%
    arrange(as.numeric(as.character(Chromosome))) %>%
    mutate(offset = lag(cumsum(as.numeric(chr_length)), default = 0))
  
  # add offsets
  data_offset <- data_long %>%
    left_join(chr_lengths, by = "Chromosome") %>%
    mutate(
      Start_genome = as.numeric(Start) + offset,
      End_genome = as.numeric(End) + offset
    )
  
  # define Y axis
  y_min <- if (!is.null(min_cn)) min_cn else floor(min(data_offset$CopyNumber, na.rm = TRUE)) - 0.5
  y_max <- if (!is.null(max_cn)) max_cn else ceiling(max(data_offset$CopyNumber, na.rm = TRUE)) + 0.5
  
  
  # plot
  p <- ggplot(data_offset, aes(x = Start_genome, xend = End_genome, y = CopyNumber, yend = CopyNumber, color = Sample)) +
    geom_segment(size = line_size) +
    geom_vline(xintercept = chr_lengths$offset, color = "gray80", linetype = "dashed")
  # Conditional color override
  if (!is.null(sample_colors)) {
    p <- p + scale_color_manual(values = sample_colors)
  }
  p <- p +
    scale_x_continuous(
      breaks = chr_lengths$offset + chr_lengths$chr_length / 2,
      labels = paste0(chr_lengths$Chromosome),
      expand = expansion(mult = c(0.01, 0.01))
    ) +
    scale_y_continuous(
      limits = c(y_min, y_max),
      breaks = seq(ceiling(y_min), floor(y_max), by = 0.5),  # force integer ticks
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    labs(
      title = "Whole Genome Copy Number Profile (Multiple Samples)",
      x = "Chromosomes",
      y = "Copy number value",
      color = 'Sample' # add legend or text here
    ) +
    theme_minimal(base_size = 14) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.grid.major.x = element_blank(),
      legend.title = element_text(size = 10)
    )
  
  return(p)
}


# cluster-LPAC
## INPUT
#  bin_val     : Nb x Ns size array containing relative CNA values of bins(not normalized)
#  of equal size. With Nb the number of bins, and Ns the number of multi-region samples.
#  seg_val     : Nb x Ns size array containing values of segment length.
#  Cval        : define the cut-off consine similarity value
#  Dmin        : define initial minimum absCNA distance(CNH)

## OUTPUT: 
#  LPAC_success: Binary variable with 1 if LPAC inference was succesfull and 0 if unsuccesfull
#  ploidies    : Ns x 1 size vector containing the tumor ploidies
#  purities    : Ns x 1 size vector containing the sample purities
#  CNH         : Ns x 1 size vector containing the copy number heterogeneity
#  best_sample : Ns x 1 size binary vector in which the best sample used for multi-region 
#  inference is identified
cluster_LPAC <- function(bin_val,seg_len=NULL, n_group=1, CNHout=NA, ploidies=NA, purities=NA, sds=NULL, Cval) {
  # Determine sample size and number of bins
  Nb <- nrow(bin_val)
  Ns <- ncol(bin_val)
  
  # define default segment length
  if (is.null(seg_len)){seg_len <- data.frame(matrix(1, nrow = Nb, ncol = Ns))}
  
  # Check if data consists of at least two samples
  if (Ns < 2) {
    return(list(LPAC_success=FALSE, ploidies=ploidies, purities=purities, CNH=CNHout, best_sample=NA, cluster=rep(n_group,Ns)))
  }
  
  # Check if data is normalized
  for (i in 1:Ns) {
    if (abs(mean(bin_val[, i]) - 1) > 0.02) {
      stop(paste("Input bin values (column", i, ") are not normalized. Each sample should have a mean of ~1."))
    }
  }
  
  # Initialize output variables
  cluster <- rep(n_group,Ns)
  success_list <- rep(FALSE,Ns)
  
  # calculate the CNH, ploidies and purities if they are absent
  if ( any(is.na(CNHout)) | any(is.na(ploidies)) | any(is.na(purities)) ){
    # Initialize output variables
    CNHout <- rep(0, Ns)
    ploidies <- rep(0, Ns)
    purities <- rep(0, Ns)
    sds <- rep(0,Ns) # standard deviation list
    # Single-sample inference of absolute CNAs
    for (i in 1:Ns) {
      result <- CNH(bin_val[,i], seg_len[,i][!is.na(seg_len[, i])], NULL, NULL)
      CNHout[i] <- result$CNH_out
      ploidies[i] <- result$ploidy_out
      purities[i] <- result$purity_out
      sds[i] <- sdr(bin_val[,i],seg_len[,i])
    }
  }


  # Step 2: Criterion for success
  step1_success <- !(ploidies < 2.5 & purities == 1)
  
  if (sum(step1_success) == 0) {
    LPAC_success <- FALSE
    best_sample <- rep(FALSE, Ns)
    return(list(LPAC_success=success_list, ploidies=ploidies, purities=purities, CNH=CNHout, best_sample=best_sample, cluster=cluster))
  } else {
    LPAC_success <- TRUE
  }
  
  # Step 3: Identify best single sample
  if (LPAC_success) {
    sd_cutoff <- 0.015
    filter_conditon <- step1_success & (sds > sd_cutoff) # select the success sample which SD is high
    id <- which(CNHout == min(CNHout[filter_conditon]))[1] # smallest CNH in success samples
    best_sample <- rep(FALSE, Ns)
    best_sample[id] <- TRUE
    
    success_list <- rep(TRUE, Ns) # update LPAC_success
    
    # Step 4: Infer CNAs by alignment
    a1 <- (purities[id] * ploidies[id] + (1 - purities[id]) * 2) / purities[id]
    a2 <- -2 * (1 - purities[id]) / purities[id]
    cna_abs_ref <- a1 * bin_val[, id] + a2
    
    # Define search grid
    alphas <- seq(0.1, 1, by=0.01)
    taus <- seq(1.5, 5, by=0.05)
    
    # Loop over non-best samples
    ids <- which(!best_sample)
    for (i in ids) {
      Dmin <- 1
      D_neverless <- TRUE # define if the distance is never less than Dmin
      bin_val_test <- bin_val[, i]
      seg_len_test <- seg_len[, i]
      
      # Compute cosine similarity of relative copy number
      Ctest <- cosine_similarity(bin_val_test,bin_val[,best_sample],seg_len_test)
      
      # same cluster judgement
      if ((Ctest >= Cval) & (sds[i] > sd_cutoff) ){ 
        for (alpha in alphas) {
          for (tau in taus) {
            a1 <- (alpha * tau + (1 - alpha) * 2) / alpha
            a2 <- -2 * (1 - alpha) / alpha
            bin_val_test_abs <- a1 * bin_val_test + a2
            
            # Compute distance
            Dtest <- sum((abs(bin_val_test_abs - cna_abs_ref))*seg_len_test)/sum(seg_len_test)
            
            # Update params if this sample meets all requirement
            if (Dtest < Dmin) {
              Dmin <- Dtest # update minimum distance
              purities[i] <- alpha
              ploidies[i] <- tau
              D_neverless <- FALSE
            }
          }
        }
        
        # if distance is less than Dmin, refine result and compute CNH
        if (!D_neverless) {
          refine_result <- refine_search(bin_val_test,seg_len_test,ploidies[i],purities[i])
          ploidies[i] <- refine_result$ploidy
          purities[i] <- refine_result$purity
          CNHout[i] <- refine_result$CNH
        }
      }
      
      # update cluster judgement
      if (((Ctest < Cval) | (D_neverless)) & (sds[i] > sd_cutoff)) {cluster[i] <- n_group+1} 
      
      # useless sample judgement
      if (sds[i] < sd_cutoff){
        success_list[i] <- FALSE
        cluster[i] <- NA
      }
    }
      
    # get index
    new_cluster_id <- which(cluster!=n_group)
    
    # new LPAC for new cluster
    if (length(new_cluster_id) != 0) {
      new_result <- cluster_LPAC(bin_val =  bin_val[,new_cluster_id,drop = FALSE], # make sure df format
                         seg_len = seg_len[,new_cluster_id,drop = FALSE],
                         n_group=n_group+1,
                         CNHout=CNHout[new_cluster_id], 
                         ploidies=ploidies[new_cluster_id], 
                         purities=purities[new_cluster_id],
                         sds=sds[new_cluster_id],
                         Cval = Cval)
      success_list[new_cluster_id] <- new_result$LPAC_success
      ploidies[new_cluster_id] <- new_result$ploidies
      purities[new_cluster_id] <- new_result$purities
      CNHout[new_cluster_id] <- new_result$CNH
      best_sample[new_cluster_id] <- new_result$best_sample
      cluster[new_cluster_id] <- new_result$cluster
    }
  }

  
  # Return results
  return(list(LPAC_success=success_list, 
              ploidies=ploidies, 
              purities=purities, 
              CNH=CNHout, 
              best_sample=best_sample, 
              cluster=cluster))
}


# CNH
## INPUT
#  seg_val: vector of relative copy number per segment
#  seg_len: vector of segment lengths
#  ploidy: vector of ploidy search range
#  purity: vector of purity search range
## OUTPUT
#  CNH_out: inferred CNH
#  ploidy_out: inferred ploidy
#  purity_out: inferred purity
CNH <- function(seg_val, seg_len, ploidy = NULL, purity = NULL) {
  # Check input validity
  if (!(is.vector(seg_val) && is.vector(seg_len) && length(seg_val) == length(seg_len))) {
    stop("Segment values and segment lengths must be column vectors of equal length.")
  }
  
  if (!(is.null(ploidy) || (is.numeric(ploidy) && length(ploidy) == 1 && ploidy > 0))) {
    stop("Ploidy must be a positive scalar or NULL.")
  }
  
  if (!(is.null(purity) || (is.numeric(purity) && length(purity) == 1 && purity > 0 && purity <= 1))) {
    stop("Purity must be a scalar between 0 and 1, or NULL.")
  }
  
  # Default range for grid search if ploidy is empty
  if (is.null(ploidy)) {
    ploidy <- seq(1.5, 5, by = 0.01)  # Tumor ploidy
  }
  
  # Default range for grid search if purity is empty
  if (is.null(purity)) {
    purity <- seq(0.1, 1, by = 0.01)  # Tumor purity
  }
  
  # Number of ploidies and purities in grid search
  Nploidy <- length(ploidy)
  Npurity <- length(purity)
  
  # Initialize vectors a1, a2 for transformations
  a1 <- numeric(Nploidy * Npurity)
  a2 <- numeric(Nploidy * Npurity)
  purity_all <- numeric(Nploidy * Npurity)
  ploidy_all <- numeric(Nploidy * Npurity)
  
  # Compute transformation coefficients
  for (i in 1:Nploidy) {
    start_idx <- (i - 1) * Npurity + 1
    end_idx <- i * Npurity
    a1[start_idx:end_idx] <- (purity * ploidy[i] + 2 * (1 - purity)) / purity
    a2[start_idx:end_idx] <- -2 * (1 - purity) / purity
    purity_all[start_idx:end_idx] <- purity
    ploidy_all[start_idx:end_idx] <- rep(ploidy[i], Npurity)
  }
  
  # Initialize output variables
  CNH_out <- 1
  purity_out <- 0
  ploidy_out <- 0
  
  # Grid search to infer CNH
  for (i in 1:(Nploidy * Npurity)) {
    # Transform relative CN values to absolute values
    q <- a1[i] * seg_val + a2[i]
    
    # Compute distance to nearest integer
    q_dist_down <- q %% 1
    q_dist_up <- 1 - q_dist_down
    q_dist_min <- pmin(q_dist_up, q_dist_down)
    
    # Compute weighted mean CNH
    CNHnew <- sum(q_dist_min * seg_len) / sum(seg_len)
    
    # Update CNH output if better than previous best
    if (CNHnew < CNH_out) {
      CNH_out <- CNHnew
      purity_out <- purity_all[i]
      ploidy_out <- ploidy_all[i]
    }
  }
  
  # refine result
  refine_result <- refine_search(seg_val, seg_len, ploidy_out, purity_out)
  
  # Return results as a named list
  return(list(CNH_out = refine_result$CNH, 
              ploidy_out = refine_result$ploidy, 
              purity_out = refine_result$purity))
}


# function to refine the result
refine_search <- function(seg_val, seg_len, ploidy_out, purity_out, refine_range=0.1, refine_step=0.001){
  # Define refined search ranges
  purity_min <- max(0.1, purity_out - refine_range)
  purity_max <- min(1.0, purity_out + refine_range)
  ploidy_min <- max(1.5, ploidy_out - refine_range)
  ploidy_max <- min(5.0, ploidy_out + refine_range)
  
  # Create refined grids
  purity_refined <- seq(purity_min, purity_max, by = refine_step)
  ploidy_refined <- seq(ploidy_min, ploidy_max, by = refine_step)
  
  # Number of refined grid points
  Nploidy_refined <- length(ploidy_refined)
  Npurity_refined <- length(purity_refined)
  
  # Initialize refined vectors
  a1_refined <- numeric(Nploidy_refined * Npurity_refined)
  a2_refined <- numeric(Nploidy_refined * Npurity_refined)
  purity_all_refined <- numeric(Nploidy_refined * Npurity_refined)
  ploidy_all_refined <- numeric(Nploidy_refined * Npurity_refined)
  
  # Initialize output variables
  CNH_out <- 1
  
  # Compute refined transformation coefficients
  for (i in 1:Nploidy_refined) {
    start_idx <- (i - 1) * Npurity_refined + 1
    end_idx <- i * Npurity_refined
    a1_refined[start_idx:end_idx] <- (purity_refined * ploidy_refined[i] + 2 * (1 - purity_refined)) / purity_refined
    a2_refined[start_idx:end_idx] <- -2 * (1 - purity_refined) / purity_refined
    purity_all_refined[start_idx:end_idx] <- purity_refined
    ploidy_all_refined[start_idx:end_idx] <- rep(ploidy_refined[i], Npurity_refined)
  }
  
  # Refined grid search
  for (i in 1:(Nploidy_refined * Npurity_refined)) {
    # Transform relative CN values to absolute values
    q <- a1_refined[i] * seg_val + a2_refined[i]
    
    # Compute distance to nearest integer
    q_dist_down <- q %% 1
    q_dist_up <- 1 - q_dist_down
    q_dist_min <- pmin(q_dist_up, q_dist_down)
    
    # Compute weighted mean CNH
    CNHnew <- sum(q_dist_min * seg_len) / sum(seg_len)
    
    # Update CNH output if better than previous best
    if (CNHnew < CNH_out) {
      CNH_out <- CNHnew
      purity_out <- purity_all_refined[i]
      ploidy_out <- ploidy_all_refined[i]
    }
  }
  
  # Return results
  return(list(
    CNH = CNH_out,
    ploidy = ploidy_out,
    purity = purity_out
  ))
}


# function to create a uniform segment length for the inconsistent segment length
create_sequence_df <- function(df, chrom, max_iterations = 500) {
  # 输入验证
  if(!all(c("start", "end") %in% names(df))) {
    stop("The data frame must contain ‘start’ and ‘end’ columns.")
  }
  
  # 初始化结果数据框
  result_df <- data.frame(chromosome = integer(),start = numeric(0), end = numeric(0))
  
  # 辅助函数：找到列中大于x且最接近x的值
  find_closest_greater <- function(x, column) {
    valid_values <- column[!is.na(column)]
    greater_values <- valid_values[valid_values > x]
    if(length(greater_values) == 0) return(NA)
    return(min(greater_values))
  }
  
  # 步骤1：找到start列中的最小值
  current_start <- min(df$start, na.rm = TRUE)
  
  # 开始循环
  for(i in 1:max_iterations) {
    # find the closest 'end' value
    closest_end <- find_closest_greater(current_start,df$end)
    
    # add to result
    result_df <- rbind(result_df, data.frame(chromosome = chrom, start = current_start, end = closest_end))
    
    # find closest 'start' value
    next_start <- find_closest_greater(closest_end, df$start)
    
    # 如果找不到新的start，结束循环
    if(is.na(next_start)) {
      break
    }
    
    current_start <- next_start
  }
  
  return(result_df)
}


# function to map segment values, df1 is template, df2 contains values
map_values_vectorized <- function(df1, df2, shrink_percent=0) {
  sample <- colnames(df2)[length(colnames(df2))]
  # 初始化结果列
  df1[[sample]] <- NA
  
  # 对每个df1的行找到匹配的df2行
  for (i in 1:nrow(df1)) {
    # 计算区间长度
    interval_length <- df1$end[i] - df1$start[i]
    shrink_amount <- interval_length * shrink_percent
    
    # 计算df1的调整后区间（向内缩小）
    adjusted_start <- df1$start[i] + shrink_amount
    adjusted_end <- df1$end[i] - shrink_amount
    
    # 找到所有包含当前df1区间的df2行
    matches <- which(df2$start <= adjusted_start & 
                     df2$end >= adjusted_end)
    
    # 如果有匹配，取第一个
    if (length(matches) > 0) {
      df1[i,sample] <- df2[matches[1],sample]
    }
  }
  
  # 删除临时列
  df1$adjusted_start <- NULL
  df1$adjusted_end <- NULL
  
  return(df1)
}


# function to create dataframe of segment value and length for the inconsistent input data (e.g. each sample has unique segment length)
## INPUT
# splited_data    : a sparse dataframe that contains chromosome, start, end, feature/length and sample column of one patient
create_segment_val <- function(splited_data){
  result_df <- NULL
  # check sample number
  if(ncol(splited_data)==5){
    result_df <- na.omit(splited_data)
    seg_len <- result_df[,5,drop=FALSE]
    seg_len[,1] <- result_df[,4]
    return(list(seg_val = result_df[,5,drop=FALSE],seg_len = seg_len, seg_info = result_df))
  }
  for(chrom in c(1:24)){
    # select chromosome
    chrom_data <- splited_data[splited_data$chromosome==chrom,]
    # remove all na rows
    chrom_data <- chrom_data[rowSums(!is.na(chrom_data[, 5:ncol(chrom_data)])) > 0, ]
    # create a template segment
    temp_df <- create_sequence_df(chrom_data,chrom)
    # copy template column name to result df
    if (is.null(result_df)){result_df <- temp_df[0,]}
    
    data_val <- chrom_data[,5:ncol(chrom_data)]
    samples <- colnames(data_val)
    # map each sample to template
    for (sample in samples){
      # remove na rows
      sample_val <- (chrom_data[!is.na(data_val[,sample]),c("chromosome","start","end",sample)])
      # remap segment
      temp_df <- map_values_vectorized(temp_df,sample_val,shrink_percent = 0.1)
    }
    result_df <- rbind(result_df,temp_df)
  }
  # remove NA rows
  result_df <- na.omit(result_df)
  # extract segment value and length
  seg_val <- result_df[,4:ncol(result_df),drop=FALSE]
  seg_len <- seg_val
  len <- result_df$end - result_df$start
  for (col in colnames(seg_val)){
    seg_len[[col]] <- len
  }
  # add length col
  result_df <- cbind(result_df[1:(4 - 1)], length = len, result_df[4:ncol(result_df)])
  return(list(seg_val = seg_val,seg_len = seg_len, seg_info = result_df))
}
```

```{r}
library(dplyr)
# function to apply LPAC to dataset that contains multi-patients
## INPUT
# bin_val   : a Nb x Ns dataframe that contains relative copy number of each sample, if segment is TRUE, then bin_val should be raw data
# seg_len   : a Nb x Ns dataframe that contains segments length of each sample
# data_id   : a 2 x Ns dataframe that contains samplenames and corresponding patient_id
# segment   : TRUE means use the create_segment_val function to handle inconsistent input data
# path      : output file path
## OUTPUT
# result_df : a dataframe that contains result of LPAC
run_multi_patient <- function(bin_val,data_id,path,seg_len=NULL,segment=FALSE, Cval){
  # Create an empty data frame to store results
  results_df <- data.frame(patient_id = integer(),
                           sample_name = character(),
                           LPAC_success = numeric(),
                           ploidies = numeric(),
                           purities = numeric(),
                           CNH = numeric(),
                           best_sample = character(),
                           cluster = numeric(),
                           stringsAsFactors = FALSE)
  
  # loop all patients
  id_list <- sort(unique(data_id$patient_id))
  
  for (id in id_list){
    print(paste0('Proccessing ',id,' ...'))
    # get sample names that belong to specific id
    selected_samples <- data_id %>%
    filter(patient_id == id) %>%
    pull(samplenames)

    if(segment){
      # Subset data(bin_val) to selected samples and front columns
      splited_data <- cbind(bin_val[,1:4],bin_val[, selected_samples])
      # remap the segments
      result <- create_segment_val(splited_data)
      seg_val <- result$seg_val
      splited_len <- result$seg_len
      # Normalize each column (sample) to have mean ~1
      splited_data <- sweep(seg_val, 2, colMeans(seg_val), "/")
      }
    else{
      if(is.null(seg_len)){splited_len <- NULL}
      else{splited_len<-seg_len[,selected_samples]}
      # Subset bin_val to select columns corresponding to the selected sample names
      splited_data <- bin_val[,selected_samples,drop=FALSE]
      }
  
    # Use tryCatch to handle errors
    result <- tryCatch({
      # Apply the LPAC function
      cluster_LPAC(bin_val = splited_data, seg_len = splited_len, Cval = Cval)
    }, error = function(e) {
      # In case of error, return NULL or some default value
      message(paste("Error with patient_id:", id, "-", conditionMessage(e)))  # Print message without stopping
      return(NULL)
    })
  
    # If result is not NULL (no error), store it
    if (!is.null(result)) {
      # Store the results in the data frame
      results_df <- rbind(results_df, data.frame(
        patient_id = rep(id,length(selected_samples)),
        sample_name = selected_samples,
        LPAC_success = result$LPAC_success,
        ploidies = result$ploidies,
        purities = result$purities,
        CNH = result$CNH,
        best_sample = result$best_sample,
        cluster = result$cluster
      ))
    } else {
      # Handle cases where LPAC fails for a specific patient (optional)
      results_df <- rbind(results_df, data.frame(
        patient_id = id,
        sample_name = selected_samples,
        LPAC_success = FALSE,
        ploidies = NA,
        purities = NA,
        CNH = NA,
        best_sample = NA,
        cluster = NA
      ))
    }
  }
  
  # Output result
  write.csv(results_df, path, row.names = FALSE)
}
path <- "LPAC_v6_EAC_93.csv"
run_multi_patient(bin_val = bin_val, data_id = data_id, segment = FALSE, path = path, Cval = 0.93)
```

# Run
```{r}
# Apply LPAC to one patient
library(dplyr)

# Get sample names that belong to specific id
selected_samples <- data_id %>%
  filter(patient_id == 'C519') %>%
  pull(samplenames)

# # Subset bin_val to select columns corresponding to the selected sample names
# splited_data <- cbind(data[,1:4],data[, selected_samples])
# result <- create_segment_val(splited_data)
# seg_val <- result$seg_val
# seg_len <- result$seg_len
# # Normalize each column (sample) to have mean ~1
# seg_val <- sweep(seg_val, 2, colMeans(seg_val), "/")

# Apply method
result <- cluster_lpac_main(bin_val =  bin_val[,selected_samples,drop=FALSE],plot = TRUE, seg_info = data[,0:4],img_path = "test")
mat <- do.call(cbind,result)

```

```{r}
# Transform tracex data into remap data
patient_id <- sort(unique(data_ref$patient_id))
for (id in patient_id){
  # Get sample names that belong to specific id
  selected_samples <- data_id %>%
    filter(patient_id == id) %>%
    pull(samplenames)
  
  # Subset bin_val to select columns corresponding to the selected sample names
  splited_data <- cbind(data[,1:4],data[, selected_samples])
  # remap
  result <- create_segment_val(splited_data)
  output_df <- result$seg_info
  # save
  path <- paste0("nscls_segments/",id,".csv")
  write.csv(output_df, path, row.names = FALSE)
}
```


# Test area
```{r}
test_list <- result$test
test_D <- unlist(test_list[1])
test_id <- unlist(test_list[2])

test_D
test_id
```

```{r}

data_1 <- read.table(file.choose(), header = TRUE, sep = "\t")
```

```{r}
# test the CNH output

one_sample_data <- bin_val[,'CF11675']

test_result_1 <- CNH(one_sample_data, rep(1, length(one_sample_data)), NULL, NULL)
t_seq <- seq(1.5, 5, by = 0.01)
a_seq <- seq(0.1, 1, by = 0.01)
test_result_2 <- compute_best_params_for_sample(one_sample_data,rep(1, length(one_sample_data)),a_seq,t_seq)
```

```{r}
# purities <- data_1$purities_multi[data_1$samplenames %in% selected_samples]
# ploidies <- data_1$ploidies_multi[data_1$samplenames %in% selected_samples]
purities <- result$purities
ploidies <- result$ploidies

# check if sum and guess ploidy is identical
check_ploidies <- function(purities,ploidies,splited_data) {
  a1 <- (purities * ploidies + (1 - purities) * 2) / purities
  a2 <- -2 * (1 - purities) / purities
  cna <- a1 * splited_data + a2
  print(colMeans(cna))
  print(ploidies)
}

check_ploidies(purities,ploidies,splited_data)
```
