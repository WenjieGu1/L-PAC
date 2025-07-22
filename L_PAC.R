```{r}
# Load the dataset
data <- read.table(file.choose(), header = TRUE, sep = "\t", comment.char = "#")

data_ref <- read.table(file.choose(), header = TRUE, sep = "\t")
# get ids and sample names
data_id <- data_ref[,1:2]

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
bin_val <- sweep(bin_val, 2, colMeans(bin_val), "/")

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

```{r}
# LPAC
## INPUT
#  bin_val     : Nb x Ns size array containing relative CNA values of bins
#  of equal size. With Nb the number of bins, and Ns the number of multi-region samples.
## OUTPUT: 
#  LPAC_success: Binary variable with 1 if LPAC inference was succesfull and 0 if unsuccesfull
#  ploidies    : Ns x 1 size vector containing the tumor ploidies
#  purities    : Ns x 1 size vector containing the sample purities
#  CNH         : Ns x 1 size vector containing the copy number heterogeneity
#  best_sample : Ns x 1 size binary vector in which the best sample used for multi-region 
#  inference is identified
# **added a Test to Dtest
# **changed Dmin value, changed CNH input
LPAC <- function(bin_val,seg_len=NULL) {
  # Determine sample size and number of bins
  Nb <- nrow(bin_val)
  Ns <- ncol(bin_val)
  
  # define default segment length
  if (is.null(seg_len)){seg_len <- data.frame(matrix(1, nrow = Nb, ncol = Ns))}
  
  # Check if data consists of at least two samples
  if (Ns < 2) {
    stop("Input should consist of multiple samples (columns).")
  }
  
  # Check if data is normalized
  for (i in 1:Ns) {
    if (abs(mean(bin_val[, i]) - 1) > 0.02) {
      stop(paste("Input bin values (column", i, ") are not normalized. Each sample should have a mean of ~1."))
    }
  }
  
  # Initialize output variables
  CNHout <- rep(0, Ns)
  ploidies <- rep(0, Ns)
  purities <- rep(0, Ns)
  
  # Single-sample inference of absolute CNAs
  for (i in 1:Ns) {
    result <- CNH(bin_val[, i], seg_len[,i], NULL, NULL) # CNH generalized
    CNHout[i] <- result$CNH_out
    ploidies[i] <- result$ploidy_out
    purities[i] <- result$purity_out
  }

  # Step 2: Criterion for success
  step1_success <- !(ploidies < 2.5 & purities == 1)
  
  if (sum(step1_success) == 0) {
    LPAC_success <- FALSE
    best_sample <- rep(FALSE, Ns)
    return(list(LPAC_success=LPAC_success, ploidies=ploidies, purities=purities, CNH=CNHout, best_sample=best_sample))
  } else {
    LPAC_success <- TRUE
  }
  
  # Step 3: Identify best single sample
  if (LPAC_success) {
    id <- which(CNHout == min(CNHout[step1_success]))[1] # smallest CNH in success samples
    best_sample <- rep(FALSE, Ns)
    best_sample[id] <- TRUE

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
      bin_val_test <- bin_val[, i]
      seg_len_test <- seg_len[, i]
      
      
      for (alpha in alphas) {
        for (tau in taus) {
          a1 <- (alpha * tau + (1 - alpha) * 2) / alpha
          a2 <- -2 * (1 - alpha) / alpha
          bin_val_test_abs <- a1 * bin_val_test + a2
          
          # Compute distance
          Dtest <- sum((abs(bin_val_test_abs - cna_abs_ref))*seg_len_test)/sum(seg_len_test)
          
          # Update best match
          if (Dtest < Dmin) {
            Dmin <- Dtest
            purities[i] <- alpha
            ploidies[i] <- tau
          }
        }
      }
      
      # Compute CNH
      a1 <- (purities[i] * ploidies[i] + (1 - purities[i]) * 2) / purities[i]
      a2 <- -2 * (1 - purities[i]) / purities[i]
      cna_abs_test <- a1 * bin_val_test + a2
      CNHout[i] <- sum((pmin(abs(cna_abs_test %% 1), abs(1 - cna_abs_test %% 1))) * seg_len_test) / sum(seg_len_test)
    }
  }
  
  # Return results
  return(list(LPAC_success=LPAC_success, ploidies=ploidies, purities=purities, CNH=CNHout, best_sample=best_sample))
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
  
  # Return results as a named list
  return(list(CNH_out = CNH_out, ploidy_out = ploidy_out, purity_out = purity_out))
}

# -----below is substitute function for CNH, you can ignore it-----
# Function to compute CNH for a single sample
compute_CNH <- function(purity, ploidy, n_values, seg_len) {
  a1 <- (purity * ploidy + 2 * (1 - purity)) / purity
  a2 <- -2 * (1 - purity) / purity
  q <- a1 * n_values + a2  # Calculate q

  # Calculate the distance to the nearest integer
  q_dist_down <- q %% 1
  q_dist_up <- 1 - q_dist_down
  q_dist_min <- pmin(q_dist_up, q_dist_down)

  # Calculate CNH
  CNH <- sum(q_dist_min * seg_len) / sum(seg_len)
  return(CNH)
}

# Function to compute best CNH for a sample and its params
compute_best_params_for_sample <- function(n_values, seg_len, a_seq=NULL, t_seq=NULL) {
  best_CNH <- Inf
  best_purity <- NA
  best_ploidy <- NA
  
    # Default range for grid search if ploidy is empty
  if (is.null(t_seq)) {
    t_seq <- seq(1.5, 5, by = 0.01)  # Tumor ploidy
  }
  
  # Default range for grid search if purity is empty
  if (is.null(a_seq)) {
    a_seq <- seq(0.1, 1, by = 0.01)  # Tumor purity
  }
  
  # grid search
  for (purity in a_seq) {
    for (ploidy in t_seq) {
      CNHnew <- compute_CNH(purity, ploidy, n_values, seg_len)
      
      # update params
      if (CNHnew < best_CNH) {
        best_CNH <- CNHnew
        best_purity <- purity
        best_ploidy <- ploidy
      }
    }
  }

  return(list(purity_out=best_purity, ploidy_out=best_ploidy, CNH_out=best_CNH))
}


# function to create the segment information
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

# function to create dataframe of segment value and length
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
  seg_val <- result_df[,4:ncol(result_df)]
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
segment <- TRUE
bin_val <- data
library(dplyr)
# Create an empty data frame to store results
results_df <- data.frame(patient_id = integer(),
                         sample_name = character(),
                         LPAC_success = numeric(),
                         ploidies = numeric(),
                         purities = numeric(),
                         CNH = numeric(),
                         best_sample = character(),
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
    splited_data <- bin_val[,selected_samples]
    }

  # Use tryCatch to handle errors
  result <- tryCatch({
    # Apply the LPAC function
    LPAC(splited_data)
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
      LPAC_success = rep(result$LPAC_success,length(selected_samples)),
      ploidies = result$ploidies,
      purities = result$purities,
      CNH = result$CNH,
      best_sample = result$best_sample
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
      best_sample = NA
    ))
  }
}

# Output result
write.csv(results_df, "LPAC_tracex_result.csv", row.names = FALSE)
```

```{r}
# Apply LPAC to one patient
library(dplyr)
# Get sample names that belong to specific id
selected_samples <- data_id %>%
  filter(patient_id == 15) %>%
  pull(samplenames)

# Subset bin_val to select columns corresponding to the selected sample names
splited_data <- bin_val[, selected_samples]

# Apply method
result <- LPAC(splited_data)
mat <- do.call(cbind,result)

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
