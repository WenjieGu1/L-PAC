# compare reference and result
```{r}
## Load data

# reference is standard LPAC,  result is new LPAC
reference_table <- read.csv(file.choose())
# reference_table <- read.table(file.choose(), header = TRUE, sep = "\t")
result_table <- read.csv(file.choose())

# names(reference_table)[names(reference_table) == "samplenames"] <- "sample_name" # change colname if ref and result colname are different
# names(reference_table)[names(reference_table) == "CNH"] <- "CNH_1"
# names(reference_table)[names(reference_table) == "CNH_multi"] <- "CNH" 
# names(reference_table)[names(reference_table) == "ploidies_multi"] <- "ploidies"
# names(reference_table)[names(reference_table) == "purities_multi"] <- "purities" 
```

```{r}
# array bar plot
df <- result_table[result_table$LPAC_success==TRUE,]

library(ggplot2)
ggplot(df, aes(x = sample_name, y = ploidies, fill = patient_id)) +
  geom_bar(stat = "identity") +
  facet_wrap(~ patient_id, scales = "free_x", ncol = 5) +
  labs(title = "Ploidies per Sample for Each Patient",
       x = "Sample",
       y = "Ploidy") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text.x = element_text(angle = 45, hjust = 1))
  
ggsave("patient_array_barplot.png", width = 12, height = 8, dpi = 300)

```


```{r}
library(tidyverse)
## Plot metric for each patient

# Add source labels and reshape to long format
ref_long <- reference_table %>%
  mutate(source = "reference") %>%
  pivot_longer(cols = c(ploidies, purities, CNH), names_to = "metric", values_to = "value") # set the metrics

res_long <- result_table %>%
  mutate(source = "cluster_LPAC") %>%
  pivot_longer(cols = c(ploidies, purities, CNH), names_to = "metric", values_to = "value") # set the metrics

# Combine both
df <- bind_rows(ref_long, res_long)

# Get unique patient IDs
all_ids <- unique(na.omit(df$patient_id))

# Generate and save plots
for (id in all_ids) {
  df_sub <- df %>% filter(patient_id == id)

  p <- ggplot(df_sub, aes(x = sample_name, y = value, fill = source)) +
    geom_bar(stat = "identity", position = "dodge") +
    facet_wrap(~ metric, scales = "free_y") +
    labs(
      title = paste("Patient", id, "- Comparison of Metrics"),
      x = "Sample", y = "Value"
    ) +
    theme_minimal(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))

  ggsave(paste0("patient_", id, "_comparison_plot.png"), plot = p, width = 10, height = 4)
}
```

```{r}
## Plot all patients in one image (compare ref and new result)
library(dplyr)
library(ggplot2)
library(RColorBrewer)

# Add source label
dfA <- reference_table %>% filter(LPAC_success == TRUE) %>% select(patient_id, sample_name, purities, ploidies,CNH) # set the metrics
dfB <- result_table %>% filter(LPAC_success == TRUE) %>% select(patient_id, sample_name, purities, ploidies,CNH,cluster) # set the metrics

# Select cluster > 1 sample
selected_sample <- dfB[dfB$cluster > 1, "patient_id"]
dfA <- dfA %>% filter(patient_id %in% selected_sample)
dfB <- dfB %>% filter(patient_id %in% selected_sample)

# Add source labels
dfA$source <- "LPAC"
dfB$source <- "cluster_LPAC"

# Preserve the order of sample_name from dfA
sample_order <- dfA$sample_name

# Combine the datasets
df_combined <- bind_rows(dfA, dfB) %>%
  mutate(
    sample_name = factor(sample_name, levels = sample_order),  # Set factor levels using dfA order
    patient_id = as.factor(patient_id)  # patient_id as factor for color
  )


# Set colors for 24 IDs
id_colors <- brewer.pal(n = 12, name = "Paired")
id_colors <- rep(id_colors, length.out = length(levels(df_combined$patient_id)))

# Plotting function
plot_metric_bar <- function(data, metric_name) {
  data <- data %>% filter(!is.na(.data[[metric_name]]))  # Remove NAs, patient id = 0

  ggplot(data, aes(x = sample_name, y = .data[[metric_name]], fill = patient_id)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.7) +
    facet_wrap(~source) +
    scale_fill_manual(values = id_colors) +
    labs(title = paste("Comparison of", metric_name), y = metric_name, x = NULL) +
    theme_minimal() +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      legend.title = element_blank()
    )
}
# Save the plots with wider horizontal dimensions
ggsave("EAC_v4_90_002_purity2.png", plot = plot_metric_bar(df_combined, "purities"), width = 16, height = 6, dpi = 300)
ggsave("EAC_v4_90_002_ploidy2.png", plot = plot_metric_bar(df_combined, "ploidies"), width = 16, height = 6, dpi = 300)
ggsave("EAC_v4_90_002_CNH2.png", plot = plot_metric_bar(df_combined, "CNH"), width = 16, height = 6, dpi = 300)

```




# compare provided and findings
```{r}
# Load data
# reference is provided data, result is standard LPAC result
# reference_table <- read.table(file.choose(), header = TRUE, sep = "\t")
reference_table <- read.csv(file.choose())
result_table <- read.csv(file.choose())

```

```{r}
ref_col <- 'purities' # name in reference df
result_col <- 'purities' # name in result df
metric <- ''

# # Sort by patient_id, then samplename; only for EAC ref file
# reference_table <- reference_table[order(reference_table$patient_id, reference_table$samplenames), ]

# Delete NA rows of reference_table, and delete the same rows for result_table
na_samples <- reference_table[is.na(reference_table[[ref_col]]),"sample_name"]
reference_table_nona <- reference_table[!(reference_table$sample_name %in% na_samples),]
result_table_nona <- result_table[!(result_table$sample_name %in% na_samples),]
```

```{r}
## scatter plot for multi patient comparison
library(ggplot2)
library(RColorBrewer)

# combine data
df_combined <- data.frame(
  col1 = reference_table_nona[[ref_col]],
  col2 = result_table_nona[[result_col]],
  group = as.factor(reference_table_nona$patient_id)
)

# Generate enough colors
num_groups <- length(unique(df_combined$group))
palette_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_groups)

# Plot
ggplot(df_combined, aes(x = col1, y = col2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = palette_colors) +
  theme_minimal() +
  labs(x = "Purity", y = "new Purity", color = "Patient")
```

```{r}
# scatter plot but for one patient comparison

# Load data
# reference is provided data, result is standard LPAC result
# reference_table <- read.table(file.choose(), header = TRUE, sep = "\t")
df1 <- read.csv(file.choose())
df2 <- read.csv(file.choose())
```
```{r}
library(ggplot2)
library(RColorBrewer)

ref_col <- 'purities' # name in reference df
result_col <- 'purities' # name in result df

df_combined <- data.frame(
  col1 = df1[[ref_col]],
  col2 = df2[[result_col]],
  group = as.factor(df1$sample_name),
  patient_id = as.factor(df1$patient_id)
)

df_combined <- df_combined[df_combined$patient_id == 1,]

# Generate enough colors
num_groups <- length(unique(df_combined$group))
palette_colors <- colorRampPalette(brewer.pal(10, "Set1"))(num_groups)
# Plot
ggplot(df_combined, aes(x = col1, y = col2, color = group)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = palette_colors) +
  theme_minimal() +
  labs(x = "LPAC purities", y = "cluster-LPAC purities", color = "patient 1")

```


```{r}
## bar plot
library(ggplot2)

# Calculate the difference between ploidies
difference <- reference_table_nona[[ref_col]] - result_table_nona[[result_col]]

# Combine the data into a new data frame
df_combined <- data.frame(
  difference = difference,
  group = as.factor(reference_table_nona$patient_id),
  samplenames = as.factor(reference_table_nona$samplenames)
)

# Ensure samplenames are ordered as in df_combined
df_combined$samplenames <- factor(df_combined$samplenames, levels = df_combined$samplenames)

# Generate enough colors
num_groups <- length(unique(df_combined$group))
palette_colors <- colorRampPalette(brewer.pal(12, "Set3"))(num_groups)

# Create bar plot without x-axis labels
ggplot(df_combined, aes(x = samplenames, y = difference, fill = group)) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = palette_colors) +  # Use the custom colors
  theme_minimal() +
  labs(x = "Samples", y = "Difference (Purity from Provide - Purity from Find)", fill = "Patient Group") +
  theme(axis.text.x = element_blank())  # Remove x-axis labels completely


```

```{r}
# Filter high difference samples
diff_df <- cbind(result_table_nona,difference,reference_table_nona[[ref_col]])
colnames(diff_df)[ncol(diff_df)] <- ref_col
summary(reference_table_nona$ploidies_multi)
summary(diff_df$difference)
df_h <- diff_df[abs(diff_df$difference)>0.5,c("patient_id","sample_name",ref_col,result_col)]
unique(df_h$patient_id)
```


```{r}
# a stat boxplot to find the diverse ploidy
# Boxplot
library(hrbrthemes)
library(viridis)

# select success samples
result_table_success <- result_table[result_table$LPAC_success == 'TRUE',]

result_plot_data <- result_table_success %>% select(c('patient_id','ploidies'))
result_plot_data$patient_id <- as.factor(result_plot_data$patient_id)

# Plot
result_plot_data %>%
  ggplot( aes(x=patient_id, y=ploidies, fill=patient_id)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("Finding's Success Ploidies boxplot") +
    xlab("")

```



# evalutate result(mse)
```{r}
# Load data
# df_data <- read.table(file.choose(), header = TRUE, sep = "\t", comment.char = "#") # only copy number transform need it

# df_ref <- read.csv(file.choose(), header = TRUE) # for EPICC or NSCLS csv file
df_ref <- read.table(file.choose(), header = TRUE, sep = "\t") # only for EAC table

# Sort by patient_id, then samplename; only for EAC and NSCLS ref file
df_ref <- df_ref[order(df_ref$patient_id, df_ref$samplenames), ]

df_result <- read.csv(file.choose(), header = TRUE)
```

```{r}
# inference success rate and purity mse
mse <- function(actual, predicted) {
  mean((actual - predicted)^2,na.rm = TRUE)
}
library(dplyr)

# EAC
success_rate <- nrow(df_result[df_result$LPAC_success ==TRUE,])/nrow(df_result)
samples_im <- df_ref[!is.na(df_ref$purities_im),"samplenames"]

samples_evaluate <- df_result %>%
  filter(LPAC_success == TRUE, sample_name %in% samples_im) %>%
  pull(sample_name)

purity_mse <- mse(df_ref[df_ref$samplenames %in% samples_evaluate,"purities_im"],df_result[df_result$sample_name %in% samples_evaluate,"purities"])

# # EPICC or NSCLS
# success_rate <- nrow(df_result[df_result$LPAC_success ==TRUE,])/nrow(df_result)
# samples_im <- df_ref[!is.na(df_ref$purities),"samplenames"] # sample that have reference
# 
# samples_evaluate <- df_result %>%
#   filter(LPAC_success == TRUE, sample_name %in% samples_im) %>%
#   pull(sample_name)
# 
# purity_mse <- mse(df_ref[df_ref$samplenames %in% samples_evaluate,"purities"],df_result[df_result$sample_name %in% samples_evaluate,"purities"])
```

```{r}
# mse of cluster > 1
df_lpac <- read.csv(file.choose(), header = TRUE) # load standard LPAC result

# select successful and cluster > 1 and has reference purity sample
cluster_over1 <- df_result %>%
  filter(LPAC_success == TRUE, cluster > 1, sample_name %in% samples_im) %>%
  pull(sample_name)


out_mse <- mse(df_ref[df_ref$samplenames %in% cluster_over1,"purities_im"],df_result[df_result$sample_name %in% cluster_over1,"purities"]) # result of cluster lpac
out_mse2 <- mse(df_ref[df_ref$samplenames %in% cluster_over1,"purities_im"],df_lpac[df_lpac$sample_name %in% cluster_over1,"purities"]) # result of standard lpac
```

```{r}
# mse of ploidy
# NSCLS
# inference success rate and purity mse
mse <- function(actual, predicted) {
  mean((actual - predicted)^2,na.rm = TRUE)
}
library(dplyr)
samples_im <- df_ref[!is.na(df_ref$ploidies),"samplenames"]

samples_evaluate <- df_result %>%
  filter(LPAC_success == TRUE, sample_name %in% samples_im) %>%
  pull(sample_name)

ploidy_mse <- mse(df_ref[df_ref$samplenames %in% samples_evaluate,"ploidies"],df_result[df_result$sample_name %in% samples_evaluate,"ploidies"])
```
```{r}
# mse of cluster > 1
df_lpac <- read.csv(file.choose(), header = TRUE) # load standard LPAC result

# select successful and cluster > 1 and has reference purity sample
cluster_over1 <- df_result %>%
  filter(LPAC_success == TRUE, cluster > 1, sample_name %in% samples_im) %>%
  pull(sample_name)


out_mse <- mse(df_ref[df_ref$samplenames %in% cluster_over1,"ploidies"],df_result[df_result$sample_name %in% cluster_over1,"ploidies"]) # result of cluster lpac
out_mse2 <- mse(df_ref[df_ref$samplenames %in% cluster_over1,"ploidies"],df_lpac[df_lpac$sample_name %in% cluster_over1,"ploidies"]) # result of standard lpac
```



```{r}
# define a function to transform log2 relative copy number to absolute copy number
convert_rcn <- function(value, purity, ploidy){
  value <- 2^value
  a1 <- (purity*ploidy+2*(1-purity))/purity
  a2 <- 2*(1-purity)/purity
  q <- a1*value - a2
  return(q)
}
# define Mean Square Error function
mse <- function(actual, predicted) {
  mean((actual - predicted)^2)
}

# set ploidy and purity as NA if LPAC fail
df_result$ploidies[!df_result$LPAC_success] <- NA
df_result$purities[!df_result$LPAC_success] <- NA

samplenames <- df_result$sample_name
for (i in samplenames){
  ploidy_ref <- df_ref[df_ref$samplenames==i,"ploidies"]
  purity_ref <- df_ref[df_ref$samplenames==i,"purities"]
  ploidy_pred <- df_result[df_result$sample_name==i,"ploidies"]
  purity_pred <- df_result[df_result$sample_name==i,"purities"]
  q_ref <- convert_rcn(df_data[,i],purity_ref,ploidy_ref)
  q_pred <- convert_rcn(df_data[,i],purity_pred,ploidy_pred)
  q_mse <- mse(q_ref,q_pred)
  df_result$q_MSE[which(df_result$sample_name == i)] <- q_mse
}

# Output result
write.csv(df_result, "LPAC_cluster_EPICC_90_02_result.csv", row.names = FALSE)
```
```{r}
df_result <- read.csv(file.choose(), header = TRUE)
```
```{r}
mean(df_result$q_MSE, na.rm = TRUE)

# ratio of LPAC success
sum(df_result$LPAC_success==TRUE)/nrow(df_result)

# groups' mse boxplot
library(hrbrthemes)
library(viridis)
library(dplyr)
library(ggplot2)
# select success samples
result_table_success <- df_result[df_result$LPAC_success == 'TRUE',]

result_plot_data <- result_table_success %>% select(c('patient_id','q_MSE'))
result_plot_data$patient_id <- as.factor(result_plot_data$patient_id)

# Plot
boxplot <- result_plot_data %>%
  ggplot( aes(x=patient_id, y=q_MSE, fill=patient_id)) +
    geom_boxplot() +
    scale_fill_viridis(discrete = TRUE, alpha=0.6) +
    geom_jitter(color="black", size=0.4, alpha=0.9) +
    theme_ipsum() +
    theme(
      legend.position="none",
      plot.title = element_text(size=11)
    ) +
    ggtitle("LPAC success MSE boxplot") +
    xlab("")

ggsave("EPICC_90_02_mse_boxplot.png", plot = boxplot, width = 14, height = 5, dpi = 300)
```

```{r}
library(ggplot2)
library(dplyr)
library(tidyr)
library(patchwork)

# 数据
df <- data.frame(
  data_method = c("L-PAC on EAC", "cluster-LPAC on EAC",
                  "L-PAC on EPICC", "cluster-LPAC on EPICC",
                  "L-PAC on NSCLC", "cluster-LPAC on NSCLS"),
  dataset = c("EAC", "EAC", "EPICC", "EPICC", "NSCLC", "NSCLC"),
  method = c("L-PAC", "cluster-LPAC", "L-PAC", "cluster-LPAC", "L-PAC", "cluster-LPAC"),
  successful_inference_rate = c(0.85875, 0.77401, 0.75971, 0.70671, 0.82838, 0.73927),
  pu_mse_successful = c(0.09267, 0.07123, 0.04434, 0.02653, 0.08625, 0.04029),
  pu_mse_excluded = c(0.40684, 0.13761, 0.02966, 0.02609, 0.07018, 0.05406),
  pl_mse_successful = c(NA, NA, 0.58769, 0.48735, 1.53469, 1.55203),
  pl_mse_excluded = c(NA, NA, 0.6735, 0.49550, 1.48014, 1.52152)
)

# 转成长格式
df_long <- df %>%
  pivot_longer(
    cols = c(successful_inference_rate, pu_mse_successful, pu_mse_excluded,
             pl_mse_successful, pl_mse_excluded),
    names_to = "metric",
    values_to = "value"
  )

# Thesis-friendly theme
theme_thesis <- function() {
  theme_minimal(base_family = "serif", base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16, hjust = 0.5),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12),
      axis.text.x = element_text(angle = 30, hjust = 1, vjust = 1),
      legend.title = element_text(face = "bold", size = 13),
      legend.text = element_text(size = 12),
      panel.grid.major = element_line(color = "grey85"),
      panel.grid.minor = element_blank(),
      plot.margin = margin(10, 10, 10, 10)
    )
}

# 画图函数
plot_metric <- function(metric_name, custom_title) {
  df_plot <- df_long %>% filter(metric == metric_name)

  ggplot(df_plot, aes(x = dataset, y = value, fill = method)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.7), width = 0.6) +
    labs(title = custom_title, x = "Dataset", y = "Value") +
    scale_fill_manual(values = c("#1B9E77", "#D95F02")) + # Muted thesis colors
    theme_thesis()
}

# 绘图
plot1 <- plot_metric("successful_inference_rate", "Successful Inference Rate")
plot2 <- plot_metric("pu_mse_successful", "Purity MSE (Successful Inference)")
plot3 <- plot_metric("pu_mse_excluded", "Purity MSE (Non-main Samples)")
plot4 <- plot_metric("pl_mse_successful", "Ploidy MSE (Successful Inference)")
plot5 <- plot_metric("pl_mse_excluded", "Ploidy MSE (Non-main Samples)")

# 组合
combined_plot <- (plot1 + plot2 + plot3 + plot4 + plot5) +
  plot_layout(ncol = 3, guides = "collect") &
  theme(legend.position = "bottom")

# 保存
ggsave("metric_barplots3.png", plot = combined_plot,
       width = 12, height = 6, dpi = 300)



```

```{r}
# 输入数据
data <- data.frame(
  CSR_cutoff = c(0.92, 0.93, 0.94, 0.95, 0.96),
  successful_inference_rate = c(0.77401, 0.54237, 0.53107, 0.51412, 0.49717),
  MSE_purity_inferred = c(0.07123, 0.06046, 0.06081, 0.06816, 0.07151),
  MSE_purity_excluded = c(0.13761, 0.05697, 0.05697, 0.07033, 0.07856)
)

library(ggplot2)

# 比例系数
scaleFactor <- max(data$successful_inference_rate) / max(c(data$MSE_purity_inferred, data$MSE_purity_excluded))

line_plot <- ggplot(data, aes(x = CSR_cutoff)) +
  # 左轴
  geom_line(aes(y = successful_inference_rate, color = "Successful inference rate", linetype = "solid"), size = 1.5) +
  geom_point(aes(y = successful_inference_rate, color = "Successful inference rate", shape = "Successful inference rate"), size = 4) +
  
  # 右轴
  geom_line(aes(y = MSE_purity_inferred * scaleFactor, color = "MSE purity (successful inference)", linetype = "dashed"), size = 1.5) +
  geom_point(aes(y = MSE_purity_inferred * scaleFactor, color = "MSE purity (successful inference)", shape = "MSE purity (successful inference)"), size = 4) +
  
  geom_line(aes(y = MSE_purity_excluded * scaleFactor, color = "MSE purity (non-main samples)", linetype = "dashed"), size = 1.5) +
  geom_point(aes(y = MSE_purity_excluded * scaleFactor, color = "MSE purity (non-main samples)", shape = "MSE purity (non-main samples)"), size = 4) +
  
  scale_y_continuous(
    name = "Successful inference rate",
    sec.axis = sec_axis(~./scaleFactor, name = "MSE purity")
  ) +
  scale_color_manual(values = c(
    "Successful inference rate" = "#1f77b4",  # 柔和蓝
    "MSE purity (successful inference)" = "#d62728",      # 柔和红
    "MSE purity (non-main samples)" = "#2ca02c"       # 柔和绿
  )) +
  scale_shape_manual(values = c(
    "Successful inference rate" = 16,   # 实心圆
    "MSE purity (successful inference)" = 17,       # 实心三角
    "MSE purity (non-main samples)" = 15        # 实心方块
  )) +
  scale_linetype_identity() +
  theme_classic(base_size = 16) +  # 论文风格
  theme(
    axis.line = element_line(size = 1),
    axis.ticks = element_line(size = 1),
    legend.key = element_blank(),
    legend.position = "top",
    legend.title = element_blank()
  ) +
  labs(
    title = "CSR cutoff vs Metrics",
    x = "CSR cutoff"
  )

ggsave("CSR_cutoff.png", plot = line_plot, width = 8, height = 5, dpi = 300)
```


