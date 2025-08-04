```{r}
# Load required libraries
library(readxl)
library(dplyr)
df <- read_excel("nejmoa1616288_appendix_2.xlsx", sheet = "TableS10")
colnames(df)=df[10,]
df <- df[11:nrow(df),]
df<-as.data.frame(df)
df$nAraw <- as.numeric(df$nAraw)
df$nBraw <- as.numeric(df$nBraw)
df$Ploidy <- as.numeric(df$Ploidy)
df$ACF <- as.numeric(df$ACF)
colnames(df)[colnames(df)=="sample"] <- "samplenames"
colnames(df)[colnames(df)=="endpos"] <- "end"
colnames(df)[colnames(df)=="startpos"] <- "start"
colnames(df)[colnames(df)=="chr"] <- "chromosome"
colnames(df)[colnames(df)=="Ploidy"] <- "ploidies"
colnames(df)[colnames(df)=="ACF"] <- "purities"

# create abs copy number column
df$abs_cn <- df$nAraw+df$nBraw

# transform abs to relative copy number
abs_to_relative <- function(value, purity, ploidy){
  a1 <- value*purity + 2*(1-purity)
  a2 <- purity*ploidy + 2*(1-purity)
  q <- a1/a2
  return(q)
}
df$relative_cn <- abs_to_relative(df$abs_cn,df$purities,df$ploidies)

# extract patient id
df$patient_id <- gsub("-.*", "", df$samplenames)

# create length column
df$start <- as.numeric(df$start)
df$end <- as.numeric(df$end)
df$length <- df$end - df$start
```
```{r}

library(tidyr)
# Transform to wide format
seg_val <- df %>%
  pivot_wider(
    names_from = samplenames,
    values_from = relative_cn,
    id_cols = c(chromosome, start, end, length)
  )

seg_len <- df %>%
  pivot_wider(
    names_from = samplenames,
    values_from = length,
    id_cols = c(chromosome, start, end)
  )

seg_val$chromosome <- as.numeric(seg_val$chromosome)
seg_len$chromosome <- as.numeric(seg_len$chromosome)

# reorder
seg_val <- seg_val %>%
  arrange(chromosome, start)
seg_len <- seg_len %>%
  arrange(chromosome, start)

# delete all NA row
seg_val <- seg_val[rowSums(!is.na(seg_val[, 5:ncol(seg_val)])) > 0, ]

# write.csv(seg_val, "tracerx_seg_vals.csv", row.names = FALSE)
# write.csv(seg_len, "tracerx_seg_len.csv", row.names = FALSE)

```

```{r}
df_ref <- df %>%
  distinct(samplenames, .keep_all = TRUE)

df_ref <- df_ref[,c("patient_id","samplenames","ploidies","purities")]
df_ref$samplenames <- gsub("-", ".", df_ref$samplenames)
write.csv(df_ref, "tracerx_ref.csv", row.names = FALSE)
```

```{r}
count_df <- data.frame(patient_id=character(),samplenames=character(),length=integer())

patient_list <- unique(df$patient_id)

for(i in patient_list){
  samples <- unique(df[df$patient_id==i,"samplenames"])
  for (j in samples){
    count_df <- rbind(count_df,data.frame(
      patient_id = i,
      samplenames = j,
      length = length(df[df$samplenames==j,"samplenames"])
    ))
  }
}
write.csv(count_df, "tracerx_seg_count.csv", row.names = FALSE)
```

