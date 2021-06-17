setwd("~/Desktop/Jiali/TAMU/proteomics pipeline/")
## install DEP, only need to run once.
if (!requireNamespace("BiocManager", quietly=TRUE))
  install.packages("BiocManager")
BiocManager::install("DEP")

library("DEP")
run_app("LFQ") # open an interactive window to run the analysis with GUI.

## run tutoral dataset
library("dplyr")

data <- UbiLength
data <- filter(data, Reverse != "+", Potential.contaminant != "+")
data$Gene.names %>% duplicated() %>% any()

data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
# give each ID a unique name
data$name %>% duplicated() %>% any()
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

# Generate a SummarizedExperiment object using an experimental design
LFQ_columns <- grep("LFQ.", colnames(data_unique)) # get LFQ column numbers
experimental_design <- UbiLength_ExpDesign
data_se <- make_se(data_unique, LFQ_columns, experimental_design)

# Let's have a look at the SummarizedExperiment object
data_se

# Plot a barplot of the protein identification overlap between samples
plot_frequency(data_se)
# Filter for proteins that are identified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 1)
plot_numbers(data_filt)

# Normalize the data
data_norm <- normalize_vsn(data_filt)
plot_normalization(data_filt, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt)

### There are different methods can be used to impute missing values, details see https://bioconductor.org/packages/3.13/bioc/vignettes/MSnbase/inst/doc/v01-MSnbase-demo.html#81_Data_imputation
# MNAR features should ideally be imputed with a left-censor method, for MAR, it is recommended to use knn or MLE

# Impute missing data using random draws from a Gaussian distribution centered around a minimal value (for MNAR)
data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)

# Impute missing data using random draws from a manually defined left-shifted Gaussian distribution (for MNAR)
data_imp_man <- impute(data_norm, fun = "man", shift = 1.8, scale = 0.3)

# Impute missing data using the k-nearest neighbour approach (for MAR)
data_imp_knn <- impute(data_norm, fun = "knn", rowmax = 0.9)

plot_imputation(data_norm, data_imp_knn)

# differential expression analysis
# Test every sample versus control
data_diff <- test_diff(data_imp, type = "control", control = "Ctrl")
dep <- add_rejections(data_diff, alpha = 0.05, lfc = log2(1))
plot_pca(dep, x = 1, y = 2, n = 2185, point_size = 4)

# Plot a heatmap of all significant proteins with the data centered per protein
plot_heatmap(dep, type = "centered", kmeans = TRUE, 
             k = 6, col_limit = 4, show_row_names = FALSE,
             indicate = c("condition", "replicate"))

plot_volcano(dep, contrast = "Ubi1_vs_Ctrl", label_size = 2, add_names = TRUE)

df <- dep@elementMetadata
