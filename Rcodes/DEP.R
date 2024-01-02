setwd("~/Desktop/Jiali/TAMU/xiaohan/")
library("DEP")
library(dplyr)
library(ggplot2)

# read count data and design table
mydata <- read.csv("Irpex_ms/raw_count.csv", header = T, stringsAsFactors = F)
mydata$ID %>% duplicated() %>% any() # check if protein IDs are unique. 
# output FALSE means there is no duplicated, and it is good to go.
mydata$name <- gsub("\\|.*","",mydata$ID)

mydataDesign <- read.csv("Irpex_ms/Design table.csv", header = T, stringsAsFactors = F)
mydataDesign$label <- colnames(mydata)[2:9]

# Generate a SummarizedExperiment object using an experimental design
data_se <- make_se(mydata, c(2:9), mydataDesign)

# filter genes with missing value
# Filter for proteins that are quantified in all replicates of at least one condition
data_filt <- filter_missval(data_se, thr = 0)
# Filter for proteins that are quantified in at least 1/3 of the samples.
data_filt <- filter_proteins(data_se, "fraction", min = 0.25)
protein_numb_tab <- plot_numbers(data_filt, plot = F)
# plot the number of protein in each run with ggplot2
protein_numb_tab$sample <- factor(protein_numb_tab$sample, levels = c("CNF_1","CNF_2","CL_1","CL_2","CML_1","CML_2","CML.PFAS_1","CML.PFAS_2"))
protein_numb_tab$condition <- factor(protein_numb_tab$condition, levels = c("CNF","CL","CML","CML.PFAS"))
ggplot(protein_numb_tab, aes(x=sample, y=proteins_in_sample, fill=condition))+
  geom_bar(stat = "identity") + #scale_fill_brewer(palette = "Dark2") + 
  geom_text(aes(label=proteins_in_sample), vjust=1.6, color="black",
            position = position_dodge(0.9), size=3.5)+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5))+
  labs(y = "Number of proteins")

#-------impute missing value----------------
# Extract protein names with missing values 
# in all replicates of at least one condition
proteins_MNAR <- get_df_long(data_filt) %>%
  group_by(name, condition) %>%
  dplyr::summarize(NAs = all(is.na(intensity))) %>% 
  filter(NAs) %>% 
  pull(name) %>% 
  unique()

# Get a logical vector
MNAR <- names(data_filt) %in% proteins_MNAR

# Perform a mixed imputation
mixed_imputation <- DEP::impute(
  data_filt, 
  fun = "mixed",
  randna = !MNAR, # we have to define MAR which is the opposite of MNAR
  mar = "knn",# imputation function for MAR
  mnar = "MinProb") # imputation function for MNAR

df <- get_df_wide(mixed_imputation)
write.csv(df, "Irpex_ms/imputed_counts.csv", row.names = F)

## Differential expression
# Normalize the data
data_norm <- normalize_vsn(mixed_imputation)
meanSdPlot(data_norm)

##------------------differential expression test-----------------------
# Test all possible comparisons of samples
data_diff_all_contrasts <- test_diff(data_norm, type = "all")

# Plot the first and second principal components
pcaData <- plot_pca(data_diff_all_contrasts, x = 1, y = 2, n =4867, point_size = 4, plot=T)
# draw PCA with ggplot2
pcaData$condition <- factor(pcaData$condition, levels = c("CNF","CL","CML","CML.PFAS"))
ggplot(pcaData, aes(PC1, PC2, color=condition, shape=replicate)) +
  geom_point(size=4) + #geom_text(aes(label=name),hjust=0, vjust=0)+
  xlab("PC1: 20.7% variance") +
  ylab("PC2: 17% variance") + 
  coord_fixed() +
  scale_color_discrete(breaks=c("CNF","CL","CML","CML.PFAS"))+
  theme_bw(base_size = 14)
ggsave("Irpex_ms/PCA.png", width = 5.5, height = 5.5)


dep <- add_rejections(data_diff_all_contrasts, alpha = 0.1, lfc = log2(1))


# Generate a results table
data_results <- get_results(dep)

# Number of significant proteins
data_results %>% filter(significant) %>% nrow()
# 105 proteins Significant total
length(data_results$CNF_vs_CL_p.adj[data_results$CNF_vs_CL_p.adj< 0.1]) # 22
length(data_results$CNF_vs_CML_p.adj[data_results$CNF_vs_CML_p.adj< 0.1]) # 5
length(data_results$CML_vs_CML.PFAS_p.adj[data_results$CML_vs_CML.PFAS_p.adj< 0.1]) # 36
length(data_results$CNF_vs_CML.PFAS_p.adj[data_results$CNF_vs_CML.PFAS_p.adj< 0.1]) # 65

# output DE protein tables
CNFvCL <- data_results[data_results$CNF_vs_CL_p.adj< 0.1,c(1,12,25)] 
CNFvCML <- data_results[data_results$CNF_vs_CML_p.adj< 0.1,c(1,13,26)]
CMLvPFAS <- data_results[data_results$CML_vs_CML.PFAS_p.adj< 0.1,c(1,11,24)]
CNFvPFAS <- data_results[data_results$CNF_vs_CML.PFAS_p.adj< 0.1,c(1,14,27)]

write.csv(CNFvCL,"DE_CNFvCL.csv", row.names = F)
write.csv(CNFvCML,"DE_CNFvsCML.csv", row.names = F)
write.csv(CMLvPFAS,"DE_CMLvPFAS.csv", row.names = F)

intersect(CNFvCL$name, CNFvCML$name)

## barplot of DE protein numbers
length(CNFvCL$name[CNFvCL$CNF_vs_CL_ratio<0])
length(CNFvCML$name[CNFvCML$CNF_vs_CML_ratio<0])
length(CMLvPFAS$name[CMLvPFAS$CML_vs_CML.PFAS_ratio<0])

x <- data.frame("condition" = c("CNFvsCL", "CNFvsCML", "CMLvsCML+PFAS"), "up" = c(20, 5, 24), "down" = c(2, 0,12))
x_melt <- reshape2::melt(x)
x_melt$condition <- factor(x_melt$condition, levels = c("CNFvsCL", "CNFvsCML", "CMLvsCML+PFAS"))
ggplot(x_melt, aes(x = condition, y=value, fill = variable)) +
  geom_bar(stat = "identity") + theme_classic(base_size = 16)+ 
  labs(fill="", x = "", y="Number of proteins")+
  theme(axis.text.x = element_text(angle = 45,vjust=0.9,hjust=1))

## export DE protein list and functions
DE_all <- data_results[data_results$name %in% c(CNFvCL$name, CNFvCML$name,CMLvPFAS$name,CNFvPFAS$name),]
DE_all_annot <- join(DE_all, KEGG[,c(1,3)], by = "name", match="first")
DE_all_annot <- join(DE_all_annot, KOG[,c(2,4,5)], by="name",match="first")
DE_all_annot <- join(DE_all_annot, IPS[,c(1,3)], by="name", match="first")
write.csv(DE_all_annot[,c(c(1,11:14,24:35))],"DEP_all_update.csv", row.names = F)
