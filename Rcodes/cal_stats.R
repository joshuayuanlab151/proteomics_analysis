setwd("~/Desktop/Jiali/TAMU/proteomics pipeline/")

# read the spetral counts from different methods
meth1 <- read.csv("prolucid_pattlab.csv", header = T, stringsAsFactors = F)
meth1_target <- meth1[!(meth1$Locus %in% meth1$Locus[grep("Reverse",meth1$Locus)]),] # remove the reversed IDs
meth2 <- read.table("prolucid_percolator.txt", header = T, sep = "\t")
meth3 <- read.table("crux_lowRes.txt", header = T, sep = "\t")
meth4 <- read.table("msgf_lowRes.txt", header = T, sep = "\t")

# Venn diagram
library(VennDiagram)
library(RColorBrewer)
myCol <- brewer.pal(4, "Dark2")

venn.diagram(
  x = list( meth1_target$Locus, meth2$protein.id, meth3$protein.id, meth4$protein.id),
  category.names = c("Prolu+Pattern" , "Prolu+percolator" , "Crux+percolator","MSGF+percolator"),
  filename = 'venn.png',
  output = TRUE ,
  imagetype="png" ,
  height = 580 , 
  width = 630 , 
  resolution = 300,
  compression = "lzw",
  lwd = 1,
  col=myCol,
  fill = myCol,
  cex = 0.5,
  fontfamily = "sans",
  cat.cex = 0.3,
  cat.default.pos = "outer"
)

## Correlation matrix
# merge tables
mydata <- merge(meth1_target[,c(1,9)], meth2, by.x = "Locus", by.y = "protein.id", all = T)
mydata <- merge(mydata, meth3 , by.x = "Locus", by.y = "protein.id", all = T)
mydata <- merge(mydata, meth4 , by.x = "Locus", by.y = "protein.id", all = T)
names(mydata) <- c("Locus","prolu+patternlab","prolu+percolator","crux+percolator","msgf+percolator")

mydata[is.na(mydata)] <- 0 # convert NA to zero

rownames(mydata) <- mydata$Locus
mydata_mat <- as.matrix(mydata[,-1])

# install.packages("PerformanceAnalytics") # install package, only run it once.

library("PerformanceAnalytics")
chart.Correlation(mydata_mat, histogram=T, pch=12)
