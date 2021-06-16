setwd("~/Desktop/Jiali/TAMU/xiaohan")

## convert the xlsx file to csv in excel.

# read the annotation file
annot <- read.csv("Annotation-06032021.csv", skip = 4, header = F, stringsAsFactors = F) # it still reads the first two rows, next delete them
colnames(annot) <- annot[3,]
annot <- annot[-c(1:3),]

# read data file
data <- read.csv("New Microsoft Excel Worksheet.csv", header = F, stringsAsFactors = F)
data$V1 <- gsub("seu","Pseu", data$V1)

# merge data and annotation files
data_annot <- merge(data, annot[,c(1:5,18,19,21,22)], by.x = "V1", by.y = "Query", all.x = T)

# write output into table. Cannot write into csv because the contents have comma to separate names.
write.table(data_annot, "results_with_annotation.txt", row.names = F, quote = F, sep = "\t")

# open the table with excel and save as xlsx