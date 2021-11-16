setwd("~/Desktop/Jiali/TAMU/xiaohan/")

## read and merge count files
filename <- list.files("Irpex_ms", pattern="*count*", full.names = T)
file_list <- filename %>% 
  lapply(read.csv, sep = "\t", stringsAsFactors=F)

protein.id <- do.call(rbind,file_list)[,1]
protein.id <- protein.id[!duplicated(protein.id)] # get all the IDs from the runs and remove the duplicated ones.
protein.id <- protein.id[!(protein.id %in% protein.id[grep("contaminant",protein.id)])] # remove the contaminant IDs

data_df <- data.frame(ID = protein.id)
for (i in 1:length(file_list)) {
  sample.id <- gsub("Irpex_ms/|.count.txt","",filename[i])
  data_df <- merge(data_df, file_list[i], by.x = "ID", by.y = "protein.id", all.x=T)
  colnames(data_df)[i+1] <- sample.id
}
# output the merged counts to a csv file
write.csv(data_df, "Irpex_ms/raw_count.csv", row.names = F)


## create a experiment design table
mydesign <- data.frame(label = colnames(data_df)[-1], condition = "", replicate = "")
write.csv(mydesign, "Irpex_ms/Design table.csv", row.names = F)
# complete the design table in excel