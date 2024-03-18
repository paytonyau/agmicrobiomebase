# Step 1: Load the table from the receipt `Webin` after the plant Checklist uploaded to the server
Webin <- read.table("Webin-accessions-2023-12-07T15_42_52.222Z_OR.txt", header = TRUE)

# Remove the last column of the dataframe - `Webin`
Webin <- head(Webin, -1)

# Step 2: Load the data from the Crop Check-list
Check.list <- read.delim('Checklist_GSC-MIxS_16Samplicons_OR_TESTv1.tsv', header = FALSE, sep = "\t")

# Format the crop checklist to fit the fastq checklist requirement and the mapping/merging process.
colnames(Check.list) <- Check.list[2,]
Check.list <- Check.list[-c(1:3),]
Check.list_2 <- Check.list[,c(3:4)]

# Merge 'Webin' and 'Check.list_2' on the columns “ALIAS” and “sample_alias”, respectively, to create ‘merged_df’
merged_df <- merge(Webin, Check.list_2, by.x = "ALIAS", by.y = "sample_alias")
merged_df$sample_title <- sub("-([0-9])([^0-9]|$)", "-0\\1\\2", merged_df$sample_title)

# Step 3: Load the data from the “md5.txt” file
md5 <- read.table("md5.txt", sep = "")
md5 <- md5[, c("V2", "V1")]
md5.R1 <- md5[grep("_R1_001.fastq.gz", md5$V2), ]
md5.R1$v3 <- sapply(strsplit(as.character(md5.R1$V2), "_"), `[`, 1)
md5.R1$v3 <- sub("(\\d)$", "0\\1", md5.R1$v3)
md5.R1$v3 <- sub("-([1-9])-", "-0\\1-", md5.R1$v3)
colnames(md5.R1) <- c("forward_file_name","forward_file_md5", "sample_title")
md5.R2 <- md5[grep("_R2_001.fastq.gz", md5$V2), ]
md5.R2$v3 <- sapply(strsplit(as.character(md5.R2$V2), "_"), `[`, 1)
md5.R2$v3 <- sub("(\\d)$", "0\\1", md5.R2$v3)
md5.R2$v3 <- sub("-([1-9])-", "-0\\1-", md5.R2$v3)
colnames(md5.R2) <- c("reverse_file_name", "reverse_file_md5", "sample_title")
merged.md5 <- merge(md5.R1, md5.R2, by="sample_title")

# Step 4: Merge "merged_df" and "merged.md5" on the "sample_title" column to create "merged_all"
merged_all <- merge(merged_df, merged.md5, by ="sample_title")
merged_all$study <- rep("PRJEB58189", nrow(merged_all))
merged_all$instrument_model <- rep("Illumina MiSeq", nrow(merged_all))
merged_all$library_name <- rep("Nextera XT v2", nrow(merged_all))
merged_all$library_source <- rep("METAGENOMIC", nrow(merged_all))
merged_all$library_selection <- rep("PCR", nrow(merged_all))
merged_all$library_strategy <- rep("AMPLICON", nrow(merged_all))
merged_all$library_layout <- rep("PAIRED", nrow(merged_all))
colnames(merged_all)[which(colnames(merged_all) == "sample")] <- "ACCESSION"
colnames(merged_all)[colnames(merged_all) == 'ACCESSION'] <- 'sample'
merged_all <- merged_all[, c("sample", "study", "instrument_model", "library_name", 
                             "library_source", "library_selection", "library_strategy", 
                             "library_layout", "forward_file_name", "forward_file_md5", 
                             "reverse_file_name", "reverse_file_md5")]

# Step 5: Final modification for the fastq checklist
new_row <- setNames(data.frame(matrix(ncol = ncol(merged_all), nrow = 1)), colnames(merged_all))
new_row[1, c("sample", "study", "instrument_model")] <- c("FileType", "fastq", "Read submission file type")
col_names <- colnames(new_row) 
new_row <- rbind(new_row, col_names)
merged_all <- rbind(new_row, merged_all)
colnames(merged_all) <- NULL

# Step 6: Export the table
# write.table(merged_all, file = "fastq2_template_16Samplicons_OR_TEST_with_mapping.tsv", 
#             sep = "\t", row.names = FALSE, quote = FALSE, na = "")

# Session Info
sessionInfo()
