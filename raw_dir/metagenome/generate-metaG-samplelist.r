
library(tidyverse); library(reshape2)

# Import path to detected fastq.gz files
paths <- as.data.frame(list.files(pattern = ".fastq.gz", full.names = FALSE, recursive=TRUE, no.. = TRUE))
colnames(paths)[1]<-"FASTQ"

# Add column with full path
path_to_files <- getwd()
paths$PATH <- path_to_files

samplelist <- colsplit(paths$FASTQ, "_", c("SAMPLEID", "suffix"))
paths_wfastq_wsample <- data.frame(paths, samplelist)
paths_wfastq_wsample$suffix <- NULL
paths_wfastq_wsample$READ <- ifelse(grepl("R1_001.fastq.gz|_1.fastq.gz|R1.fastq.gz", paths_wfastq_wsample$FASTQ), "R1", "R2")

# Widen
paths_wide <- dcast(paths_wfastq_wsample, SAMPLEID + PATH ~ READ, value.var = "FASTQ")

# Assign meta'omic designation and create empty 'assembly group' column
paths_wide$OMIC <- "METAGENOME"
paths_wide$ASSEMBLY_GROUPING <- ""

# write to table
write.table(paths_wide, "samplelist-metaG.txt", quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")


