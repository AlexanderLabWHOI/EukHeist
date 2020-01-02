
#Generate assembly group files for metagenomic and metatranscriptomic sample lists

library(dplyr);library(reshape2)
# Import:
metag_in <- read.delim("samplelist-metaG-wgroups.txt", header=T)
metat_in <- read.delim("samplelist-metaT-wgroups.txt", header=T)
metag_in

metag_assemble <- metag_in %>%
    group_by(ASSEMBLY_GROUPING) %>%
    summarize(SAMPLE_LIST=paste(unique(SAMPLEID), collapse=",")) %>%
    as.data.frame
write.table(metag_assemble, "assembly-list-metaG.txt",quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")

metat_assemble <- metat_in %>%
    group_by(ASSEMBLY_GROUPING) %>%
    summarize(SAMPLE_LIST=paste(unique(SAMPLEID), collapse=",")) %>%
    as.data.frame
write.table(metat_assemble, "assembly-list-metaT.txt",quote=FALSE,col.names=TRUE,row.names=FALSE,sep="\t")
