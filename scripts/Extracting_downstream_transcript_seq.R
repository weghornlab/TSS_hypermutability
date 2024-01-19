library("dplyr")
library("data.table")
library("stringr")
args = commandArgs(trailingOnly=TRUE)
# args[1] = "chr22"
options(scipen=999)

set.seed(8)
     
### Dataset ###
all_chrs = fread("/users/dweghorn/dcastellano/muCOV/data/Promoterome/5BFANTOMCAT5DRobustgene.osc", header = T)
# all_chrs$length = transcript$`eedb:end` - transcript$`eedb:start.0base`
# hist(all_chrs$length, breaks = 100)
# summary(all_chrs$length)

# for (CHR in c("chr1", "chr2", "chr3", "chr4", "chr5", "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21", "chr22")) 
# {

  CHR = args[1]  
  transcript = subset(all_chrs, `eedb:chrom` == CHR)
  
  for (i in seq(nrow(transcript))) 
  {
  #i = 529
    options(scipen=999)
    row = slice(transcript, i)
    results = data.frame()
  #i
    row$transcript_type = str_split(row$`eedb:name`, "\\|")[[1]][3]
    row$geneID = str_split(row$`eedb:name`, "\\|")[[1]][4]
    
    if(row$`eedb:strand` == "+") {
      results = data.frame(seq(row$`eedb:start.0base`, row$`eedb:end`, 1))
      colnames(results) <- c("position")
      results$chr = row$`eedb:chrom`
      results$start = results$position - 1
      results$stop = results$position + 0
      results$geneID = row$geneID
      results$transcript = row$transcript_type
      results$strand = row$`eedb:strand`
      results$bin1000 = rep(1:ceiling(nrow(results)/as.numeric(1000)), each=as.numeric(1000), length.out=nrow(results))
      results$bin100  = rep(1:ceiling(nrow(results)/as.numeric(100)),  each=as.numeric(100),  length.out=nrow(results))
      results$bin10  = rep(1:ceiling(nrow(results)/as.numeric(10)),  each=as.numeric(10),  length.out=nrow(results))
      
      results1 = select(results, "chr", "start", "stop", "transcript", "geneID", "strand", "bin1000", "bin100", "bin10")
      write.table(results1, file = paste0("../results/BEDS/", CHR, "_transcripts_downstream_FANTOM5.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)
    }    
    
  #i = 32181
    
    if(row$`eedb:strand` == "-") {
      
      results = data.frame(seq(row$`eedb:start.0base`, row$`eedb:end`, 1))
      colnames(results) <- c("position")
      results$chr = row$`eedb:chrom`
      results$start = results$position - 1
      results$stop = results$position + 0
      results$geneID = row$geneID
      results$transcript = row$transcript_type
      results$strand = row$`eedb:strand`
      results$bin1000 = rev(rep(1:ceiling(nrow(results)/as.numeric(1000)), each=as.numeric(1000), length.out=nrow(results)))
      results$bin100  = rev(rep(1:ceiling(nrow(results)/as.numeric(100)),  each=as.numeric(100),  length.out=nrow(results)))
      results$bin10  = rev(rep(1:ceiling(nrow(results)/as.numeric(10)),  each=as.numeric(10),  length.out=nrow(results)))
      
      results1 = select(results, "chr", "start", "stop", "transcript", "geneID", "strand", "bin1000", "bin100", "bin10")
      write.table(results1, file = paste0("../results/BEDS/", CHR, "_transcripts_downstream_FANTOM5.bed"), quote = F, row.names = F, col.names = F, sep = "\t", append = T)
    }    
  }
# }
