library("dplyr")
library("data.table")
args = commandArgs(trailingOnly=TRUE)

# ## Options
# args <- character(3)
# args[1] = "chr1" #chr
# args[2] = "G_ERV" #donors list
# args[3] = "CRG36"

if (args[3] == "all") {filter = "all"} 
if (args[3] == "CRG36") {filter = "CRG36"} 
	if (args[3] == "CRG36_maskExons") {filter = "CRG36_maskExons"} ## correlate CRG36_maskExons vs CRG36_Exons mutation matrices
	if (args[3] == "CRG36_Exons") {filter = "CRG36_Exons"} 
	# if (args[3] == "CRG36_maskExons" | args[3] == "CRG36_Exons") {filter = "CRG36"} 

if(args[2] == 'UKBB'){
	output_dir <- 'BEDS'
}else{
	output_dir <- 'BEDS_otherdatasets'
}

### Files ###
mutations <- fread(paste0("../results/",output_dir,"/", args[1], "_mutations_transcripts_and_co_FANTOM5_", filter, "_strand5mer.txt"), header = F)
# colnames(mutations) <- c("chr", "start", "end", "REF", "ALT", "cancer_type", "donor", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer", "strand_kmer", "strand_ALT", "strand3mer_strandALT") 
colnames(mutations) <- c("chr", "start", "end", "REF", "ALT", "cancer_type", "donor", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer", "strand_kmer", "strand_ALT") # The files now don't seem to have the last column 

if(args[2] == 'UKBB'){
	mutations_from_list_wo_hyperm <- mutations
}else{
	list_of_donors <- fread(paste0("/users/dweghorn/dcastellano/muCOV/results/Donors_Lists/", args[2],".txt"), header = F)
	colnames(list_of_donors) <- c("donor")
	hypermutated_donors <- fread("/users/dweghorn/dcastellano/muCOV/results/Donors_Lists/hypermutated.txt", header = F)
	colnames(hypermutated_donors) <- c("donor")
	mutations_from_list = inner_join(mutations, list_of_donors, by = c("donor"))
	mutations_from_list_wo_hyperm = anti_join(mutations_from_list, hypermutated_donors, by = c("donor"))
}

write.table(mutations_from_list_wo_hyperm, paste0("../results/",output_dir,"/", args[1], "_", args[2], "_mutations_transcripts_and_co_FANTOM5_", args[3], "_strand5mer_CpG.txt"), col.names = T, quote = F, row.names = F, append = F, sep = "\t") 


### Mutation and site counts by kmer ###
context_mutation = data.frame(table(select(mutations_from_list_wo_hyperm, strand_kmer, strand_ALT))) # <- list of combinations with mutations
context_mutation = subset(context_mutation, Freq > 0)
colnames(context_mutation) <- c("strand_kmer", "strand_ALT","mutation_counts")

# sites <- fread(paste0("../results/",output_dir,"/", args[1], "_transcripts_and_co_FANTOM5_", filter, "_strand5mer_unmapped.txt"), header = F) 
# colnames(sites) <- c("chr", "start", "end", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer", "strand_kmer")

# # Let's create another "sites" with the filtered unmappable regions out
# mappable <- fread(paste0("/nfs/users/dweghorn/mcortes/projects/mutable_TSS_regs/results/goodSitesLiftover/",args[1],"_hg19_lifted.bed"), header = F)
# colnames(mappable) <- c("chr", "start", "end", "ref_hg38", "ref_hg19")
# mappable <- select(mappable, "chr", "start", "end")
# sites = inner_join(sites, mappable, by = c("chr", "start", "end"))
# rm(mappable)
# write.table(sites, paste0("../results/",output_dir,"/", args[1],"_transcripts_and_co_FANTOM5_", filter, "_strand5mer.txt"), col.names = F, quote = F, row.names = F, append = F, sep = "\t") 

# mike comment: just load the good sites already (computed in step 3 of the pipeline)
sites <- fread(paste0("../results/",output_dir,"/", args[1], "_transcripts_and_co_FANTOM5_", filter, "_strand5mer.txt"), header = F) 
colnames(sites) <- c("chr", "start", "end", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer", "strand_kmer")

# Proceed
context_site = data.frame(table(select(sites, strand_kmer))) 
colnames(context_site) <- c("strand_kmer", "site_counts")
context_site$chr = args[1]

preglobal_mutation_matrix = full_join(context_mutation, context_site, by = c("strand_kmer"))


### Output ###
write.table(preglobal_mutation_matrix, paste0("../results/Mutation_Matrices/", args[1], "_", args[2], "_", args[3], ".txt"), col.names = F, quote = F, row.names = F, append = F, sep = "\t") 
