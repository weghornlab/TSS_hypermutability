rm(list=ls())
library("dplyr")
library("data.table")
library("tidyr")
library("stringr")
args = commandArgs(trailingOnly=TRUE)

### Options
# args <- c()
#  args[1] = "chr22"
#  args[2] = "no"
#  args[3] = 1000
#  args[4] = "G_ERV"
#  args[5] = "CRG36"

CHR = args[1]
include_CpGs = args[2]

if (args[5] == "all") {filter = "all"} 
if (args[5] == "CRG36") {filter = "CRG36"} 
if (args[5] == "CRG36_maskExons") {filter = "CRG36_maskExons"} ## correlate CRG36_maskExons vs CRG36_Exons mutation matrices
if (args[5] == "CRG36_Exons") {filter = "CRG36_Exons"} 

mutations_dataset = args[4]
mutation_matrix = args[4]

if(args[4] == 'UKBB'){
	output_dir <- 'BEDS'
}else{
	output_dir <- 'BEDS_otherdatasets'
}

### Files ###
transcript_mutations <- fread(paste0("../results/",output_dir,"/",CHR,"_",mutations_dataset,"_mutations_transcripts_and_co_FANTOM5_", args[5], "_strand5mer_CpG.txt"), header = T) 
transcript_mutations$mutID = paste0(transcript_mutations$chr, "_", transcript_mutations$end, "_", transcript_mutations$REF, ">", transcript_mutations$ALT)
transcript_mutations$strand3mer = str_sub(transcript_mutations$strand_kmer,2,-2) # had to add this column
transcript_mutations$strand3mer_strandALT = paste0(transcript_mutations$strand3mer, ">", transcript_mutations$strand_ALT) # had to add this column

all_sites50kb <- fread(paste0("../results/",output_dir,"/", args[1], "_transcripts_and_co_FANTOM5_", args[5], "_strand5mer.txt"), header = F) 
colnames(all_sites50kb) <- c("chr", "start", "end", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer", "strand_kmer")
# nrow(all_sites50kb)
all_sites50kb$strand3mer = str_sub(all_sites50kb$strand_kmer,2,-2)
all_sites50kb$strand1mer = str_sub(all_sites50kb$strand_kmer,3,-3)

CpG = data.frame(rbind("ACG>T", "TCG>T", "CCG>T", "GCG>T", "CGA>A", "CGT>A", "CGC>A", "CGG>A"))
colnames(CpG) <- c("strand3mer_strandALT")
CpG$CpG = "Yes"
CpG$strand3mer = c("ACG", "TCG", "CCG", "GCG", "CGA", "CGT", "CGC", "CGG")

if (args[3] == 1000) {
	transcript_mutations$window_rank = transcript_mutations$bin1000
	all_sites50kb$window_rank = all_sites50kb$bin1000
}

if (args[3] == 100) {
	transcript_mutations$window_rank = transcript_mutations$bin100
	all_sites50kb$window_rank = all_sites50kb$bin100
	transcript_mutations = subset(transcript_mutations, bin1000 <= 5)
	all_sites50kb = subset(all_sites50kb, bin1000 <= 5)
}

if (args[3] == 10) {
	transcript_mutations$window_rank = transcript_mutations$bin10
	all_sites50kb$window_rank = all_sites50kb$bin10
	transcript_mutations = subset(transcript_mutations, bin1000 <= 5)
	all_sites50kb = subset(all_sites50kb, bin1000 <= 5)
}


### Merge ###
# transcript_mutations2 = left_join(transcript_mutations, CpG, by=join_by("strand3mer.x"=="strand3mer")) # had to modify this
transcript_mutations = left_join(transcript_mutations, CpG, by = c("strand3mer_strandALT", "strand3mer")) # had to modify this
transcript_mutations$CpG[is.na(transcript_mutations$CpG)] <- "No"

all_sites50kb = left_join(all_sites50kb, CpG, by = c("strand3mer"))
all_sites50kb$CpG[is.na(all_sites50kb$CpG)] <- "No"

# Mutation matrix
preglobal_mutation_matrix = fread(paste0("../results/Mutation_Matrices/", mutation_matrix,"_", filter, ".txt"), sep = "\t", header = F) 
colnames(preglobal_mutation_matrix) = c("strand_kmer", "strand_ALT","mutation_counts", "site_counts", "chr")
preglobal_mutation_matrix$strand_3mer = substr(preglobal_mutation_matrix$strand_kmer, 2, 4)
preglobal_mutation_matrix$strand3mer_strandALT = paste0(preglobal_mutation_matrix$strand_3mer, ">", preglobal_mutation_matrix$strand_ALT)

# nrow(preglobal_mutation_matrix)
preglobal_mutation_matrix =  subset(preglobal_mutation_matrix, site_counts > 1) #to rmv kmer with N's
# nrow(preglobal_mutation_matrix)
preglobal_mutation_matrix = full_join(preglobal_mutation_matrix, CpG, by = c("strand3mer_strandALT"))
preglobal_mutation_matrix$CpG[is.na(preglobal_mutation_matrix$CpG)] <- "No"
# nrow(preglobal_mutation_matrix)
preglobal_mutation_matrix_all_kmers = preglobal_mutation_matrix #? to count CpG kmers that are not CpG>TpG mutations?


if (include_CpGs == "all") {
	transcript_mutations = transcript_mutations
	preglobal_mutation_matrix = preglobal_mutation_matrix
	# all_sites50kb = all_sites50kb
}

if (include_CpGs == "no") {
	transcript_mutations =  subset(transcript_mutations, CpG == "No")
	preglobal_mutation_matrix =  subset(preglobal_mutation_matrix, CpG == "No")	
	## all_sites50kb = subset(all_sites50kb, CpG == "No")
}

if (include_CpGs == "only") {
	transcript_mutations =  subset(transcript_mutations, CpG == "Yes")
	preglobal_mutation_matrix =  subset(preglobal_mutation_matrix, CpG == "Yes")
	all_sites50kb = subset(all_sites50kb, CpG == "Yes")
}

# total_mutations_mutation_matrix = sum(preglobal_mutation_matrix$mutation_counts)

# scaling_factor = total_mutations_mutation_dataset/total_mutations_mutation_matrix ##this produces NA when G_ERV


# CHR
# include_CpGs

mutations_wide = data.frame(tapply(preglobal_mutation_matrix$mutation_counts, list(preglobal_mutation_matrix$strand_kmer, preglobal_mutation_matrix$strand_ALT), sum))
mutations_wide$strand_kmer = row.names(mutations_wide)
mutations_long <- gather(mutations_wide, strand_ALT, mutation_counts, A:T, factor_key=TRUE)
###mutations_long <-  na.omit(mutations_long)
mutations_long$mutation_counts[is.na(mutations_long$mutation_counts)] <- 0

## table(preglobal_mutation_matrix$strand_kmer) #there are combinations kmer>mutations that do not occur in all chromosomes and even kmers that are not found in some chromosomes
# preglobal_mutation_matrix_unique         = unique(select(preglobal_mutation_matrix,           "strand_kmer", "site_counts", "chr"))
preglobal_mutation_matrix_all_kmers_unique = unique(select(preglobal_mutation_matrix_all_kmers, "strand_kmer", "site_counts", "chr"))

sites = data.frame(tapply(preglobal_mutation_matrix_all_kmers_unique$site_counts, list(preglobal_mutation_matrix_all_kmers_unique$strand_kmer), sum))
# sites = data.frame(tapply(preglobal_mutation_matrix_unique$site_counts, preglobal_mutation_matrix_unique$strand_kmer, sum))
# head(sites)
sites$strand_kmer = row.names(sites)
colnames(sites) <- c("site_counts", "strand_kmer")
## subset(sites, strand_kmer == "CGAAC")
## subset(sites, strand_kmer == "TCGCG")

global_mutation_matrix = inner_join(mutations_long, sites, by = c("strand_kmer"))
global_mutation_matrix$mutation_rate = global_mutation_matrix$mutation_counts / global_mutation_matrix$site_counts
head(global_mutation_matrix)


# global_mutation_matrix$mutation_rate = global_mutation_matrix$mutation_rate*scaling_factor
# print(scaling_factor)

#recurrent = data.frame(table(transcript_mutations$mutID))
#subset(transcript_mutations, mutID == "chr22_17299461_A>G")


#######
all_sites50kb[, c("start", "end", "strand", "bin1000", "bin100", "bin10", "transcript_size", "kmer", "strand3mer", "strand3mer_strandALT", "CpG", "transcript_type") := NULL] # release some memory
pmuts <- list("C" = c("A", "G", "T"), "T" = c("A", "C", "G"), "G" = c("A", "C", "T"), "A" = c("C", "G", "T"))[all_sites50kb$strand1mer]
all_sites50kb <- all_sites50kb[rep(1:nrow(all_sites50kb), each = 3)]
all_sites50kb$strand_ALT <- unlist(pmuts)
rm(pmuts)
all_sites50kb[global_mutation_matrix, "mutation_rate" := i.mutation_rate, on = c("strand_kmer" = "strand_kmer", "strand_ALT" = "strand_ALT")]

# X_Y_mutability
all_sites50kb$A_G_mutability <- ifelse(all_sites50kb$strand1mer == "A" & all_sites50kb$strand_ALT == "G", all_sites50kb$mutation_rate, 0)
all_sites50kb$A_T_mutability <- ifelse(all_sites50kb$strand1mer == "A" & all_sites50kb$strand_ALT == "T", all_sites50kb$mutation_rate, 0)
all_sites50kb$C_G_mutability <- ifelse(all_sites50kb$strand1mer == "C" & all_sites50kb$strand_ALT == "G", all_sites50kb$mutation_rate, 0)
all_sites50kb$C_T_mutability <- ifelse(all_sites50kb$strand1mer == "C" & all_sites50kb$strand_ALT == "T", all_sites50kb$mutation_rate, 0)
all_sites50kb$G_A_mutability <- ifelse(all_sites50kb$strand1mer == "G" & all_sites50kb$strand_ALT == "A", all_sites50kb$mutation_rate, 0)
all_sites50kb$G_C_mutability <- ifelse(all_sites50kb$strand1mer == "G" & all_sites50kb$strand_ALT == "C", all_sites50kb$mutation_rate, 0)
all_sites50kb$G_T_mutability <- ifelse(all_sites50kb$strand1mer == "G" & all_sites50kb$strand_ALT == "T", all_sites50kb$mutation_rate, 0)
all_sites50kb$T_C_mutability <- ifelse(all_sites50kb$strand1mer == "T" & all_sites50kb$strand_ALT == "C", all_sites50kb$mutation_rate, 0)
all_sites50kb$T_G_mutability <- ifelse(all_sites50kb$strand1mer == "T" & all_sites50kb$strand_ALT == "G", all_sites50kb$mutation_rate, 0)
all_sites50kb$A_C_mutability <- ifelse(all_sites50kb$strand1mer == "A" & all_sites50kb$strand_ALT == "C", all_sites50kb$mutation_rate, 0)
all_sites50kb$C_A_mutability <- ifelse(all_sites50kb$strand1mer == "C" & all_sites50kb$strand_ALT == "A", all_sites50kb$mutation_rate, 0)
all_sites50kb$T_A_mutability <- ifelse(all_sites50kb$strand1mer == "T" & all_sites50kb$strand_ALT == "A", all_sites50kb$mutation_rate, 0)

# X_content
all_sites50kb$G_content <- ifelse(all_sites50kb$strand1mer == "G", 1/3, 0) # we have to put 1/3 because each position is repeated 3 times, so when we sum we have the correct number
all_sites50kb$C_content <- ifelse(all_sites50kb$strand1mer == "C", 1/3, 0)
all_sites50kb$T_content <- ifelse(all_sites50kb$strand1mer == "T", 1/3, 0)
all_sites50kb$A_content <- ifelse(all_sites50kb$strand1mer == "A", 1/3, 0)

all_sites50kb_AGG <- all_sites50kb[, .(transcript_mutability = sum(mutation_rate), A_G_mutability = sum(A_G_mutability), A_T_mutability = sum(A_T_mutability), C_G_mutability = sum(C_G_mutability),
								C_T_mutability = sum(C_T_mutability), G_A_mutability = sum(G_A_mutability), G_C_mutability = sum(G_C_mutability),
								G_T_mutability = sum(G_T_mutability), T_C_mutability = sum(T_C_mutability), T_G_mutability = sum(T_G_mutability),
								A_C_mutability = sum(A_C_mutability), C_A_mutability = sum(C_A_mutability), T_A_mutability = sum(T_A_mutability),
								G_content = sum(G_content), C_content = sum(C_content), T_content = sum(T_content), A_content = sum(A_content)), by = .(chr, geneID, position, window_rank)]

# Now the mutations
transcript_mutations$strand_REF = substr(transcript_mutations$strand_kmer, 3, 3)
transcript_mutations$mutation_strand = paste0(transcript_mutations$strand_REF,">",transcript_mutations$strand_ALT)
transcript_mutations[, c("start", "end", "strand", "bin1000", "bin100", "bin10", "transcript_size", "kmer", "strand3mer", "strand3mer_strandALT", "CpG", "cancer_type", "strand_REF", "strand_ALT", "strand_kmer", "REF", "ALT", "transcript_type") := NULL] # release some memory

# X_Y_mutations
transcript_mutations$A_G_mutations <- ifelse(transcript_mutations$mutation_strand == "A>G", 1, 0)
transcript_mutations$A_T_mutations <- ifelse(transcript_mutations$mutation_strand == "A>T", 1, 0)
transcript_mutations$C_G_mutations <- ifelse(transcript_mutations$mutation_strand == "C>G", 1, 0)
transcript_mutations$C_T_mutations <- ifelse(transcript_mutations$mutation_strand == "C>T", 1, 0)
transcript_mutations$G_A_mutations <- ifelse(transcript_mutations$mutation_strand == "G>A", 1, 0)
transcript_mutations$G_C_mutations <- ifelse(transcript_mutations$mutation_strand == "G>C", 1, 0)
transcript_mutations$G_T_mutations <- ifelse(transcript_mutations$mutation_strand == "G>T", 1, 0)
transcript_mutations$T_C_mutations <- ifelse(transcript_mutations$mutation_strand == "T>C", 1, 0)
transcript_mutations$T_G_mutations <- ifelse(transcript_mutations$mutation_strand == "T>G", 1, 0)
transcript_mutations$A_C_mutations <- ifelse(transcript_mutations$mutation_strand == "A>C", 1, 0)
transcript_mutations$C_A_mutations <- ifelse(transcript_mutations$mutation_strand == "C>A", 1, 0)
transcript_mutations$T_A_mutations <- ifelse(transcript_mutations$mutation_strand == "T>A", 1, 0)

transcript_mutations_AGG <- transcript_mutations[, .(total_mutations = .N, donors_mutated = paste0(donor,collapse=", "),
								A_G_mutations = sum(A_G_mutations), A_T_mutations = sum(A_T_mutations), C_G_mutations = sum(C_G_mutations),
								C_T_mutations = sum(C_T_mutations), G_A_mutations = sum(G_A_mutations), G_C_mutations = sum(G_C_mutations),
								G_T_mutations = sum(G_T_mutations), T_C_mutations = sum(T_C_mutations), T_G_mutations = sum(T_G_mutations),
								A_C_mutations = sum(A_C_mutations), C_A_mutations = sum(C_A_mutations), T_A_mutations = sum(T_A_mutations)), by = .(chr, geneID, position, window_rank)]


sites_mutations <- left_join(all_sites50kb_AGG, transcript_mutations_AGG, by = c("chr", "geneID", "position", "window_rank"))
sites_mutations <- setcolorder(sites_mutations, c('chr', 'geneID', 'window_rank', 'transcript_mutability', 'total_mutations', 'donors_mutated', 'position', 'C_T_mutations', 'C_T_mutability', 'G_A_mutations', 'G_A_mutability', 'A_G_mutations', 'A_G_mutability', 'T_C_mutations', 'T_C_mutability', 'G_T_mutations', 'G_T_mutability', 'C_A_mutations', 'C_A_mutability', 'C_G_mutations', 'C_G_mutability', 'G_C_mutations', 'G_C_mutability', 'A_C_mutations', 'A_C_mutability', 'T_G_mutations', 'T_G_mutability', 'A_T_mutations', 'A_T_mutability', 'T_A_mutations', 'T_A_mutability', 'G_content', 'C_content', 'T_content', 'A_content') )

sites_mutations[is.na(sites_mutations)] <- 0

write.table(sites_mutations, file = paste0("../results/Expected_mutations_per_bin/Chr_parallel/",CHR,"_",mutations_dataset,"_obs_and_exp_mutations_",args[3],"_bin_size_",args[2],"_CpGs_", args[5],".txt"), col.names = F, row.names = F, quote = F, sep = "\t", append = F)