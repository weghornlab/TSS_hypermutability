library(dplyr)
library(data.table)
library(tidyr)
library(stringr)
args <- commandArgs(trailingOnly=TRUE)
inMutsFile <- args[1]
inTargetFile <- args[2]
inMatrixFile <- args[3]
include_CpGs <- args[4]
outFile <- args[5]

# debug
inMutsFile <- "./exampleData/outputMutations/chr22.tsv.gz"
inTargetFile <- "./exampleData/inputTarget/chr22.tsv.gz"
inMatrixFile <- "./exampleData/outputMatrix/all.tsv.gz"
include_CpGs <- "no"
outFile <- "./exampleData/outputDensity/chr22.tsv.gz"

transcript_mutations <- fread(inMutsFile, nThread = 1)
transcript_mutations$mutID = paste0(transcript_mutations$chr, "_", transcript_mutations$end, "_", transcript_mutations$REF, ">", transcript_mutations$ALT)
transcript_mutations$strand3mer = str_sub(transcript_mutations$strand_kmer,2,-2) # had to add this column
transcript_mutations$strand3mer_strandALT = paste0(transcript_mutations$strand3mer, ">", transcript_mutations$strand_ALT) # had to add this column

all_sites50kb <- fread(inTargetFile, nThread = 1)
all_sites50kb$strand3mer = str_sub(all_sites50kb$strand_kmer,2,-2)
all_sites50kb$strand1mer = str_sub(all_sites50kb$strand_kmer,3,-3)

CpG = data.frame(rbind("ACG>T", "TCG>T", "CCG>T", "GCG>T", "CGA>A", "CGT>A", "CGC>A", "CGG>A"))
colnames(CpG) <- c("strand3mer_strandALT")
CpG$CpG = "Yes"
CpG$strand3mer = c("ACG", "TCG", "CCG", "GCG", "CGA", "CGT", "CGC", "CGG")

### Merge ###
transcript_mutations = left_join(transcript_mutations, CpG, by = c("strand3mer_strandALT", "strand3mer")) # had to modify this
transcript_mutations$CpG[is.na(transcript_mutations$CpG)] <- "No"

all_sites50kb = left_join(all_sites50kb, CpG, by = c("strand3mer"))
all_sites50kb$CpG[is.na(all_sites50kb$CpG)] <- "No"


# mut matrix
preglobal_mutation_matrix <- fread(inMatrixFile, nThread = 1)
colnames(preglobal_mutation_matrix) = c("strand_kmer", "strand_ALT","mutation_counts", "site_counts")
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
	#transcript_mutations = transcript_mutations
	#preglobal_mutation_matrix = preglobal_mutation_matrix
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

global_mutation_matrix <- preglobal_mutation_matrix[
	,
	list("mutation_counts" = sum(mutation_counts), "site_counts" = sum(site_counts)),
	by = c("strand_kmer", "strand_ALT")
]
global_mutation_matrix[, "mutation_rate" := mutation_counts / site_counts]

# clean some memory before continuing
rm(list = ls()[-which(ls() %in% c("outFile", "all_sites50kb", "transcript_mutations", "global_mutation_matrix", "CHR", "args", "mutations_dataset"))]); gc()
all_sites50kb[, c("start", "end", "strand", "kmer", "strand3mer", "CpG") := NULL] # release some memory
all_sites50kb <- all_sites50kb[!grepl("N", strand_kmer)]
transcript_mutations <- transcript_mutations[!grepl("N", strand_kmer)]
transcript_mutations$strand_REF = substr(transcript_mutations$strand_kmer, 3, 3)
transcript_mutations$mutation_strand = paste0(transcript_mutations$strand_REF,">",transcript_mutations$strand_ALT)
transcript_mutations[, c("start", "end", "strand", "kmer", "strand3mer", "strand3mer_strandALT", "CpG", "strand_REF", "strand_ALT", "strand_kmer", "REF", "ALT") := NULL] # release some memory

#######
pmuts <- list("C" = c("A", "G", "T"), "T" = c("A", "C", "G"), "G" = c("A", "C", "T"), "A" = c("C", "G", "T"))[all_sites50kb$strand1mer]
pmuts <- unlist(pmuts, recursive = FALSE, use.names = FALSE)
all_sites50kb <- all_sites50kb[rep(1:nrow(all_sites50kb), each = 3)]
all_sites50kb[, strand_ALT := pmuts] # add column in a more memory efficient way
rm(pmuts); gc()
all_sites50kb[global_mutation_matrix, "mutation_rate" := i.mutation_rate, on = c("strand_kmer" = "strand_kmer", "strand_ALT" = "strand_ALT")]
rm(global_mutation_matrix); gc()
all_sites50kb[is.na(mutation_rate), mutation_rate := 0.0]

# X_Y_mutability
sdf <- expand.grid(strand1mer = c("A", "C", "G", "T"), strand_ALT = c("A", "C", "G", "T"))
sdf <- sdf[sdf$strand1mer != sdf$strand_ALT, ]
sdf$cname <- paste0(sdf$strand1mer, "_", sdf$strand_ALT, "_mutability")
for (i in 1:nrow(sdf)) {

	r <- sdf[i, "strand1mer"]
	a <- sdf[i, "strand_ALT"]
	.c <- sdf[i, "cname"]

	all_sites50kb[, X_Y_mutability := 0.0]
	all_sites50kb[strand1mer == (r) & strand_ALT == (a), X_Y_mutability := mutation_rate]
	tmpdt <- all_sites50kb[, list("X_Y_mutability" = sum(X_Y_mutability)), by = c("chr", "geneID", "position", "bin")]
	
	if (i == 1) {

		xydt <- tmpdt

	} else {

		xydt[tmpdt, "X_Y_mutability" := i.X_Y_mutability, on = c("chr", "geneID", "position", "bin")]

	}
	names(xydt)[ncol(xydt)] <- .c

	rm(tmpdt); gc()

}
all_sites50kb[, "X_Y_mutability" := NULL]

# X content
for (x in c("A", "C", "G", "T")) {
	
	all_sites50kb[, X_content := 0L]
	all_sites50kb[strand1mer == (x), X_content := 1L]
	tmpdt <- all_sites50kb[, list("X_content" = sum(X_content) / 3L), by = c("chr", "geneID", "position", "bin")]
	
	xydt[tmpdt, "X_content" := i.X_content, on = c("chr", "geneID", "position", "bin")]
	names(xydt)[ncol(xydt)] <- paste0(x, "_content")

	rm(tmpdt); gc()

}
rm(all_sites50kb); gc()

# total mutability
xydt[, "transcript_mutability" := C_A_mutability + G_A_mutability + T_A_mutability + A_C_mutability + G_C_mutability + T_C_mutability + A_G_mutability + C_G_mutability + T_G_mutability + A_T_mutability + C_T_mutability + G_T_mutability]
all_sites50kb_AGG <- xydt

# Now the mutations

# X_Y_mutations
transcript_mutations$A_G_mutations <- ifelse(transcript_mutations$mutation_strand == "A>G", 1L, 0L)
transcript_mutations$A_T_mutations <- ifelse(transcript_mutations$mutation_strand == "A>T", 1L, 0L)
transcript_mutations$C_G_mutations <- ifelse(transcript_mutations$mutation_strand == "C>G", 1L, 0L)
transcript_mutations$C_T_mutations <- ifelse(transcript_mutations$mutation_strand == "C>T", 1L, 0L)
transcript_mutations$G_A_mutations <- ifelse(transcript_mutations$mutation_strand == "G>A", 1L, 0L)
transcript_mutations$G_C_mutations <- ifelse(transcript_mutations$mutation_strand == "G>C", 1L, 0L)
transcript_mutations$G_T_mutations <- ifelse(transcript_mutations$mutation_strand == "G>T", 1L, 0L)
transcript_mutations$T_C_mutations <- ifelse(transcript_mutations$mutation_strand == "T>C", 1L, 0L)
transcript_mutations$T_G_mutations <- ifelse(transcript_mutations$mutation_strand == "T>G", 1L, 0L)
transcript_mutations$A_C_mutations <- ifelse(transcript_mutations$mutation_strand == "A>C", 1L, 0L)
transcript_mutations$C_A_mutations <- ifelse(transcript_mutations$mutation_strand == "C>A", 1L, 0L)
transcript_mutations$T_A_mutations <- ifelse(transcript_mutations$mutation_strand == "T>A", 1L, 0L)

transcript_mutations_AGG <- transcript_mutations[, .(total_mutations = .N,
								A_G_mutations = sum(A_G_mutations), A_T_mutations = sum(A_T_mutations), C_G_mutations = sum(C_G_mutations),
								C_T_mutations = sum(C_T_mutations), G_A_mutations = sum(G_A_mutations), G_C_mutations = sum(G_C_mutations),
								G_T_mutations = sum(G_T_mutations), T_C_mutations = sum(T_C_mutations), T_G_mutations = sum(T_G_mutations),
								A_C_mutations = sum(A_C_mutations), C_A_mutations = sum(C_A_mutations), T_A_mutations = sum(T_A_mutations)), by = .(chr, geneID, position, bin)]


sites_mutations <- left_join(all_sites50kb_AGG, transcript_mutations_AGG, by = c("chr", "geneID", "position", "bin"))
sites_mutations <- setcolorder(sites_mutations, c('chr', 'geneID', 'bin', 'transcript_mutability', 'total_mutations', 'position', 'C_T_mutations', 'C_T_mutability', 'G_A_mutations', 'G_A_mutability', 'A_G_mutations', 'A_G_mutability', 'T_C_mutations', 'T_C_mutability', 'G_T_mutations', 'G_T_mutability', 'C_A_mutations', 'C_A_mutability', 'C_G_mutations', 'C_G_mutability', 'G_C_mutations', 'G_C_mutability', 'A_C_mutations', 'A_C_mutability', 'T_G_mutations', 'T_G_mutability', 'A_T_mutations', 'A_T_mutability', 'T_A_mutations', 'T_A_mutability', 'G_content', 'C_content', 'T_content', 'A_content') )

sites_mutations[is.na(total_mutations), total_mutations := 0L]
sites_mutations[is.na(C_T_mutations), C_T_mutations := 0L]
sites_mutations[is.na(G_A_mutations), G_A_mutations := 0L]
sites_mutations[is.na(A_G_mutations), A_G_mutations := 0L]
sites_mutations[is.na(T_C_mutations), T_C_mutations := 0L]
sites_mutations[is.na(G_T_mutations), G_T_mutations := 0L]
sites_mutations[is.na(C_A_mutations), C_A_mutations := 0L]
sites_mutations[is.na(C_G_mutations), C_G_mutations := 0L]
sites_mutations[is.na(G_C_mutations), G_C_mutations := 0L]
sites_mutations[is.na(A_C_mutations), A_C_mutations := 0L]
sites_mutations[is.na(T_G_mutations), T_G_mutations := 0L]
sites_mutations[is.na(A_T_mutations), A_T_mutations := 0L]
sites_mutations[is.na(T_A_mutations), T_A_mutations := 0L]

fwrite(sites_mutations, outFile, nThread = 1)
