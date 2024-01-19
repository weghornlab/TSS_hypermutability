library("dplyr")
library("data.table")
library("tidyr")
library("ggplot2")
library("Biostrings")
args = commandArgs(trailingOnly=TRUE)

# args[1] = "chr21"
# args[2] = "CRG36_Exons"

### Datasets ###
dataset <- args[3]
# dataset <- 'UKBB'
# dataset <- 'others'
cat("LOADING FILES...\n")
if(dataset == 'UKBB'){
    UKBB_mutations <- fread(paste0("../data/Mutations/UKBB/UKBB_hg19.bed"), header = F) 
    colnames(UKBB_mutations) <- c("chr", "start", "end", "REF", "ALT", "interval", "mutation_type", "unknown_col")
    nrow(UKBB_mutations)
    UKBB_mutations = subset(UKBB_mutations, chr == args[1])
    UKBB_mutations$cancer_type = "Polymorphisms"
    UKBB_mutations$donor = "Rare"
    UKBB_mutations_short = select(UKBB_mutations, "chr", "start", "end", "REF", "ALT", "cancer_type", "donor") 
    ### Merge Datasets ###
    mutations_short = UKBB_mutations_short
    output_dir = 'BEDS'
}else{
    somatic_mutations  <- fread("../data/Mutations/PCAWG_mutations.txt", header = T) 
    colnames(somatic_mutations) <- c("chr", "start", "end", "Variant_Classification", "Hugo_Symbol", "REF", "Tumor_Seq_Allele1", "ALT", "i_NumCallers", "t_alt_count", "t_ref_count", "cancer_type", "donor", "TCGA")
    somatic_mutations$chr = paste0("chr", somatic_mutations$chr)
    nrow(somatic_mutations)
    somatic_mutations = subset(somatic_mutations, chr == args[1])
    nrow(somatic_mutations)
    somatic_mutations$start = somatic_mutations$start - 1

    germline_mutations <- fread("../data/Mutations/germinal_ultimate_dataset.bed", header = F) 
    colnames(germline_mutations) <- c("chr", "start", "end", "REF", "ALT", "dataset", "mutation_type", "condition")
    nrow(germline_mutations)
    germline_mutations = subset(germline_mutations, chr == args[1])
    nrow(germline_mutations)

    # I needed to copy this and get chr1 file from David because of permission's issues
    ERV_mutations <- fread(paste0("../data/Mutations/interval1_hg19/interval1_mutations_hg19.bed"), header = F) 
    colnames(ERV_mutations) <- c("chr", "start", "end", "REF", "ALT", "interval", "mutation_type", "unknown_col")
    nrow(ERV_mutations)
    ERV_mutations = subset(ERV_mutations, chr == args[1])
    nrow(ERV_mutations)

    Common_mutations <- fread(paste0("../data/Mutations/interval5_hg19/interval5_mutations_hg19.bed"), header = F) 
    colnames(Common_mutations) <- c("chr", "start", "end", "REF", "ALT", "interval", "mutation_type", "unknown_col")
    nrow(Common_mutations)
    Common_mutations = subset(Common_mutations, chr == args[1])
    nrow(Common_mutations)

    all_mutations <- fread(paste0("../data/Mutations/all_intervals/gnomAD_hg19.bed"), header = F) ## except interval 1 and 5
    colnames(all_mutations) <- c("chr", "start", "end", "REF", "ALT", "interval", "mutation_type", "unknown_col")
    nrow(all_mutations)
    all_mutations = subset(all_mutations, chr == args[1])
    nrow(all_mutations)

    germline_mutations$cancer_type = "Germline"
    germline_mutations$donor = "Trios"

    ERV_mutations$cancer_type = "Polymorphisms"
    ERV_mutations$donor = ERV_mutations$interval

    Common_mutations$cancer_type = "Polymorphisms"
    Common_mutations$donor = Common_mutations$interval

    all_mutations$cancer_type = "Polymorphisms"
    all_mutations$donor = all_mutations$interval

    ### Although here common, all and ERV are loaded independently in fact all does not include common and ERV. So in practice I'm only loading all and later on I select the specific freq interval.

    somatic_mutations_short = unique(select(somatic_mutations, "chr", "start", "end", "REF", "ALT", "cancer_type", "donor")) #there are somatic mutations sequenced multiple times in one patient, maybe through different sequencing platforms
    germline_mutations_short = select(germline_mutations, "chr", "start", "end", "REF", "ALT", "cancer_type", "donor") #there are germline mutations that are recurrent across trios
    ERV_mutations_short = select(ERV_mutations, "chr", "start", "end", "REF", "ALT", "cancer_type", "donor") 
    Common_mutations_short = select(Common_mutations, "chr", "start", "end", "REF", "ALT", "cancer_type", "donor") 
    all_mutations_short = select(all_mutations, "chr", "start", "end", "REF", "ALT", "cancer_type", "donor") 

    ### Merge Datasets ###
    mutations_short = rbind(somatic_mutations_short, germline_mutations_short, ERV_mutations_short, Common_mutations_short, all_mutations_short)  
    output_dir = 'BEDS_otherdatasets'
}

# sites <- fread(paste0("../results/BEDS/", args[1], "_transcripts_and_co_FANTOM5_", args[2], "_5mer.txt"), header = F, fill=TRUE)
# colnames(sites) <- c("chr", "start", "end", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer")

sites <- fread(paste0("../results/BEDS_otherdatasets_old/", args[1], "_transcripts_and_co_FANTOM5_", args[2], "_strand5mer.txt"), header = F, fill=TRUE)
colnames(sites) <- c("chr", "start", "end", "transcript_type", "geneID", "strand", "bin1000", "bin100", "bin10", "transcript_size", "position", "kmer", "strand_kmer")
## here I added the output of this script as input because I deleted the input before running with the script with the rest of common and all gnomAD mutations.

# # Get "mappable" sites obtained from liftover
# mappable <- fread(paste0("/nfs/users/dweghorn/mcortes/projects/mutable_TSS_regs/results/goodSitesLiftover/",args[1],"_hg19_lifted.bed"), header = F)
# colnames(mappable) <- c("chr", "start", "end", "ref_hg38", "ref_hg19")
# mappable <- select(mappable, "chr", "start", "end")

# Get bad liftover sites
bad <- fread(paste0("/nfs/users/dweghorn/mcortes/projects/mutable_TSS_regs/results/badSitesLiftover/",args[1],"_hg19_lifted.bed"), header = F)
colnames(bad) <- c("chr", "start", "end")

cat("JOINING DATA...\n")
# mutations_sites = inner_join(mutations_short, sites, by = c("chr", "start", "end"))
# mutations_sites = inner_join(mutations_sites, mappable, by = c("chr", "start", "end"))
# comment by mike: this is causing problems, lets go for a more transparent filtering in base R
sites <- sites[!(end %in% bad$end)]
rm(bad); gc()
mutations_short[
    sites,
    `:=` (
        transcript_type	= i.transcript_type,
        geneID = i.geneID,
        strand = i.strand,
        bin1000	= i.bin1000,
        bin100 = i.bin100,
        bin10 = i.bin10,
        transcript_size	= i.transcript_size,
        position = i.position,
        kmer = i.kmer
    ),
    on = "end"
]
mutations_short <- mutations_short[!is.na(kmer)]

mutations_sites <- mutations_short
### Strand Wise Kmers & Mutations ###
mutations_sites_negative_strand <- subset(mutations_sites, strand == "-")
sites_negative_strand <- subset(sites, strand == "-")
mutations_sites_positive_strand <- subset(mutations_sites, strand == "+")
sites_positive_strand <- subset(sites, strand == "+")

#mutations_sites_negative_strand$strand_kmer <- sapply(mutations_sites_negative_strand$kmer, function(x) as.character(reverseComplement(DNAString(x))))
#mutations_sites_negative_strand$strand_ALT <- sapply(mutations_sites_negative_strand$ALT, function(x) as.character(reverseComplement(DNAString(x))))
#sites_negative_strand$strand_kmer <- sapply(sites_negative_strand$kmer, function(x) as.character(reverseComplement(DNAString(x))))

cat("MAGIC...\n")
mutations_sites_negative_strand$strand_kmer <- as.character(reverseComplement(DNAStringSet(mutations_sites_negative_strand$kmer)))
mutations_sites_negative_strand$strand_ALT <- as.character(reverseComplement(DNAStringSet(mutations_sites_negative_strand$ALT)))
sites_negative_strand$strand_kmer <- as.character(reverseComplement(DNAStringSet(sites_negative_strand$kmer)))

mutations_sites_positive_strand$strand_kmer <- mutations_sites_positive_strand$kmer
mutations_sites_positive_strand$strand_ALT <- mutations_sites_positive_strand$ALT
sites_positive_strand$strand_kmer <- sites_positive_strand$kmer

mutations_sites2 = rbind(mutations_sites_negative_strand, mutations_sites_positive_strand)
sites2 = rbind(sites_negative_strand, sites_positive_strand)

# # We remove the mutations where REF doesn't match the middle base from the trimer
# mutations_sites2$kmer_ref <- substr(mutations_sites2[[16]],3,3)
# mutations_sites2 = mutations_sites2[kmer_ref == REF]
# mutations_sites2[, ("kmer_ref") := NULL]

### Outputs ###
cat("WRITING RESULT...\n")
write.table(sites2, paste0("../results/",output_dir,"/", args[1], "_transcripts_and_co_FANTOM5_", args[2], "_strand5mer.txt"), col.names = F, quote = F, row.names = F, sep = "\t", append = F) 
write.table(mutations_sites2, file = paste0("../results/",output_dir,"/", args[1], "_mutations_transcripts_and_co_FANTOM5_", args[2], "_strand5mer.txt"), quote = F, row.names = F, col.names = F, sep = "\t", append = F)

cat("DONE\n")

# table(mutations_sites$transcript_type)
# inspected_muts = subset(mutations_sites, transcript_type == "coding_mRNA")
# inspected_sites = subset(sites, transcript_type == "coding_mRNA") 

# mut_bin_pos = data.frame(tapply(inspected_muts$start, list(inspected_muts$bin1000, inspected_muts$position), length))
# mut_bin_pos$bin <- rownames(mut_bin_pos)
# mut_bin_pos_long <- gather(mut_bin_pos, position, counts, downstream:upstream, factor_key=TRUE)
# mut_bin_pos_long$bin = as.numeric(as.character(mut_bin_pos_long$bin))

# sites_bin_pos = data.frame(tapply(inspected_sites$start, list(inspected_sites$bin1000, inspected_sites$position), length))
# sites_bin_pos$bin <- rownames(sites_bin_pos)
# sites_bin_pos_long <- gather(sites_bin_pos, position, counts, downstream:upstream, factor_key=TRUE)
# sites_bin_pos_long$bin = as.numeric(as.character(sites_bin_pos_long$bin))

# sites_muts_bin_pos_long = inner_join(mut_bin_pos_long, sites_bin_pos_long, by = c("bin", "position")) 
# sites_muts_bin_pos_long$mu = sites_muts_bin_pos_long$counts.x / sites_muts_bin_pos_long$counts.y

#ggplot(sites_muts_bin_pos_long, aes(x=bin, y=mu, shape=position, color=position)) +
#  geom_point()



