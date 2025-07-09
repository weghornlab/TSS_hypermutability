library(data.table)
library(Biostrings)
library(R.utils)
args = commandArgs(trailingOnly = TRUE)
mutIn <- args[1]
targetIn <- args[2]
matrixOut <- args[3]
mutOut <- args[4]

# # debug
# mutIn <- "./exampleData/inputMutations/chr22.tsv.gz"
# targetIn <- "./exampleData/inputTarget/chr22.tsv.gz"
# matrixOut <- "./exampleData/outputMatrix/chr22.tsv.gz"
# mutOut <- "./exampleData/outputMutations/chr22.tsv.gz"

# load
sites <- fread(targetIn, nThread = 1)
mutations_short <- fread(mutIn, nThread = 1)

# annotate mutations
mutations_short[
    sites,
    ':=' (
        geneID = i.geneID,
        strand = i.strand,
        bin	= i.bin,
        position = i.position,
        kmer = i.kmer,
        strand_kmer = i.strand_kmer
    ),
    on = c("chr", "start", "end")
]
mutations_short <- mutations_short[!is.na(kmer)]
mutations_short[, "strand_ALT" := ALT]
mutations_short[strand == "-", "strand_ALT" := as.character(reverseComplement(DNAStringSet(strand_ALT)))]

# mutational matrix
pkmers <- apply(
    gtools::permutations(4, 5, c("A", "C", "G", "T"), repeats.allowed = TRUE),
    1,
    paste,
    collapse = ""
)
mdt <- CJ(strand_kmer = pkmers, strand_ALT = c("A", "C", "G", "T"), unique = TRUE)
mdt[, "strand_REF" := substr(strand_kmer, 3, 3)]
mdt <- mdt[strand_REF != strand_ALT]
mdt[
    mutations_short[, list("n" = .N), by = c("strand_kmer", "strand_ALT")],
    "n" := i.n,
    on = c("strand_kmer", "strand_ALT")
]
mdt[is.na(n), "n" := 0L]
mdt[
    sites[, list("abundance" = .N), by = "strand_kmer"],
    "abundance" := i.abundance,
    on = "strand_kmer"
]
mdt[is.na(abundance), "abundance" := 0L]
mdt[, ':=' ("strand_REF" = NULL)]

# save
fwrite(mdt, matrixOut, sep = "\t", nThread = 1, compress = "gzip")
fwrite(mutations_short, mutOut, sep = "\t", compress = "gzip")

