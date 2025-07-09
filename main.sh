#!/usr/bin/env bash
inMutsDir="${1}"
inTargetDir="${2}"
outMatrixDir="${3}"
outMutsDir="${4}"
include_CpGs="${5}"
outDensityDir="${6}"
ncores="${7}"

# debug
inMutsDir=./exampleData/inputMutations/
inTargetDir=./exampleData/inputTarget/
outMatrixDir=./exampleData/outputMatrix/
outMutsDir=./exampleData/outputMutations/
include_CpGs="no"
outDensityDir=./exampleData/outputDensity/
ncores=1

# show arguments
echo "inMutsDir: ${inMutsDir}"
echo "inTargetDir: ${inTargetDir}"
echo "outMatrixDir: ${outMatrixDir}"
echo "include_CpGs: ${include_CpGs}"
echo "inDensityDir: ${outDensityDir}"
echo "ncores: ${ncores}"

# get chromosomes with mutations
chrs=$(find "${inMutsDir}" -type f -exec basename {} \;)

# caculate mutational matrix and annotate mutations, both by chromosome
echo "MUTATIONAL MATRIX..."
inMutsFiles=($(echo "${chrs}" | sed 's,^,'"${inMutsDir}"','))
inTargetFiles=($(echo "${chrs}" | sed 's,^,'"${inTargetDir}"','))
outMatrixFiles=($(echo "${chrs}" | sed 's,^,'"${outMatrixDir}"','))
outMutsFiles=($(echo "${chrs}" | sed 's,^,'"${outMutsDir}"','))
parallel --link --progress -j"${ncores}" Rscript mutmat.r ::: \
    "${inMutsFiles[@]}" ::: \
    "${inTargetFiles[@]}" ::: \
    "${outMatrixFiles[@]}" ::: \
    "${outMutsFiles[@]}"

# join mutational matrices across chromosomes
find "${outMatrixDir}" -type f -name "chr*.tsv.gz" -exec sh -c "zcat {} | head -n1" \; | uniq > "${outMatrixDir}"all.tsv
find "${outMatrixDir}" -type f -name "chr*.tsv.gz" -exec sh -c "zcat {} | tail --quiet -n +2" \; >> "${outMatrixDir}"all.tsv
gzip -f "${outMatrixDir}"all.tsv

# caculate mutational densities by window and gene, by chromosome
echo "MUTATION DENSITY..."
outDensityFiles=($(echo "${chrs}" | sed 's,^,'"${outDensityDir}"','))
parallel --link --progress -j"${ncores}" Rscript density.r ::: \
    "${outMutsFiles[@]}" ::: \
    "${inTargetFiles[@]}" ::: \
    "${outMatrixDir}"all.tsv.gz ::: \
    "${include_CpGs}" ::: \
    "${outDensityFiles[@]}"

# join mutational matrices across chromosomes
find "${outDensityFiles}" -type f -name "chr*.tsv.gz" -exec sh -c "zcat {} | head -n1" \; | uniq > "${outDensityDir}"all.tsv
find "${outDensityFiles}" -type f -name "chr*.tsv.gz" -exec sh -c "zcat {} | tail --quiet -n +2" \; >> "${outDensityDir}"all.tsv
gzip -f "${outDensityDir}"all.tsv

echo "DONE..."
