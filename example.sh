# all directories should exist beforehand, including the output directories
# arguments in order:
# input mutations
# input target bins
# output directory for mutation matrix files
# output directory for annotated mutations with bin data
# CpG2TpG consideration, can be "no", "yes" or "only"
# main output directory
# number of cores

./main.sh ./exampleData/inputMutations/ ./exampleData/inputTarget/ ./exampleData/outputMatrix/ \
    ./exampleData/outputMutations/ "no" ./exampleData/outputDensity/ 1

# the main output will be in ./exampleData/outputDensityDir/all.tsv.gz
