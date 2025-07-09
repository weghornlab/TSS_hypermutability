# all directories should exist beforehand, including the output directories
./main.sh \
    ./exampleData/inputMutations/ # input mutations
    ./exampleData/inputTargetDir/ # input target bins
    ./exampleData/outputMatrix/ # output directory for mutation matrix files
    ./exampleData/outputMutations/ # output directory for annotated mutations with bin data
    "no" # CpG2TpG consideration, can be "no", "yes" or "only"
    ./exampleData/outputDensityDir/ # main output directory
    1

# the main output will be in ./exampleData/outputDensityDir/all.tsv.gz
