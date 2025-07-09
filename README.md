# TSS_hypermutability

Scripts to execute the main analysis presented in 10.1101/2025.02.04.635982
Check the "exampleData" directory for the format of the input files. The columns may be in any order but all of them are required. The "main.sh" script takes the input mutations and bin definitions for any region of interest and outputs a file with the number of mutations observed at each bin as well as the expected amount.

For usage refer to "example.sh". We provide the bin definitions we have used for the main analysis at https://zenodo.org/records/15846952 as well as the meta-cohorts with DNMs, early mosaic mutations and late mosaic mutations. The rest of the mutation dataset are available elsewhere from their original sources (gnomAD, UKBB, and PCAWG). 
