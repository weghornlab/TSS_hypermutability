# TSS_hypermutability

Scripts from the TSS hypermutability project. The main pipelin to generate the bins, intersect the mutations and compute the mutabilities consists of the following 7 steps:

Step 1: 1_Seq_context.sh. Creates the 1 kbp and 100 bp bins for the analysis.

Step 2: 2_Mutations_preparation.sh. Generates files to compute in parallel the combination of different datasets, filters, bin sizes and chromosomes. 

Step 3: 3_Strandwise_kmer_per_chromosome_and_filtering_parallel.sh. Intersects the mutations with the different bins. 

Step 4: 4_Mutation_matrix_per_chromosome_and_donor_lists_parallel.sh. Computes the mutation matrices per each dataset and filter.

Step 5: 5_Mutation_matrix_wraper.sh. Step to put the parallel results together.

Step 6: 6_Mutability_per_transcript_chromosome_and_donor_lists_parallel.sh. Computes the mutability and number of mutations for each dataset in each bin.

Step 7: 7_Mutability_wraper.sh. Step to put the parallel results together.
