#!/bin/sh
#$ -N Mutations_input
#$ -cwd
#$ -j y
#$ -q short-centos79,long-centos79
#$ -l h_rt=1000
#$ -l virtual_free=1G
#$ -o Cluster/2_Mutations_preparation.out


rm ../results/parallel_jobs_bin_mutability.txt ../results/parallel_jobs_mutation_matrix.txt ../results/parallel_jobs_strandwise_kmers.txt


## Mutations extraction and kmer strand oriented ###
for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 #
	do
		for filtering in all CRG36 CRG36_maskExons CRG36_Exons
			do
				echo ${chr}-${filtering} >> ../results/parallel_jobs_strandwise_kmers.txt
			done	
	done


### List of parameters to run in parallel ###
for donor_lists in UKBB G_ERV G_interval2 G_interval3 G_interval4 G_Common G_interval6 G_Trios Pancancer G_All 
	do
		for filtering in CRG36_maskExons CRG36_Exons all CRG36 
			do
				for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
					do
						echo ${chr}-${donor_lists}-${filtering} >> ../results/parallel_jobs_mutation_matrix_all.txt
						
						for window_length in 1000 100
							do
								for CpG in no only
									do
										echo ${chr}-${CpG}-${window_length}-${donor_lists}-${filtering} >> ../results/parallel_jobs_bin_mutability_all.txt
									done
							done		
					done
			done		
	done		
