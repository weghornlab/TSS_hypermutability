#!/bin/sh
#$ -N Matrix_wraper
#$ -cwd
#$ -j y
#$ -q short-centos79
#$ -l h_rt=1800
#$ -l virtual_free=1G
#$ -o Cluster/5_Mutation_matrix.out

		
for donor_lists in UKBB G_ERV G_Trios G_interval2 G_interval3 G_interval4 G_Common G_interval6 G_All Pancancer
	do
		for filtering in CRG36_maskExons CRG36_Exons CRG36 
			do
				rm ../results/Mutation_Matrices/${donor_lists}_${filtering}.txt
				cat ../results/Mutation_Matrices/*_${donor_lists}_${filtering}.txt > ../results/Mutation_Matrices/${donor_lists}_${filtering}.txt
				rm ../results/Mutation_Matrices/*_${donor_lists}_${filtering}.txt
				echo ${donor_lists} - ${filtering}
				awk '{print $5}' ../results/Mutation_Matrices/${donor_lists}_${filtering}.txt | sort | uniq | grep "chr" | wc -l
			done
	done
	