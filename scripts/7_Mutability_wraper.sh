#!/bin/sh
#$ -N Mutabi_wraper
#$ -cwd
#$ -j y
#$ -q short-centos79
#$ -l h_rt=18000
#$ -l virtual_free=1G
#$ -o Cluster/7_Mutability_wraper_new_downsampled.out

		
for donor_lists in UKBB G_interval2 G_interval3 G_interval4 G_Common G_interval6 G_Trios G_ERV Pancancer G_All 
	do
		for filtering in CRG36 CRG36_maskExons CRG36_Exons
			do
				for bin_L in 1000 100
					do
						cat ../results/Expected_mutations_per_bin/Chr_parallel/*_${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_only_CpGs_${filtering}.txt > ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_only_CpGs_${filtering}.txt
						cat ../results/Expected_mutations_per_bin/Chr_parallel/*_${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_no_CpGs_${filtering}.txt > ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_no_CpGs_${filtering}.txt
						echo ${donor_lists} - ${filtering} - ${bin_L}
						awk '{print $1}' ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_no_CpGs_${filtering}.txt | sort | uniq | grep "chr" | wc -l
						awk '{print $1}' ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_only_CpGs_${filtering}.txt | sort | uniq | grep "chr" | wc -l
					done
			done
	done


# For downsampled ERVs
# for donor_lists in G_ERV 
# 	do
# 		for filtering in CRG36 #CRG36_maskExons CRG36_Exons all
# 			do
# 				for bin_L in 1000 #100 #10
# 					do
# 						cat ../results/Expected_mutations_per_bin/Chr_parallel_downsampled/*_${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_only_CpGs_${filtering}_downsampled.txt > ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_only_CpGs_${filtering}_downsampled.txt
# 						cat ../results/Expected_mutations_per_bin/Chr_parallel_downsampled/*_${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_no_CpGs_${filtering}_downsampled.txt > ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_no_CpGs_${filtering}_downsampled.txt
# 						echo ${donor_lists} - ${filtering} - ${bin_L}
# 						awk '{print $1}' ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_no_CpGs_${filtering}_downsampled.txt | sort | uniq | grep "chr" | wc -l
# 						awk '{print $1}' ../results/Expected_mutations_per_bin/${donor_lists}_obs_and_exp_mutations_${bin_L}_bin_size_only_CpGs_${filtering}_downsampled.txt | sort | uniq | grep "chr" | wc -l

# 					done
# 			done
# 	done

