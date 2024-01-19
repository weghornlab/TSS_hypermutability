#!/bin/sh
#$ -N Par_Matrix
#$ -cwd
#$ -j y
#$ -t 1-220
#$ -q short-centos79
#$ -l h_rt=00:50:00
#$ -l virtual_free=64G
#$ -o Cluster/4_Mutation_matrix/$TASK_ID.out

source /users/dweghorn/dcastellano/anaconda3/bin/activate TSS

date

jobsfile=/users/dweghorn/cserranocolome/TSS/results/parallel_jobs_mutation_matrix_all.txt
jobs=$(awk "NR==$SGE_TASK_ID" $jobsfile)
echo ${jobs}

chr=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $1}')
echo ${chr}
donor_lists=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $2}')
echo ${donor_lists}
filtering=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $3}')
echo ${filtering}

Rscript --vanilla Mutation_matrix_per_donor_list.R ${chr} ${donor_lists} ${filtering}

echo "The End"
date