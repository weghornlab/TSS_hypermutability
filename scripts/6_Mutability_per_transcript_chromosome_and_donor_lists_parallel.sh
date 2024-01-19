#!/bin/sh
#$ -N Par_Mutability
#$ -cwd
#$ -j y
#$ -t 1-1
#$ -q short-centos79
#$ -l h_rt=04:30:00
#$ -l virtual_free=64G
#$ -o Cluster/6_Mutability_transcript/$TASK_ID.out

source /users/dweghorn/dcastellano/anaconda3/bin/activate TSS
date

jobsfile=/users/dweghorn/cserranocolome/TSS/results/parallel_jobs_bin_mutability_all.txt
jobs=$(awk "NR==$SGE_TASK_ID" $jobsfile)
echo ${jobs}

# chr=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $1}')
# echo ${chr}
CpG=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $2}')
echo ${CpG}
window_length=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $3}')
echo ${window_length}
donor_lists=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $4}')
echo ${donor_lists}
filtering=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $5}')
echo ${filtering}


for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22
    do
        Rscript --vanilla Bin_mutability_parallel_ERVs_downsampling_adapted.R ${chr} ${CpG} ${window_length} ${donor_lists} ${filtering}
        Rscript --vanilla Bin_mutability_parallel.R ${chr} ${CpG} ${window_length} ${donor_lists} ${filtering}
    done

echo "The End"
date
exit





