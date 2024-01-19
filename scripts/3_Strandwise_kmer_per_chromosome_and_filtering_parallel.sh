#!/bin/sh
#$ -N Par_StrandKmer
#$ -cwd
#$ -j y
#$ -t 1-66
#$ -q short-centos79
#$ -l h_rt=02:00:00
#$ -l virtual_free=64G
#$ -o Cluster/3_Strandwise_kmers/$TASK_ID.out

source /users/dweghorn/dcastellano/anaconda3/bin/activate TSS

jobsfile=/users/dweghorn/cserranocolome/TSS/results/parallel_jobs_strandwise_kmers_all.txt
jobs=$(awk "NR==$SGE_TASK_ID" $jobsfile)
echo ${jobs}

chr=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $1}')
echo ${chr}
filtering=$(echo ${jobs} | sed 's/-/ /g' | awk '{print $2}')
echo ${filtering}

for donor_lists in UKBB others
	do
                if [ "$donor_lists" == "UKBB" ]; then
                        output_dir=BEDS       
                        echo ${output_dir}
                else
                        output_dir=BEDS_otherdatasets
                        echo ${output_dir}
                fi 
        Rscript --vanilla Germline_and_PCAWG_mutations_extraction.R ${chr} ${filtering} ${donor_lists}
        awk '{print $17, $18}' ../results/${output_dir}/${chr}_mutations_transcripts_and_co_FANTOM5_${filtering}_strand5mer.txt | sed 's/./& /g' | awk '{print $2, $3, $4, ">", $6}' | sed 's/ //g' > ../results/${output_dir}/${chr}_${filtering}_CpG.txt
        paste ../results/${output_dir}/${chr}_mutations_transcripts_and_co_FANTOM5_${filtering}_strand5mer.txt ../results/${output_dir}/${chr}_${filtering}_CpG.txt > ../results/${output_dir}/${chr}_mutations_transcripts_and_co_FANTOM5_${filtering}_strand5mer_CpG.txt
        rm ../results/BEDS_otherdatasets/${chr}_${filtering}_CpG.txt ../results/BEDS_otherdatasets/${chr}_mutations_transcripts_and_co_FANTOM5_${filtering}_strand5mer.txt ../results/BEDS_otherdatasets/${chr}_transcripts_and_co_FANTOM5_${filtering}_5mer.txt
	done

echo "The End"	