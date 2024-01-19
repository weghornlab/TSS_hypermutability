#!/bin/sh
#$ -N Extracting_transcripts_context
#$ -cwd
#$ -j y
#$ -q short-centos79,long-centos79
#$ -l h_rt=432000
#$ -l virtual_free=64G
#$ -o Cluster/1_Seq_context.out

source /users/dweghorn/dcastellano/anaconda3/bin/activate TSS

### Transcripts with mappability filters
## Exons and GERP elements to remove
zcat /users/dweghorn/dcastellano/muCOV/data/Exons/gencode.v19.annotation.gff3 | grep "exon" | awk '{print $1, $4, $5, "exon"}' | sed 's/ /\t/g' | uniq > ../results/exon_and_gerp_elements.bed 
awk '{print $1, $2, $3, "GERP"}' /users/dweghorn/dcastellano/muCOV/Conservation/GERP_sum_bins_of_interest.bed | sed 's/ /\t/g' | uniq >> ../results/exon_and_gerp_elements.bed 
sort -V -k1,1 -k2,2 ../results/exon_and_gerp_elements.bed | uniq > ../results/exon_and_gerp_elements_sorted.bed 
bedtools merge -i ../results/exon_and_gerp_elements_sorted.bed > ../results/exon_and_gerp_elements_collapsed.bed

zcat /users/dweghorn/dcastellano/muCOV/data/Exons/gencode.v19.annotation.gff3 | grep "exon" | awk '{print $1, $4, $5, "exon"}' | sed 's/ /\t/g' | uniq > ../results/exon_elements.bed 
sort -V -k1,1 -k2,2 ../results/exon_elements.bed | uniq > ../results/exon_elements_sorted.bed 
bedtools merge -i ../results/exon_elements_sorted.bed > ../results/exon_elements_collapsed.bed

for chr in chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 #
	do
		echo START ${chr}
		
		echo Extracting sites
			Rscript --vanilla Extracting_downstream_transcript_seq.R ${chr}
			# Rscript --vanilla Transcribed_sites_distribution.R ${chr} ## printing sites transcribed by more than one transcript
			Rscript --vanilla Extracting_upstream_transcript_seq.R ${chr}
			Rscript --vanilla Extracting_postTTS_transcript_seq.R ${chr}
			Rscript --vanilla Extracting_upstream_TTS_seq.R ${chr}

		echo Removing duplicated and transcribed sites
			Rscript --vanilla Rmv_transcribed_sites_postTTS.R ${chr}
			Rscript --vanilla Rmv_transcribed_sites_upstream.R ${chr}
			Rscript --vanilla Rmv_transcribed_sites_downstream.R ${chr}
			Rscript --vanilla Rmv_transcribed_sites_upstream_TTS.R ${chr}
			Rscript --vanilla Rmv_transcribed_sites.R ${chr}

		echo Removing CRG36
			bedtools intersect -a ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_all.bed -b /users/dweghorn/dcastellano/muCOV/data/Mappability_filter/wgEncodeCrgMapabilityAlign36mer_score1.bed > ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_CRG36.bed
			wait

		echo + Removing exons and GERP elements
			bedtools intersect -v -a ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_CRG36.bed -b ../results/exon_and_gerp_elements_collapsed.bed > ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_CRG36_maskExons.bed
			wait

		echo + Keeping only exons and GERP elements
			bedtools intersect -a ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_CRG36.bed -b ../results/exon_and_gerp_elements_collapsed.bed > ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_CRG36_Exons.bed
			wait	

		for file in all CRG36 CRG36_maskExons CRG36_Exons
			do
				echo Extracting the kmer
					awk '{print $1, $2-2, $3+2}' ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}.bed  | sed 's/ /\t/g' | awk '(NR>1) && ($2 >= 0 ) ' > ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}_5mer.bed
					bedtools getfasta -tab -fi /users/dweghorn/dcastellano/muCOV/data/Fastas/GRCh37.primary_assembly.genome.fa -bed ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}_5mer.bed -fo ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}_5mer.bed.fa
					awk '{print $2-2, $0}' ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}.bed  | sed 's/ /\t/g' | awk '(NR>1) && ($1 >= 0 ) ' > ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}.bed2
					paste ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}.bed2 ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}_5mer.bed.fa | awk '{print $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12, $14}' > ../results/BEDS/${chr}_transcripts_and_co_FANTOM5_${file}_5mer.txt
				
				echo Computing unique bins
					Rscript --vanilla Unique_bins.R ${chr} ${file}
			done

		rm ../results/BEDS/${chr}_*.bed ../results/BEDS/${chr}_*.fa ../results/BEDS/${chr}_*.bed2
		echo END ${chr}
	done

