#!/bin/bash

# # all raw data process by mixcr
# for i in 1 2 3 4; do
# 	analysis_name=CryoTCR_${i}
# 	input_file=CryoTCR_S1_L00${i}
# 	input_file1=${input_file}_R1_001.fastq.gz
# 	input_file2=${input_file}_R2_001.fastq.gz
# 	mixcr analyze amplicon \
# 	    -s mmu \
# 	    --starting-material rna \
# 	    --5-end no-v-primers --3-end c-primers \
# 	    --adapters adapters-present \
# 	    $input_file1 $input_file2 $analysis_name
# 	echo $analysis_name OVER!
# done

# for i in 1 2 3 4; do
# 	analysis_name=NonCryoTCR_${i}
# 	input_file=NonCryoTCR_S1_L00${i}
# 	input_file1=${input_file}_R1_001.fastq.gz
# 	input_file2=${input_file}_R2_001.fastq.gz
# 	mixcr analyze amplicon \
# 	    -s mmu \
# 	    --starting-material rna \
# 	    --5-end no-v-primers --3-end c-primers \
# 	    --adapters adapters-present \
# 	    $input_file1 $input_file2 $analysis_name
# 	echo $analysis_name OVER!
# done

# # filtered and merged data processed by cellranger vdj, then process by mixcr
# for i in CryoTCR NonCryoTCR; do
# 	input_file=/mnt/e/Cryo-TCR/data/TCR_Cryo_Cellranger/${i}/outs/filtered_contig.fastq
# 	mixcr analyze amplicon \
# 	    -s mmu \
# 	    --starting-material rna \
# 	    --5-end no-v-primers --3-end c-primers \
# 	    --adapters no-adapters \
# 	    $input_file $i
# 	echo $i OVER!
# done

# all filter data process by mixcr
for i in 1 2 3 4; do
	analysis_name=CryoTCR_${i}
	input_file=CryoTCR_${i}.filter.fastq
	mixcr analyze amplicon \
	    -s mmu \
	    --starting-material rna \
	    --5-end no-v-primers --3-end c-primers \
	    --adapters adapters-present \
	    $input_file ./res/$analysis_name
	echo $analysis_name OVER!
done

for i in 1 2 3 4; do
	analysis_name=NonCryoTCR_${i}
	input_file=NonCryoTCR_${i}.filter.fastq
	mixcr analyze amplicon \
	    -s mmu \
	    --starting-material rna \
	    --5-end no-v-primers --3-end c-primers \
	    --adapters adapters-present \
	    $input_file ./res/$analysis_name
	echo $analysis_name OVER!
done