#!/bin/bash

cd /mnt/e/Cryo-TCR/data/TCR_data

for i in TCR_Raw_mixcr TCR_FilterCellranger_mixcr scRNA2_Raw_mixcr; do
	TCR_analysis_pipeline -i $i/ -o $i/myTCR_res/ -t ALL -m ./metadata.txt -c &> $i/log
	echo $i Over!
done

for i in TCR_Raw_cellranger; do
	TCR_analysis_pipeline -i $i/ -o $i/myTCR_res/ -t ALL -m ./metadata_cellranger.txt &> $i/log
	echo $i Over!
done

for i in Old_scRNA2_Raw_mixcr; do
	TCR_analysis_pipeline -i $i/ -o $i/myTCR_res/ -t ALL -m ./Old_metadata.txt -c &> $i/log
	echo $i Over!
done

