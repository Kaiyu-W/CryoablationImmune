# ##Renaming
# mkdir TCR_Cryo
cd LBFC20210771/210616_A00268_0667_BH5225DSX2/
# for i in Cryo*; do
# 	index_raw=$( echo $i | sed "s/\(.*_S\)\(..\)\(_.*\)/\2/" )
# 	index=$[$index_raw-15]
# 	type=$( echo $i | sed "s/\(.*L002_\)\(..\)\(_001.*\)/\2/" )
# 	j=CryoTCR_S1_L00${index}_${type}_001.fastq.gz
# 	cp $i ../../TCR_Cryo/${j}
# done
# for i in non*; do
# 	index_raw=$( echo $i | sed "s/\(.*_S\)\(..\)\(_.*\)/\2/" )
# 	index=$[index_raw-19]
# 	type=$( echo $i | sed "s/\(.*L002_\)\(..\)\(_001.*\)/\2/" )
# 	j=NonCryoTCR_S1_L00${index}_${type}_001.fastq.gz
# 	cp $i ../../TCR_Cryo/${j}
# done

##count processing
cd ../../TCR_Cryo/
cellranger vdj --id=CryoTCR \
	--reference=/mnt/c/Cellranger_Reference/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
	--fastqs=./ \
	--sample=CryoTCR 
cellranger vdj --id=NonCryoTCR \
	--reference=/mnt/c/Cellranger_Reference/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
	--fastqs=./ \
	--sample=NonCryoTCR 

echo Finish!
