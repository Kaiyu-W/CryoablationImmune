# ##Renaming
mkdir TCR_Cryo_denovo
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
cd ../../TCR_Cryo_denovo/
cellranger vdj --id=CryoTCR \
	--fastqs=../TCR_Cryo/ \
	--denovo \
	--inner-enrichment-primers=../TCR_inner_reverse_primers.fa \
	--sample=CryoTCR 
cellranger vdj --id=NonCryoTCR \
	--denovo \
	--inner-enrichment-primers=../TCR_inner_reverse_primers.fa \
	--fastqs=../TCR_Cryo/ \
	--sample=NonCryoTCR 

echo Finish!
