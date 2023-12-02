# ##Renaming
# mkdir scRNA_Cryo
cd LBFC20210771/210611_A00679_0575_AH527KDSX2/
# for i in Cryo*; do
# 	index_raw=$( echo $i | sed "s/\(.*_S\)\(.\)\(_.*\)/\2/" )
# 	index=$[$index_raw-1]
# 	type=$( echo $i | sed "s/\(.*L004_\)\(..\)\(_001.*\)/\2/" )
# 	j=Cryo_S1_L00${index}_${type}_001.fastq.gz
# 	cp $i ../../scRNA_Cryo/${j}
# done
# for i in non*; do
# 	index_raw=$( echo $i | sed "s/\(.*_S\)\(.\)\(_.*\)/\2/" )
# 	index=$[index_raw-5]
# 	type=$( echo $i | sed "s/\(.*L004_\)\(..\)\(_001.*\)/\2/" )
# 	j=NonCryo_S1_L00${index}_${type}_001.fastq.gz
# 	cp $i ../../scRNA_Cryo/${j}
# done

##count processing
cd ../../scRNA_Cryo/
cellranger count --id=Cryo \
	--transcriptome=/mnt/c/Cellranger_Reference/refdata-gex-mm10-2020-A \
	--fastqs=./ \
	--sample=Cryo \
	--nosecondary
cellranger count --id=NonCryo \
	--transcriptome=/mnt/c/Cellranger_Reference/refdata-gex-mm10-2020-A \
	--fastqs=./ \
	--sample=NonCryo \
	--nosecondary

echo Finish!
