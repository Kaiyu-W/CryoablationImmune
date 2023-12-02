#!/bin/bash

cd /sibcb2/chenluonanlab7/wangkaiyu/202106/TCR
alias cellranger=/sibcb2/chenluonanlab7/wangkaiyu/cellranger-6.0.1/cellranger

# ##Renaming
# for i in *.fastq.gz; do
# 	index=$( echo $i | sed "s/\(.*CryoTCR_S1_L00\)\(.\)\(_.*\)/\2/" )
# 	type=$( echo $i | sed "s/\(.*\)\(_S1_L00.*\)/\1/" )
# 	suffix=$( echo $i | sed "s/\(.*_S1_L00._\)\(.*001.fastq.gz\)/\2/" )
# 	j=${type}_${index}_S1_L001_${suffix}
# 	mv $i $j
# done

##count processing

for i in *_S1_L001_I1_001.fastq.gz; do
	name=$( echo $i | sed "s/\(.*\)\(_S1_L001_I1_001.fastq.gz\)/\1/" )
	cellranger vdj --id=$name \
		--reference=/sibcb2/chenluonanlab7/wangkaiyu/refdata-cellranger-vdj-GRCm38-alts-ensembl-5.0.0 \
		--fastqs=./ \
		--sample=$name 
done

echo Finish!