#!/bin/bash

#config
#Rawdata_path="/mnt/d/Lab"
#in_file1="ZA-1_R1.fastq.gz"
#in_file2="ZA-1_R2.fastq.gz"

# if [ $# -ne 1 ]; then
#   echo    "usage:$0 <Rawdata_file_1.fastq(Absolute_Path)>"
#   echo    "e.g: $0 /mnt/d/Lab/ZA_01/ZA-1_R1.fastq.gz  "
#   exit 1
# fi

workingdirectory=/sibcb1/chenluonanlab6/wangkaiyu/Ji_lab/bulk
index_fa=/sibcb1/chenluonanlab6/wangkaiyu/musculus/Mus_musculus.GRCm38.dna.primary_assembly.fa
index_gtf=/sibcb1/chenluonanlab6/wangkaiyu/musculus/Mus_musculus.GRCm38.101.gtf
index_transcriptome=/sibcb1/chenluonanlab6/wangkaiyu/musculus/Mus_musculus.GRCm38.cdna.all.fa
symbollist_gene=/sibcb1/chenluonanlab6/wangkaiyu/musculus/geneid2symbol.tsv
symbollist_transcript=/sibcb1/chenluonanlab6/wangkaiyu/musculus/transcriptid2symbol.tsv


cd $workingdirectory

# Fastqc
mkdir -p ${workingdirectory}/process/fastqc_output
fastqc -o ${workingdirectory}/process/fastqc_output -t 16 ./data/*.fastq.gz

echo "QC OVER! QC output saved in ${workingdirectory}/process/fastqc_output"


# STAR
mkdir -p ${workingdirectory}/process/star_index
mkdir -p ${workingdirectory}/process/star_output
mkdir -p ${workingdirectory}/process/samtools_output_STAR
mkdir -p ${workingdirectory}/process/res_counts

# build index:
index_dir=${workingdirectory}/process/star_index
if [ ! -d $index_dir ]; then
    STAR --runMode genomeGenerate \
         --genomeDir $index_dir \
         --runThreadN 40 \
         --genomeFastaFiles $index_fa \
         --sjdbGTFfile $index_gtf \
         --limitGenomeGenerateRAM 100000000000
fi

echo STAR index saved in $index_dir


# align and samtools process
cd ${workingdirectory}/data
name=$( ls *R1.fastq.gz | sort -u | sed "s/_R1.fastq.gz//" )
for i in $name; do

    STAR \
        --genomeDir ${workingdirectory}/process/star_index \
        --runThreadN 20 \
        --readFilesIn ./${i}_R1.fastq.gz ./${i}_R2.fastq.gz \
        --readFilesCommand zcat \
        --outFileNamePrefix ${workingdirectory}/process/star_output/${i}_ \
        --outSAMtype BAM SortedByCoordinate \
        --outBAMsortingThreadN 10 \
        --outSAMstrandField intronMotif \
        --quantMode TranscriptomeSAM GeneCounts # transcriptome index, no
        # --bamRemoveDuplicatesType
        
    # samtools
    cd ${workingdirectory}/process/samtools_output_STAR

    samtools view -b -q 20 -o ${i}_q20.bam ${workingdirectory}/process/star_output/${i}_Aligned.toTranscriptome.out.bam
    samtools fixmate -m ${i}_q20.bam ${i}_q20_fixmate.bam
    samtools sort ${i}_q20_fixmate.bam --output-fmt BAM -o ${i}_q20_sorted.bam
    samtools markdup -r ${i}_q20_sorted.bam ${i}_q20_markdup.bam
    samtools index ${i}_q20_markdup.bam

    samtools flagstat ${i}_q20_markdup.bam --threads 8 > ${i}_flagstat.txt
    samtools idxstats ${i}_q20_markdup.bam > ${i}_idxstats.txt

    # copy count results
    # cp ${i}_idxstats.txt ${workingdirectory}/process/res_counts
    # cp ${workingdirectory}/process/star_output/${i}_ReadsPerGene.out.tab ${workingdirectory}/process/res_counts
    gawk 'BEGIN{print "Gene""\t""Length""\t""Mapped Counts""\t""Unmapped Counts"}{print $0}' ${i}_idxstats.txt > ${workingdirectory}/process/res_counts/${i}_idxstats.tsv
    sed "1,4 d" ${workingdirectory}/process/star_output/${i}_ReadsPerGene.out.tab | gawk 'BEGIN{print "Gene""\t""Unstrand""\t""First strand""\t""Second strand"}{print $0}' > ${workingdirectory}/process/res_counts/${i}_ReadsPerGene.tsv
    
    cd ${workingdirectory}/data
done

echo STAR alignment and samtools process OVER!


# merge all results
cd ${workingdirectory}/process/res_counts

Rscript ${workingdirectory}/merge_idxstats.r ${workingdirectory}/process/res_counts
Rscript ${workingdirectory}/merge_ReadsPerGene.r ${workingdirectory}/process/res_counts
sed -n '/^\*/ !p' idxstats_merge.tsv > temp.tsv
mv temp.tsv idxstats_merge.tsv

echo Merge results saved in ${workingdirectory}/process/res_counts


# statistics of mapping quality
cd ${workingdirectory}/process/star_output
echo SAMPLE INPUT_READs UNIQUE_MAPPING MULTI_MAPPING UNMAPPED_toomanymismatches UNMAPPED_tooshort UNMAPPED_other CHIMERIC > ${workingdirectory}/process/res_counts/res_stat.txt
for i in *_Log.final.out; do 
    SAMPLE=$( echo $i | sed "s/_Log.final.out//" )
    INPUT_READs=$( sed -n "/Number of input reads/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )
    UNIQUE_MAPPING=$( sed -n "/Uniquely mapped reads %/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )
    MULTI_MAPPING=$( sed -n "/% of reads mapped to multiple loci/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )
    UNMAPPED_toomanymismatches=$( sed -n "/% of reads unmapped: too many mismatches/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )
    UNMAPPED_tooshort=$( sed -n "/% of reads unmapped: too short/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )
    UNMAPPED_other=$( sed -n "/% of reads unmapped: other/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )
    CHIMERIC=$( sed -n "/% of chimeric reads/p" $i | gawk -F "|\t" '{print $2}' | sed "s/ //g" )

    echo $SAMPLE $INPUT_READs $UNIQUE_MAPPING $MULTI_MAPPING $UNMAPPED_toomanymismatches $UNMAPPED_tooshort $UNMAPPED_other $CHIMERIC >> ${workingdirectory}/process/res_counts/res_stat.txt
done

echo Statistics results saved in ${workingdirectory}/process/res_counts/res_stat.txt


# change Ensembl id to gene symbol

cd ${workingdirectory}/process/res_counts
Rscript ${workingdirectory}/Chid_idxstats.r ${workingdirectory}/process/res_counts $symbollist_transcript 
Rscript ${workingdirectory}/Chid_ReadsPerGene.r ${workingdirectory}/process/res_counts $symbollist_gene 


echo ALL process over!

