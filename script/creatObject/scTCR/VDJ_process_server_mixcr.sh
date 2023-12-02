#!/bash/bin

export PATH=/sibcb1/chenluonanlab6/wangkaiyu/MiXCR/mixcr-3.0.13:$PATH
mkdir /sibcb1/chenluonanlab6/wangkaiyu/Ji_lab/202106/TCR_mixcr
cd /sibcb1/chenluonanlab6/wangkaiyu/Ji_lab/202106/TCR_mixcr

for i in 1 2 3 4; do
	analysis_name=CryoTCR_${i}
	input_file1=/sibcb2/chenluonanlab7/wangkaiyu/202106/Cryo_S1_L00${i}_R1_001.fastq.gz
	input_file2=/sibcb2/chenluonanlab7/wangkaiyu/202106/Cryo_S1_L00${i}_R2_001.fastq.gz
	mixcr analyze shotgun \
	    -s mmu \
	    --starting-material rna \
	    $input_file1 $input_file2 ./$analysis_name
	echo $analysis_name OVER!
done

for i in 1 2 3 4; do
	analysis_name=NonCryoTCR_${i}
	input_file1=/sibcb2/chenluonanlab7/wangkaiyu/202106/NonCryo_S1_L00${i}_R1_001.fastq.gz
	input_file2=/sibcb2/chenluonanlab7/wangkaiyu/202106/NonCryo_S1_L00${i}_R2_001.fastq.gz
	mixcr analyze shotgun \
	    -s mmu \
	    --starting-material rna \
	    $input_file1 $input_file2 ./$analysis_name
	echo $analysis_name OVER!
done