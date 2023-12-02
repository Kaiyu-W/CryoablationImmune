#!/bin/bash

# configure tools
vdjtools=/mnt/c/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjtools_patch=/mnt/c/vdjtools-1.2.1/vdjtools-patch.sh

# configure parameters received from shell
n_para=0
mode=raw
Convert=FALSE
while getopts 'hfci:o:t:m:' OPT; do
	case ${OPT} in
		f)
			mode=filter
			echo -e "Run script with data of merged and filtered pattern \n"
        	;;
       	c)
			Convert=TRUE
			echo -e "First transfer the input file to vdjtools format \n"
			;;
      	i)
        	data_dir=${OPTARG}
        	n_para=$[$n_para + 1]
        	;;
      	o)
	        res_dir=${OPTARG}
	        n_para=$[$n_para + 1]
	        ;;
      	t)
	        VDJ_type=${OPTARG}
	        n_para=$[$n_para + 1]
	        ;;
      	m)
			metadata=${OPTARG}
			n_para=$[$n_para + 1]
			;;
      	h)
       		echo -e "Options:\n\tshell.sh -f <use -f if filtered(2, merged from cellranger or others), or not if 4> \n\t\t -c <use -c if transferring format of input to what vdjtools needs> \n\t\t -i <data_dir> \n\t\t -o <res_dir> \n\t\t -t <VDJ_type(ALL/TRA/...)> \n\t\t -m <metadata_absolute_path> \n\t\t -h\n"
       		exit 0
        	;;
      	*)
			echo -e "Unrecognized parameter. Please check it carefully!\n"
	        exit 1
	        ;;
    esac
done
if [ $n_para -ne 4 ]; then
	echo "Received not enough parameters!"
	echo -e "Options:\n\tshell.sh -f <use -f if filtered(2, merged from cellranger or others), or not if 4> \n\t\t -c <use -c if transferring format of input to what vdjtools needs> \n\t\t -i <data_dir> \n\t\t -o <res_dir> \n\t\t -t <VDJ_type(ALL/TRA/...)> \n\t\t -m <metadata_absolute_path> \n\t\t -h\n"
	exit 1
fi
if [ ! -d $data_dir ]; then
	echo -e "Data_dir not exist!\n"
	exit 1
fi
if [[ "ALL IGH IGK IGL TRA TRB TRD TRG" =~ "$VDJ_type" ]]; then
	Nothing=0
else
	echo -e "VDJ_type error! Only choose from 'ALL IGH IGK IGL TRA TRB TRD TRG'!\n"
	exit 1
fi
if [ ! -f $metadata ]; then
	echo -e "metadata not exist!\n"
	exit 1
fi
if [ $(echo $metadata | grep -c "^/") -ne 1 ]; then
	metadata=$(pwd)/$metadata
fi
if [ ! -d $res_dir ]; then
	mkdir -p $res_dir
fi

#####################
# >>> conda initialize >>>
# !! Contents within this block are managed by 'conda init' !!
__conda_setup="$('/home/wky/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/home/wky/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/home/wky/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/home/wky/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<
#####################

conda activate vdjtools

# process
# 1. MiXCR process
# not shown

# cd $data_dir
# mkdir vdjtools

# 2. VDJTools process and plot
# format convert
if [[ $Convert == TRUE ]]; then
	java -jar $vdjtools Convert -S mixcr $data_dir/*.clonotypes.${VDJ_type}.txt $res_dir
else
	cp $data_dir/*.clonotypes.${VDJ_type}.txt $res_dir
fi

if [ ! -f $res_dir/*.clonotypes.${VDJ_type}.txt ]; then
	echo -e "Data didn't exist! \n\tFiles should be as format of {SAMPLE}TCR.clonotypes.{VDJ_type}.txt\n"
	exit 1
fi

# process and plot
cd $res_dir

java -jar $vdjtools CalcBasicStats *.clonotypes.${VDJ_type}.txt .
java -jar $vdjtools CalcPairwiseDistances -p --low-mem *.clonotypes.${VDJ_type}.txt .
java -jar $vdjtools JoinSamples -p *.clonotypes.${VDJ_type}.txt . # only the first 5 sample at most
java -jar $vdjtools RarefactionPlot *.clonotypes.${VDJ_type}.txt .

if [[ $mode == filter ]]; then
	for j in Cryo NonCryo; do
		java -jar $vdjtools PlotQuantileStats -t 5 ${j}TCR.clonotypes.${VDJ_type}.txt ./${j}TCR
		java -jar $vdjtools PlotFancySpectratype -t 20 ${j}TCR.clonotypes.${VDJ_type}.txt ./${j}TCR
		java -jar $vdjtools PlotSpectratypeV -u ${j}TCR.clonotypes.${VDJ_type}.txt ./${j}TCR
		java -jar $vdjtools PlotSpectratypeV ${j}TCR.clonotypes.${VDJ_type}.txt ./${j}TCR
		$vdjtools_patch java -jar $vdjtools PlotFancyVJUsage ${j}TCR.clonotypes.${VDJ_type}.txt ./${j}TCR
	done
else
	for i in 1 2 3 4; do
		for j in Cryo NonCryo; do
			java -jar $vdjtools PlotQuantileStats -t 5 ${j}TCR_${i}.clonotypes.${VDJ_type}.txt ./${j}TCR_${i}
			java -jar $vdjtools PlotFancySpectratype -t 20 ${j}TCR_${i}.clonotypes.${VDJ_type}.txt ./${j}TCR_${i}
			java -jar $vdjtools PlotSpectratypeV -u ${j}TCR_${i}.clonotypes.${VDJ_type}.txt ./${j}TCR_${i}
			java -jar $vdjtools PlotSpectratypeV ${j}TCR_${i}.clonotypes.${VDJ_type}.txt ./${j}TCR_${i}
			$vdjtools_patch java -jar $vdjtools PlotFancyVJUsage ${j}TCR_${i}.clonotypes.${VDJ_type}.txt ./${j}TCR_${i}
		done
	done
fi

java -jar $vdjtools ClusterSamples -p -l . .
java -jar $vdjtools CalcSegmentUsage -p *.clonotypes.${VDJ_type}.txt .

# 3. R-DESeq2 analysis and plot
# Rscript /mnt/e/Cryo-TCR/data/TCR_data/diff_analysis.r V $(pwd) $(pwd)
# Rscript /mnt/e/Cryo-TCR/data/TCR_data/diff_analysis.r J $(pwd) $(pwd)
# R code running in jupyter script, since input data of R comes from python output.

# 4. Python plot

python /mnt/e/Cryo-TCR/data/TCR_data/stats.py -i ./ -o ./res_python/ -m $metadata -s .clonotypes.${VDJ_type}.txt

echo ${VDJ_type} OVER!
