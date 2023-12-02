#!/bin/bash
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
vdjtools=/mnt/c/vdjtools-1.2.1/vdjtools-1.2.1.jar
vdjtools_patch=/mnt/c/vdjtools-1.2.1/vdjtools-patch.sh

if [ $# == 0 ]; then
	mode=raw
elif [ $# == 1 ]; then
	mode=$1
else
	echo "Wrong parameter. Only need one 'mode'"
  	exit 1
fi

# 1. MiXCR process
# not shown

cd /mnt/e/Cryo-TCR/data/TCR_data/res
mkdir vdjtools

# 2. VDJTools process and plot
# format convert
java -jar $vdjtools Convert -S mixcr ./*.ALL.txt ./vdjtools/

# process and plot
cd vdjtools

java -jar $vdjtools CalcBasicStats *.clonotypes.ALL.txt .
java -jar $vdjtools CalcPairwiseDistances -p --low-mem *.clonotypes.ALL.txt .
java -jar $vdjtools JoinSamples -p *.clonotypes.ALL.txt . # only the first 5 sample at most
java -jar $vdjtools RarefactionPlot *.clonotypes.ALL.txt .

if [[ $mode == filter ]]; then
	for j in Cryo NonCryo; do
		java -jar $vdjtools PlotQuantileStats -t 5 ${j}TCR.clonotypes.ALL.txt ./${j}TCR
		java -jar $vdjtools PlotFancySpectratype -t 20 ${j}TCR.clonotypes.ALL.txt ./${j}TCR
		java -jar $vdjtools PlotSpectratypeV -u ${j}TCR.clonotypes.ALL.txt ./${j}TCR
		java -jar $vdjtools PlotSpectratypeV ${j}TCR.clonotypes.ALL.txt ./${j}TCR
		$vdjtools_patch java -jar $vdjtools PlotFancyVJUsage ${j}TCR.clonotypes.ALL.txt ./${j}TCR
	done
else
	for i in 1 2 3 4; do
		for j in Cryo NonCryo; do
			java -jar $vdjtools PlotQuantileStats -t 5 ${j}TCR_${i}.clonotypes.ALL.txt ./${j}TCR_${i}
			java -jar $vdjtools PlotFancySpectratype -t 20 ${j}TCR_${i}.clonotypes.ALL.txt ./${j}TCR_${i}
			java -jar $vdjtools PlotSpectratypeV -u ${j}TCR_${i}.clonotypes.ALL.txt ./${j}TCR_${i}
			java -jar $vdjtools PlotSpectratypeV ${j}TCR_${i}.clonotypes.ALL.txt ./${j}TCR_${i}
			$vdjtools_patch java -jar $vdjtools PlotFancyVJUsage ${j}TCR_${i}.clonotypes.ALL.txt ./${j}TCR_${i}
		done
	done
fi

java -jar $vdjtools ClusterSamples -p -l . .
java -jar $vdjtools CalcSegmentUsage -p *.clonotypes.ALL.txt .

# 3. R-DESeq2 analysis and plot
# Rscript /mnt/e/Cryo-TCR/data/TCR_data/diff_analysis.r V $(pwd) $(pwd)
# Rscript /mnt/e/Cryo-TCR/data/TCR_data/diff_analysis.r J $(pwd) $(pwd)
# R code running in jupyter script, since input data of R comes from python output.

# 4. Python plot
cp /mnt/e/Cryo-TCR/data/TCR_data/stats.ipynb .
runipy stats.ipynb


echo ALL OVER!
