#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -o results/out.step0.txt
#$ -e results/error.step0.txt
#$ -r y
#$ -pe orte 2
#$ -N TrimJob
#$ -q all.q

NUMCPU=16

# set SAMPLEINFO_FN="/home/yshima/Data/nov_2014_data/sample_description.txt"
# set SAMPLE_OPTION=""

# set flowcells=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
# set rawfastq1_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r1 $SAMPLE_OPTION` )
# set rawfastq2_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r2 $SAMPLE_OPTION` )
# set names=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )

# set FLOWCELL=$flowcells[$SGE_TASK_ID]
# set RAWFASTQ1_FN=$rawfastq1_fns[$SGE_TASK_ID]
# set RAWFASTQ2_FN=$rawfastq2_fns[$SGE_TASK_ID]
# set NAME=$names[$SGE_TASK_ID]

BASEDIR="/home/paxorus/test"
FASTQF1="$BASEDIR/raw/Sample_EC_C_1.fastq.gz"
FASTQF2="$BASEDIR/raw/Sample_EC_C_2.fastq.gz"

TRIM_FASTQF_DIR="$BASEDIR/results"
TRIM_FASTQF1="$TRIM_FASTQF_DIR/Sample_EC_C_1.trim.fastq.gz"
TRIM_FASTQF2="$TRIM_FASTQF_DIR/Sample_EC_C_2.trim.fastq.gz"

scratchdir="$BASEDIR/scratch/nodejobtmp-$$-$JOB_ID" # "-$SGE_TASK_ID"
TMP1_FASTQ="$scratchdir/tmp1.fastq.gz"
TMP2_FASTQ="$scratchdir/tmp2.fastq.gz"

CUTADAPT_BIN="/home/yshima/anaconda/bin/cutadapt"
CUTADAPT_OPTIONS1="-a CTGTCTCTTATACACATCT -q 30 --minimum-length 36 --paired-output $TMP2_FASTQ -o $TMP1_FASTQ $FASTQF1 $FASTQF2"
CUTADAPT_OPTIONS2="-a CTGTCTCTTATACACATCT -q 30 --minimum-length 36 --paired-output $TRIM_FASTQF1 -o $TRIM_FASTQF2 $TMP2_FASTQ $TMP1_FASTQ"

echo Program initialized.

curdir=$PWD
curhost=$HOSTNAME
curtime="$(date)"
echo "#sgejob run started on $curhost at $(date)"

echo "#scratchdir = $scratchdir"
mkdir -p $scratchdir
cd $scratchdir

# run cutadapt with options
if [ ! -e $TRIM_FASTQF1 ]; then
    set curtime=`date`
	# fastqf -> tmp_fastq
    echo "# Trimming FASTQ file - cutadapt step 1: $CUTADAPT_BIN $CUTADAPT_OPTIONS1 ($curtime)"
	eval "$CUTADAPT_BIN $CUTADAPT_OPTIONS1"
	# tmp_fastq -> trim_fastqf
    echo "# Trimming FASTQ file - cutadapt step 2: $CUTADAPT_BIN $CUTADAPT_OPTIONS2 ($curtime)"
    eval "$CUTADAPT_BIN $CUTADAPT_OPTIONS2"

	if [ -e $TMP1_FASTQ ]; then
		rm $TMP1_FASTQ
	fi
	if [ -e $TMP2_FASTQ ]; then
		rm $TMP2_FASTQ
	fi
fi

echo "#sgejob run finished on $curhost at $(date)"
