#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step0.trim
#$ -e SGEOUT.atac.step0.trim
#$ -r y
#$ -pe orte 2
#$ -t 1-14 # number of samples

set NUMCPU=16


set SAMPLEINFO_FN="/home/yshima/Data/nov_2014_data/sample_description.txt"
set SAMPLE_OPTION=""

set flowcells=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
set rawfastq1_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r1 $SAMPLE_OPTION` )
set rawfastq2_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r2 $SAMPLE_OPTION` )
set names=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )



set FLOWCELL=$flowcells[$SGE_TASK_ID]
set RAWFASTQ1_FN=$rawfastq1_fns[$SGE_TASK_ID]
set RAWFASTQ2_FN=$rawfastq2_fns[$SGE_TASK_ID]
set NAME=$names[$SGE_TASK_ID]


set BASEDIR="/home/yshima/Data/nov_2014_data"
set FASTQF1="$BASEDIR/raw/${FLOWCELL}/$RAWFASTQ1_FN"
set FASTQF2="$BASEDIR/raw/${FLOWCELL}/$RAWFASTQ2_FN"

#set ADAPTOR_FASTF="$BASEDIR/data/external/adapter_sequences/transposon.fa"
set TRIM_FASTQF_DIR="$BASEDIR/results/ATACseq_201411/trim_fastq/{$FLOWCELL}"
set TRIM_FASTQF1="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r1_trim.fastq.gz"
set TRIM_FASTQF2="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r2_trim.fastq.gz"

set scratchdir="$BASEDIR/scratch/nodejobtmp-$$-$JOB_ID-$SGE_TASK_ID"
set TMP1_FASTQ="$scratchdir/tmp1.fastq.gz"
set TMP2_FASTQ="$scratchdir/tmp2.fastq.gz"

set CUTADAPT_BIN="/home/yshima/programs/cutadapt"
set CUTADAPT_OPTIONS1="-a CTGTCTCTTATACACATCT -q 30 --minimum-length 36 --paired-output ${TMP2_FASTQ} -o ${TMP1_FASTQ} ${FASTQF1} ${FASTQF2}"
set CUTADAPT_OPTIONS2="-a CTGTCTCTTATACACATCT -q 30 --minimum-length 36 --paired-output ${TRIM_FASTQF1} -o {$TRIM_FASTQF2} ${TMP2_FASTQ} ${TMP1_FASTQ}"

foreach T_OUTDIR ( $TRIM_FASTQF_DIR )	
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"

echo "#scratchdir = $scratchdir"
mkdir -p $scratchdir
cd $scratchdir

if (! -e $TRIM_FASTQF1 ) then
   set curtime=`date`
   echo "# Trimming FASTQ file - cutadapt step 1: ${CUTADAPT_BIN} ${CUTADAPT_OPTIONS1} ($curtime)"
   ${CUTADAPT_BIN} ${CUTADAPT_OPTIONS1}
   echo "# Trimming FASTQ file - cutadapt step 2: ${CUTADAPT_BIN} ${CUTADAPT_OPTIONS2} ($curtime)"
   ${CUTADAPT_BIN} ${CUTADAPT_OPTIONS2}

   rm ${TMP1_FASTQ} ${TMP2_FASTQ}
endif

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
