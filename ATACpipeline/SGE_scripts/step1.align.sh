#!/bin/bash
#$ -S /bin/sh
#$ -cwd
#$ -o results/out.step1.txt
#$ -e results/out.step2.txt
#$ -r y
#$ -pe orte 16
#$ -l mem_free=2G
#$ -l mem_token=2G
#$ -N BowTieJob
#$ -q all.q


# NUMCPU=16

# SAMPLEINFO_FN="/home/yshima/Data/nov_2014_data/sample_description.txt"
# SAMPLE_OPTION=""

# flowcells=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
# rawfastq1_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r1 $SAMPLE_OPTION` )
# rawfastq2_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r2 $SAMPLE_OPTION` )
# names=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )

# FLOWCELL=$flowcells[$SGE_TASK_ID]
# RAWFASTQ1_FN=$rawfastq1_fns[$SGE_TASK_ID]
# RAWFASTQ2_FN=$rawfastq2_fns[$SGE_TASK_ID]
# NAME=$names[$SGE_TASK_ID]
NAME="puma"

# PATHWAY PREPARATION (hard-coded)
BASEDIR="/home/paxorus/test.step1"
TRIM_FASTQF_DIR="$BASEDIR/trim_fastq"
TRIM_FASTQF1="$TRIM_FASTQF_DIR/Sample_EC_C_1.trim.fastq.gz"
TRIM_FASTQF2="$TRIM_FASTQF_DIR/Sample_EC_C_2.trim.fastq.gz"

TRACK_DIR="$BASEDIR/tracks"

SAMTOOLS_BIN="/home/yshima/programs/samtools-1.1/bin/samtools"
MARKDUPLICATES_BIN="java -jar -Xmx20g /home/yshima/programs/picard-tools-1.124/picard.jar MarkDuplicates"
BEDTOOLS_DIR="/home/yshima/programs/bedtools-2.17.0/bin/"
WIGTOBIGWIG_BIN="/home/yshima/programs/UCSCtools/wigToBigWig"
SCALEBED_BIN="/home/yshima/programs/ATACpipeline/misc/scale_bed_rpm.pl"
CHROMSIZES_FN="/home/yshima/references/mm10_chrom_and_size.txt"

BOWTIE2_BIN="/home/yshima/programs/bowtie2-2.2.4/bowtie2"
BOWTIE2_INDEX="/home/yshima/references/mm10_bt2index/mm10"
BOWTIE2_OPTIONS="-t -X2000 --no-mixed --no-discordant -p$NUMCPU -x $BOWTIE2_INDEX -1 $TRIM_FASTQF1 -2 $TRIM_FASTQF2 -S $SAM_FN"

# OUTPUT FILES
BAM_UNSORT_FN="$TRACK_DIR/$NAME.unsorted.bam"
BAM_FN_NOSUFFIX="$TRACK_DIR/$NAME"
BAM_FN="$TRACK_DIR/$NAME.bam"
BAM_DEDUP_FN="$TRACK_DIR/$NAME.dedup.bam"
MARKDUPLICATES_METRICS_FN="$TRACK_DIR/$NAME.MarkDuplicates.metrics.txt"
BED_FN="$TRACK_DIR/$NAME.quickviz.bed"
HISTO_BED_FN="$TRACK_DIR/$NAME-coverage.quickviz.bed"
HISTO_BIGWIG_FN="$TRACK_DIR/$NAME-coverage.quickviz.bw"
SCALED10M_BED_FN="$TRACK_DIR/$NAME_scaled10M.quickviz.bed"
SCALED10M_BIGWIG_FN="$TRACK_DIR/$NAME_scaled10M.quickviz.bw"

REGIONALSIGNAL_DIR="$BASEDIR/results/ATACseq_201411/regional_signals/"
FEATUREPOS_DIR="$BASEDIR/analysis/files/feature_bed"
TSS_DEDUP_COUNTS_FN="$REGIONALSIGNAL_DIR/$NAME-tss-counts-dedup.txt"
GENEBODIES_DEDUP_COUNTS_FN="$REGIONALSIGNAL_DIR/$NAME-genebodies-counts-dedup.txt"
COUNTS5KBWIN_DEDUP_FN="$REGIONALSIGNAL_DIR/$NAME-5kbwin_2.5kbstep_counts-dedup.txt"

SAM_DIR="$BASEDIR/results/ATACseq_201411/bowtie2_align"
SAM_FN="$SAM_DIR/$NAME.sam"




for T_OUTDIR in "$SAM_DIR/$TRACK_DIR/$REGIONALSIGNAL_DIR"
do
	if [ ! -e $T_OUTDIR]; then
		mkdir -p $T_OUTDIR
	fi
done

echo Program initialized.

curdir=$PWD
curhost=$HOSTNAME
echo "#sgejob run started on $curhost at ($(date))"

if [ ! -e $SAM_FN ]; then
	echo "# Aligning: $BOWTIE2_BIN $BOWTIE2_OPTIONS ($(date))"
	eval "$BOWTIE2_BIN $BOWTIE2_OPTIONS"
fi

if [ ! -e $BAM_FN ]; then
	echo "# Making BAM file"
	echo "#Converting SAM to BAM ($(date))"
	echo "$SAMTOOLS_BIN view -@ $NUMCPU -b -S $SAM_FN -o $BAM_UNSORT_FN"
	eval "$SAMTOOLS_BIN view -@ $NUMCPU -b -S $SAM_FN -o $BAM_UNSORT_FN"
	echo

	echo "#Sorting BAM ($(date))"
	echo "$SAMTOOLS_BIN sort -@ $NUMCPU -m 20G $BAM_UNSORT_FN $BAM_FN_NOSUFFIX"
	eval "$SAMTOOLS_BIN sort -@ $NUMCPU -m 5G $BAM_UNSORT_FN $BAM_FN_NOSUFFIX"

	echo "rm $BAM_UNSORT_FN"
	rm $BAM_UNSORT_FN
	echo

	curtime="$(date)"
	NUMALNS="$SAMTOOLS_BIN flagstat $BAM_FN | grep -m 1 mapped | awk '{print $1}'"
	echo "# Number of alignments: $NUMALNS"
fi

if [ ! -e $BAM_DEDUP_FN]; then
	curtime="$(date)"
	echo "# Removing duplicates ($curtime)"

	echo "$MARKDUPLICATES_BIN METRICS_FILE=$MARKDUPLICATES_METRICS_FN INPUT=$BAM_FN OUTPUT=$BAM_DEDUP_FN REMOVE_DUPLICATES=TRUE ($curtime)"
	eval "$MARKDUPLICATES_BIN METRICS_FILE=$MARKDUPLICATES_METRICS_FN INPUT=$BAM_FN OUTPUT=$BAM_DEDUP_FN REMOVE_DUPLICATES=TRUE"
fi

if [ ! -e $BED_FN]; then
	set curtime=`date`
	echo "#Converting to BED ($curtime)"
	echo "$BEDTOOLS_DIR/bamToBed -i $BAM_DEDUP_FN > $BED_FN"
	eval "$BEDTOOLS_DIR/bamToBed -i $BAM_DEDUP_FN > $BED_FN"
	echo
fi


if [ ! -e $HISTO_BED_FN]; then
	echo "#Calculating read coverage ($(date))"
	echo "$BEDTOOLS_DIR/genomeCoverageBed -g $CHROMSIZES_FN -i $BED_FN -bg > $HISTO_BED_FN"
	eval "$BEDTOOLS_DIR/genomeCoverageBed -g $CHROMSIZES_FN -i $BED_FN -bg > $HISTO_BED_FN"
	echo
fi

if [ ! -e $HISTO_BIGWIG_FN]; then
	echo "#Creating unnormalized BigWig track ($(date))"
	echo "$WIGTOBIGWIG_BIN $HISTO_BED_FN $CHROMSIZES_FN $HISTO_BIGWIG_FN"
	eval "$WIGTOBIGWIG_BIN $HISTO_BED_FN $CHROMSIZES_FN $HISTO_BIGWIG_FN"
	echo
fi

if [ ! -e $SCALED10M_BIGWIG_FN]; then
	echo "#Creating BigWig track of counts scaled to 10M alignments(=reads) ($(date))"

	NUMALNS="$SAMTOOLS_BIN flagstat $BAM_FN | grep -m 1 mapped | awk '{print $1}'"
	echo "# Number of alignments: $NUMALNS"

	eval "$BEDTOOLS_DIR/genomeCoverageBed -bg -ibam $BAM_FN -g $CHROMSIZES_FN | perl $SCALEBED_BIN $NUMALNS 10000000 > $SCALED10M_BED_FN"
	eval "$WIGTOBIGWIG_BIN $SCALED10M_BED_FN $CHROMSIZES_FN $SCALED10M_BIGWIG_FN"
fi

echo "#sgejob run finished on $curhost at $(date)"
