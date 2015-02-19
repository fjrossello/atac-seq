#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step1.align
#$ -e SGEOUT.atac.step1.align
#$ -r y
#$ -pe orte 16
#$ -t 1-14 # number of samples
#$ -l mem_free=2G
#$ -l mem_token=2G



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
set TRIM_FASTQF_DIR="$BASEDIR/results/ATACseq_201411/trim_fastq/{$FLOWCELL}"
set TRIM_FASTQF1="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r1_trim.fastq.gz"
set TRIM_FASTQF2="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r2_trim.fastq.gz"

set TRACK_DIR="$BASEDIR/results/ATACseq_201411/tracks/${FLOWCELL}"
set SAMTOOLS_BIN="/home/yshima/programs/samtools-1.1/bin/samtools"
set MARKDUPLICATES_BIN="java -jar -Xmx20g /home/yshima/programs/picard-tools-1.124/picard.jar MarkDuplicates"

set BEDTOOLS_DIR="/home/yshima/programs/bedtools-2.17.0/bin/"
set WIGTOBIGWIG_BIN="/home/yshima/programs/UCSCtools/wigToBigWig"
set SCALEBED_BIN="/home/yshima/programs/ATACpipeline/misc/scale_bed_rpm.pl"
set CHROMSIZES_FN="/home/yshima/references/mm10_chrom_and_size.txt"

## OUTPUT FILES
set BAM_UNSORT_FN="$TRACK_DIR/${NAME}.unsorted.bam"
set BAM_FN_NOSUFFIX="$TRACK_DIR/${NAME}"
set BAM_FN="$TRACK_DIR/${NAME}.bam"
set BAM_DEDUP_FN="$TRACK_DIR/${NAME}.dedup.bam"
set MARKDUPLICATES_METRICS_FN="$TRACK_DIR/${NAME}.MarkDuplicates.metrics.txt"
set BED_FN="$TRACK_DIR/${NAME}.quickviz.bed"
set HISTO_BED_FN="$TRACK_DIR/${NAME}-coverage.quickviz.bed"
set HISTO_BIGWIG_FN="$TRACK_DIR/${NAME}-coverage.quickviz.bw"
set SCALED10M_BED_FN="$TRACK_DIR/${NAME}_scaled10M.quickviz.bed"
set SCALED10M_BIGWIG_FN="$TRACK_DIR/${NAME}_scaled10M.quickviz.bw"

set REGIONALSIGNAL_DIR="$BASEDIR/results/ATACseq_201411/regional_signals/{$FLOWCELL}"
set FEATUREPOS_DIR="$BASEDIR/analysis/files/feature_bed"
set TSS_DEDUP_COUNTS_FN="$REGIONALSIGNAL_DIR/{$NAME}-tss-counts-dedup.txt"
set GENEBODIES_DEDUP_COUNTS_FN="$REGIONALSIGNAL_DIR/{$NAME}-genebodies-counts-dedup.txt"
set COUNTS5KBWIN_DEDUP_FN="$REGIONALSIGNAL_DIR/{$NAME}-5kbwin_2.5kbstep_counts-deup.txt"


set SAM_DIR="$BASEDIR/results/ATACseq_201411/bowtie2_align/${FLOWCELL}"
set SAM_FN="$SAM_DIR/${FLOWCELL}_${NAME}.sam"
set BOWTIE2_BIN="/home/yshima/programs/bowtie2-2.2.4/bowtie2"
set BOWTIE2_INDEX="/home/yshima/references/mm10_bt2index/mm10"
set BOWTIE2_OPTIONS="-t -X2000 --no-mixed --no-discordant -p$NUMCPU -x $BOWTIE2_INDEX -1 $TRIM_FASTQF1 -2 $TRIM_FASTQF2 -S $SAM_FN"

foreach T_OUTDIR ( $SAM_DIR $TRACK_DIR $REGIONALSIGNAL_DIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end


set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"

if (! -e $SAM_FN) then
   set curtime=`date`
   echo "# Aligning: $BOWTIE2_BIN $BOWTIE2_OPTIONS ($curtime)"
   $BOWTIE2_BIN $BOWTIE2_OPTIONS
endif

if (! -e $BAM_FN) then
   echo "# Making BAM file"
   set curtime=`date`
   echo "#Converting SAM to BAM ($curtime)"
   echo "$SAMTOOLS_BIN view -@ $NUMCPU -b -S $SAM_FN -o $BAM_UNSORT_FN"
   $SAMTOOLS_BIN view -@ $NUMCPU -b -S $SAM_FN -o $BAM_UNSORT_FN
   echo

   set curtime=`date`
   echo "#Sorting BAM ($curtime)"
   echo "$SAMTOOLS_BIN sort -@ $NUMCPU -m 20G $BAM_UNSORT_FN $BAM_FN_NOSUFFIX"
   $SAMTOOLS_BIN sort -@ $NUMCPU -m 5G $BAM_UNSORT_FN $BAM_FN_NOSUFFIX

   echo "rm $BAM_UNSORT_FN"
   rm $BAM_UNSORT_FN
   echo

   set curtime=`date`
   set NUMALNS=`$SAMTOOLS_BIN flagstat $BAM_FN | grep -m 1 mapped | awk '{print $1}'`
   echo "# Number of alignments: $NUMALNS"
endif

if (! -e $BAM_DEDUP_FN) then
   set curtime=`date`
   echo "# Removing duplicates ($curtime)"

   echo "$MARKDUPLICATES_BIN METRICS_FILE=$MARKDUPLICATES_METRICS_FN INPUT=$BAM_FN OUTPUT=$BAM_DEDUP_FN REMOVE_DUPLICATES=TRUE ($curtime)"
   $MARKDUPLICATES_BIN METRICS_FILE=$MARKDUPLICATES_METRICS_FN INPUT=$BAM_FN OUTPUT=$BAM_DEDUP_FN REMOVE_DUPLICATES=TRUE
endif

if (! -e $BED_FN) then
   set curtime=`date`
   echo "#Converting to BED ($curtime)"
   echo "$BEDTOOLS_DIR/bamToBed -i $BAM_DEDUP_FN > $BED_FN"
   $BEDTOOLS_DIR/bamToBed -i $BAM_DEDUP_FN > $BED_FN
   echo
endif


if (! -e $HISTO_BED_FN) then
   set curtime=`date`
   echo "#Calculating read coverage ($curtime)"
   echo "$BEDTOOLS_DIR/genomeCoverageBed -g $CHROMSIZES_FN -i $BED_FN -bg > $HISTO_BED_FN"
   $BEDTOOLS_DIR/genomeCoverageBed -g $CHROMSIZES_FN -i $BED_FN -bg > $HISTO_BED_FN
   echo
endif

if (! -e $HISTO_BIGWIG_FN) then
   set curtime=`date`
   echo "#Creating unnormalized BigWig track ($curtime)"
   echo "$WIGTOBIGWIG_BIN $HISTO_BED_FN $CHROMSIZES_FN $HISTO_BIGWIG_FN"
   $WIGTOBIGWIG_BIN $HISTO_BED_FN $CHROMSIZES_FN $HISTO_BIGWIG_FN
   echo
endif

if (! -e $SCALED10M_BIGWIG_FN) then
   set curtime=`date`
   echo "#Creating BigWig track of counts scaled to 10M alignments(=reads) ($curtime)"

   set NUMALNS=`$SAMTOOLS_BIN flagstat $BAM_FN | grep -m 1 mapped | awk '{print $1}'`
   echo "# Number of alignments: $NUMALNS"

   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -ibam $BAM_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 10000000 > $SCALED10M_BED_FN
   $WIGTOBIGWIG_BIN $SCALED10M_BED_FN $CHROMSIZES_FN $SCALED10M_BIGWIG_FN
endif

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
