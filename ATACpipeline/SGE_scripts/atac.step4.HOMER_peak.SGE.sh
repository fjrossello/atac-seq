#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step4.HOMER_peaks
#$ -e SGEOUT.atac.step4.HOMER_peaks
#$ -r y
#$ -pe batch 2 
#$ -t 1-6

set NUMCPU=2


set SAMPLEINFO_FN="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact/data/amo_samples_ATAC.txt"
set SAMPLE_OPTION="dataset run201409-merge"

set flowcells=( `perl txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
set names=( `perl txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )



set FLOWCELL=$flowcells[$SGE_TASK_ID]
set NAME=$names[$SGE_TASK_ID]


set BASEDIR="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact"

set TRACK_DIR="$BASEDIR/results/ATACseq.20140915/tracks/${FLOWCELL}"
set VIZ_DIR="$BASEDIR/results/ATACseq.20140915/viz/${FLOWCELL}"
set SAMTOOLS_BIN="/groups/eddy/home/davisf/bin/samtools-0.1.19/samtools"
set BEDTOOLS_DIR="/groups/eddy/home/davisf/bin/bedtools-v2.15.0"
set WIGTOBIGWIG_BIN="/groups/eddy/home/davisf/software/bigWig/wigToBigWig"
set SCALEBED_BIN="/groups/eddy/home/davisf/work/genomics/fpdgenlib/fpdgenlib_work/trunk/src/scripts/scale_bed_rpm.pl"
set CHROMSIZES_FN="$BASEDIR/data/external/rsem/mm10_genome_size.tab"

set HOMER_BINDIR="/groups/eddy/home/davisf/software/HOMER/bin"

## OUTPUT FILES
set BED_FN="$VIZ_DIR/${NAME}.bed"
set BED_SUB100_FN="$VIZ_DIR/${NAME}.sub100nt.bed"
set BED_SUP150_FN="$VIZ_DIR/${NAME}.sup150nt.bed"

set HOMER_BASE_DIR="$BASEDIR/results/ATACseq.20140915/HOMER/${FLOWCELL}/${NAME}"
set ATAC_PEAKS_FN="$HOMER_BASE_DIR/${NAME}.HOMER_peaks.bed"

set HOMER_TAG_DIR="$HOMER_BASE_DIR/tags"
set HOMER_PEAKS_FN="$HOMER_TAG_DIR/peaks.txt"
set HOMER_PEAKS_BED_FN="$HOMER_TAG_DIR/peaks.bed"

foreach T_OUTDIR ( $HOMER_TAG_DIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"

if (! -e $HOMER_PEAKS_FN) then
   set curtime=`date`
   echo "# ${HOMER_BINDIR}/makeTagDirectory $HOMER_TAG_DIR -forceBED $BED_SUB100_FN  ($curtime)"
   ${HOMER_BINDIR}/makeTagDirectory $HOMER_TAG_DIR -forceBED $BED_SUB100_FN

   set curtime=`date`
   echo "# ${HOMER_BINDIR}/findPeaks $HOMER_TAG_DIR -region -size 500 -minDist 50 -o auto -tbp 0 ($curtime)"
   ${HOMER_BINDIR}/findPeaks $HOMER_TAG_DIR -region -size 500 -minDist 50 -o auto -tbp 0
endif

if (! -e $HOMER_PEAKS_BED_FN) then
   set curtime=`date`
   echo "# perl ${HOMER_BINDIR}/pos2bed.pl $HOMER_PEAKS_FN > $HOMER_PEAKS_BED_FN ($curtime)"
   perl ${HOMER_BINDIR}/pos2bed.pl $HOMER_PEAKS_FN > $HOMER_PEAKS_BED_FN
endif

if (! -e $ATAC_PEAKS_FN) then
   set curtime=`date`
   echo "${BEDTOOLS_DIR}/bedtools sort -i $HOMER_PEAKS_BED_FN | ${BEDTOOLS_DIR}/bedtools merge -i > $ATAC_PEAKS_FN ($curtime)"
   ${BEDTOOLS_DIR}/bedtools sort -i $HOMER_PEAKS_BED_FN | ${BEDTOOLS_DIR}/bedtools merge -i > $ATAC_PEAKS_FN
endif

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
