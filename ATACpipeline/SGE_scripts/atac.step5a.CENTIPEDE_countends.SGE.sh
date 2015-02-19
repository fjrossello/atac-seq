#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step5a.CENTIPEDE_countends
#$ -e SGEOUT.atac.step5a.CENTIPEDE_countends
#$ -r y
#$ -pe batch 16 
#$ -t 1-6

set NUMCPU=16

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

## INPUT
set BED_FN="$VIZ_DIR/${NAME}.bed"

## OUTPUT
set ALL_ENDS_BIGWIG_FN="$VIZ_DIR/${NAME}.rawcounts.ends.bw"
set ALL_ENDS_REALSTRAND_FW_BIGWIG_FN="$VIZ_DIR/${NAME}.rawcounts.fw_ends.bw"
set ALL_ENDS_REALSTRAND_RV_BIGWIG_FN="$VIZ_DIR/${NAME}.rawcounts.rv_ends.bw"



foreach T_OUTDIR ( $TRACK_DIR $VIZ_DIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"

## SET UP scratch directory
set scratchdir="/scratch/davisf/nodejobtmp-$$-$JOB_ID-$SGE_TASK_ID"
echo "#scratchdir = $scratchdir"
mkdir -p $scratchdir


## Calc bw track of RAW counts of end positions
if (! -e $ALL_ENDS_BIGWIG_FN || ! -e $ALL_ENDS_REALSTRAND_FW_BIGWIG_FN) then

   set TEMP_BED_FN="$scratchdir/temp.bed"
   set TEMP2_BED_FN="$scratchdir/temp2.bed"
   set TEMP3_BED_FN="$scratchdir/temp3.bed"

   echo "# Calculating ends read coverage ($curtime)"

   set curtime=`date`
   awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,($2 + 1),".",0,"+"; print $1,($3 - 1),$3,".",0,"-"}' < $BED_FN > $TEMP_BED_FN
   sort -k1,1 -k2,2n -k3,3n $TEMP_BED_FN > $TEMP2_BED_FN
   rm $TEMP_BED_FN


   if (! -s $ALL_ENDS_REALSTRAND_FW_BIGWIG_FN) then
      set curtime=`date`
      echo "#-> Create fw coverage track ($curtime)"
      echo "${BEDTOOLS_DIR}/genomeCoverageBed -strand + -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} > $TEMP3_BED_FN ($curtime)"
      ${BEDTOOLS_DIR}/genomeCoverageBed -strand + -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} > $TEMP3_BED_FN

      set curtime=`date`
      echo "#-> Convert fw coverage to bigWig ($curtime)"
      echo "$WIGTOBIGWIG_BIN $TEMP3_BED_FN $CHROMSIZES_FN $ALL_ENDS_REALSTRAND_FW_BIGWIG_FN ($curtime)"
      $WIGTOBIGWIG_BIN $TEMP3_BED_FN $CHROMSIZES_FN $ALL_ENDS_REALSTRAND_FW_BIGWIG_FN

      rm $TEMP3_BED_FN

      set curtime=`date`
      echo "#-> Create rv coverage track ($curtime)"
      echo "${BEDTOOLS_DIR}/genomeCoverageBed -strand - -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} > $TEMP3_BED_FN ($curtime)"
      ${BEDTOOLS_DIR}/genomeCoverageBed -strand - -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} > $TEMP3_BED_FN

      echo "#-> Convert rv strands to bigWig ($curtime)"
      echo "$WIGTOBIGWIG_BIN $TEMP3_BED_FN $CHROMSIZES_FN $ALL_ENDS_REALSTRAND_RV_BIGWIG_FN ($curtime)"
      $WIGTOBIGWIG_BIN $TEMP3_BED_FN $CHROMSIZES_FN $ALL_ENDS_REALSTRAND_RV_BIGWIG_FN

      rm $TEMP3_BED_FN
   endif


   if (! -s $ALL_ENDS_BIGWIG_FN ) then
      set curtime=`date`
      echo "#-> Creating BigWig track of raw end counts ($curtime)"
      echo "${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} > $TEMP_BED_FN ($curtime)"
      ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} > $TEMP_BED_FN

      set curtime=`date`
      echo "$WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $ALL_ENDS_BIGWIG_FN ($curtime)"
      $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $ALL_ENDS_BIGWIG_FN
   endif

   rm $TEMP_BED_FN $TEMP2_BED_FN
endif

rmdir $scratchdir

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
