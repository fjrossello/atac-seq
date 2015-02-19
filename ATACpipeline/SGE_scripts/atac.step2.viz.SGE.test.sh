#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step2.viz
#$ -e SGEOUT.atac.step2.viz
#$ -r y
#$ -pe orte 16 


set NUMCPU=16

set SAMPLEINFO_FN="/home/yshima/Data/nov_2014_data/sample_description.txt"
set SAMPLE_OPTION=""


set flowcells=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
set rawfastq1_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r1 $SAMPLE_OPTION` )
set rawfastq2_fns=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r2 $SAMPLE_OPTION` )
set names=( `perl /home/yshima/programs/ATACpipeline/SGE_scripts/txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )

set FLOWCELL=$flowcells[1]
set RAWFASTQ1_FN=$rawfastq1_fns[1]
set RAWFASTQ2_FN=$rawfastq2_fns[1]
set NAME=$names[1]


set BASEDIR="/home/yshima/Data/nov_2014_data"
set TRIM_FASTQF_DIR="$BASEDIR/results/ATACseq_201411/trim_fastq/{$FLOWCELL}"
set TRIM_FASTQF1="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r1_trim.fastq.gz"
set TRIM_FASTQF2="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r2_trim.fastq.gz"

set TRACK_DIR="$BASEDIR/results/ATACseq_201411/tracks/${FLOWCELL}"
set VIZ_DIR="$BASEDIR/results/ATACseq_201411/viz/${FLOWCELL}"
set SAMTOOLS_BIN="/home/yshima/programs/samtools-1.1/bin/samtools"
set BEDTOOLS_DIR="/home/yshima/programs/bedtools-2.17.0/bin/"
set WIGTOBIGWIG_BIN="/home/yshima/programs/UCSCtools/wigToBigWig"
set SCALEBED_BIN="/home/yshima/programs/ATACpipeline/misc/scale_bed_rpm.pl"
set CHROMSIZES_FN="/home/yshima/references/mm10_chrom_and_size.txt"

## OUTPUT FILES
set BAM_UNSORT_FN="$TRACK_DIR/${NAME}.unsorted.bam"
set BAM_FN_NOSUFFIX="$TRACK_DIR/${NAME}"
set BAM_FN="$TRACK_DIR/${NAME}.bam"

set BAM_DEDUP_FN="$TRACK_DIR/${NAME}.dedup.bam"

set BAM_DEDUP_PROPPAIR_FN="$TRACK_DIR/${NAME}.dedup_proppair_unsorted.bam"
set BAM_DEDUP_PROPPAIR_SORT_NOSUFFIX="$TRACK_DIR/${NAME}.dedup_proppair"
set BAM_DEDUP_PROPPAIR_SORT_FN="$TRACK_DIR/${NAME}.dedup_proppair.bam"

set BED_UNSORT_FN="$VIZ_DIR/${NAME}_unsort.bed"
set BED_FN="$VIZ_DIR/${NAME}.bed"
set BED_FN_BASE="${NAME}.bed"

set BED_SUB100_FN="$VIZ_DIR/${NAME}.sub100nt.bed"
set BED_SUP150_FN="$VIZ_DIR/${NAME}.sup150nt.bed"

set SCALED10M_ALL_BIGWIG_FN="$VIZ_DIR/${NAME}_scaled10M.bw"
set SCALED10M_SUB100_BIGWIG_FN="$VIZ_DIR/${NAME}_sub100_scaled10M.bw"
set SCALED10M_SUP150_BIGWIG_FN="$VIZ_DIR/${NAME}_sup150_scaled10M.bw"

set SCALED10M_ALL_ENDS_BIGWIG_FN="$VIZ_DIR/${NAME}_scaled10M.ends.bw"
set SCALED10M_SUB100_ENDS_BIGWIG_FN="$VIZ_DIR/${NAME}_sub100_scaled10M.ends.bw"
set SCALED10M_SUP150_ENDS_BIGWIG_FN="$VIZ_DIR/${NAME}_sup150_scaled10M.ends.bw"


set REGIONALSIGNAL_DIR="$BASEDIR/results/ATACseq_201411/regional_signals/{$FLOWCELL}"
set FEATUREPOS_DIR="$BASEDIR/analysis/files/feature_bed"

set TSS_COUNTS_SUB100NT_FN="$REGIONALSIGNAL_DIR/{$NAME}-tss-counts-sub100.txt"
set GENEBODIES_COUNTS_SUB100NT_FN="$REGIONALSIGNAL_DIR/{$NAME}-genebodies-counts-sub100.txt"
##set COUNTS5KBWIN_SUB100NT_FN="$REGIONALSIGNAL_DIR/{$NAME}-5kbwin_2.5kbstep-counts-sub100.txt"

set TSS_COUNTS_SUP150NT_FN="$REGIONALSIGNAL_DIR/{$NAME}-tss-counts-sup150.txt"
set GENEBODIES_COUNTS_SUP150NT_FN="$REGIONALSIGNAL_DIR/{$NAME}-genebodies-counts-sup150.txt"
##set COUNTS5KBWIN_SUP150NT_FN="$REGIONALSIGNAL_DIR/{$NAME}-5kbwin_2.5kbstep-counts-sup150.txt"

set RSCRIPT_BIN="/home/yshima/programs/R-3.1.2/bin/Rscript "
set RSCRIPT_LENDISTR_FN="/home/yshima/programs/ATACpipeline/SGE_scripts/atac.step2.viz.R"
set FRAGLEN_TMP_FN="$VIZ_DIR/${NAME}_fraglen.txt"
set FRAGLEN_PDF_FN="$VIZ_DIR/${NAME}_fraglen.pdf"
set FRAGLEN_LOG_PDF_FN="$VIZ_DIR/${NAME}_fraglen_log.pdf"


set TWOBIT_FN="/groups/eddy/home/davisf/work/databases/ucsc_genomes/mm10/fasta/mm10.2bit"

foreach T_OUTDIR ( $TRACK_DIR $VIZ_DIR $REGIONALSIGNAL_DIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"

## SET UP scratch directory
set scratchdir="/scratch/davisf/nodejobtmp-$$-$JOB_ID-1"
echo "#scratchdir = $scratchdir"
mkdir -p $scratchdir

# Step 1. Sort properly-paired alignments
if (! -e $BAM_DEDUP_PROPPAIR_SORT_FN) then
   set curtime=`date`
   echo "#Sorting properly-paired alignments ($curtime)"
   $SAMTOOLS_BIN view -@ $NUMCPU -bf 0x2 $BAM_DEDUP_FN -o $BAM_DEDUP_PROPPAIR_FN
   $SAMTOOLS_BIN sort -n -@ $NUMCPU -m 5G $BAM_DEDUP_PROPPAIR_FN $BAM_DEDUP_PROPPAIR_SORT_NOSUFFIX
   rm $BAM_DEDUP_PROPPAIR_FN
endif

# Step 2. Make BED files of all, sub-100, and super-150nt fragments
if (! -e $BED_FN) then
   set curtime=`date`
   echo "# Converting BAM to BED ($curtime)"

## fpd 140915_2232 -- changed from +4/-5 to +4/-4
   $BEDTOOLS_DIR/bamToBed -bedpe -i $BAM_DEDUP_PROPPAIR_SORT_FN | cut -f 1,2,6 | awk 'BEGIN{OFS="\t"} {print $1, ($2 + 4), ($3 - 4)}' > $BED_UNSORT_FN
   sort -k1,1 -k2,2n -k3,3n $BED_UNSORT_FN > $BED_FN
   awk '{if (($3 - $2) < 100) print}' $BED_FN > $BED_SUB100_FN
   awk '{if (($3 - $2) >= 150) print}' $BED_FN > $BED_SUP150_FN
   echo
endif


# Step 3. Convert to bw
if (! -e $SCALED10M_ALL_BIGWIG_FN) then
   set TEMP_BED_FN="$scratchdir/temp.bed"

   echo "# Calculating read coverage ($curtime)"

   set curtime=`date`
   echo "#-> Counting number of alignments ($curtime)"
   set NUMALNS=`wc -l $BED_FN | awk '{print $1}'`
   echo "#-> NUMALNS=$NUMALNS"

   set curtime=`date`
   echo "#-> Creating BigWig track of counts scaled to 10M alignments(=reads) ($curtime)"
   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $BED_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 10000000 > $TEMP_BED_FN
   $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $SCALED10M_ALL_BIGWIG_FN

   rm $TEMP_BED_FN
endif

## Step 3.5. Calc bw track of end positions (rather than full read coverage ala ChIP or RNAseq)
if (! -e $SCALED10M_ALL_ENDS_BIGWIG_FN) then

   set TEMP_BED_FN="$scratchdir/temp.bed"
   set TEMP2_BED_FN="$scratchdir/temp2.bed"

   echo "# Calculating ends read coverage ($curtime)"

   set curtime=`date`
   echo "#-> Counting number of alignments ($curtime)"
   set NUMALNS=`wc -l $BED_FN | awk '{print $1}'`
   echo "#-> NUMALNS=$NUMALNS"

   echo "#-> Making ends BED file ($curtime)"
   awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,($2 + 1),".",0,"+"; print $1,($3 - 1),$3,".",0,"-"}' < $BED_FN > $TEMP_BED_FN
   sort -k1,1 -k2,2n -k3,3n $TEMP_BED_FN > $TEMP2_BED_FN
   rm $TEMP_BED_FN

   set curtime=`date`
   echo "#-> Creating BigWig track of end counts scaled to 10M alignments(=reads) ($curtime)"
   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 5000000 > $TEMP_BED_FN
   $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $SCALED10M_ALL_ENDS_BIGWIG_FN
## NOTE: normalizing to 5M, because the denominator is number of fragments, while numerator is number of ends

   rm $TEMP_BED_FN $TEMP2_BED_FN
endif


# Step 4. Convert sub-100 to bw
if (! -e $SCALED10M_SUB100_BIGWIG_FN) then
   set TEMP_BED_FN="$scratchdir/temp.bed"

   set curtime=`date`
   echo "# Calculating read coverage OF SUB-100nt fragments ($curtime)"

   set curtime=`date`
   echo "#-> Counting number of alignments ($curtime)"
   set NUMALNS=`wc -l $BED_SUB100_FN | awk '{print $1}'`
   echo "#-> NUMALNS=$NUMALNS"

   set curtime=`date`
   echo "#-> Creating BigWig track of counts scaled to 10M alignments(=reads) ($curtime)"
   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $BED_SUB100_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 10000000 > $TEMP_BED_FN
   $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $SCALED10M_SUB100_BIGWIG_FN

   rm $TEMP_BED_FN
endif

## Step 4.5. Calc bw track of end positions (rather than full read coverage ala ChIP or RNAseq)
if (! -e $SCALED10M_SUB100_ENDS_BIGWIG_FN) then
   set TEMP_BED_FN="$scratchdir/temp.bed"
   set TEMP2_BED_FN="$scratchdir/temp2.bed"

   set curtime=`date`
   echo "# Calculating ends read coverage OF SUB-100nt fragments ($curtime)"

   set curtime=`date`
   echo "#-> Counting number of alignments ($curtime)"
   set NUMALNS=`wc -l $BED_SUB100_FN | awk '{print $1}'`
   echo "#-> NUMALNS=$NUMALNS"

   echo "#-> Making ends BED file ($curtime)"
   awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,($2 + 1),".",0,"+"; print $1,($3 - 1),$3,".",0,"-"}' < $BED_SUB100_FN > $TEMP_BED_FN
   sort -k1,1 -k2,2n -k3,3n $TEMP_BED_FN > $TEMP2_BED_FN
   rm $TEMP_BED_FN

   set curtime=`date`
   echo "#-> Creating BigWig track of end counts scaled to 10M alignments(=reads) ($curtime)"
   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 5000000 > $TEMP_BED_FN
   $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $SCALED10M_SUB100_ENDS_BIGWIG_FN


   rm $TEMP_BED_FN $TEMP2_BED_FN
endif

# Step 5. Convert sup-150 to bw
if (! -e $SCALED10M_SUP150_BIGWIG_FN) then
   set TEMP_BED_FN="$scratchdir/temp.bed"

   set curtime=`date`
   echo "# Calculating read coverage OF SUPER-150nt fragments ($curtime)"

   set curtime=`date`
   echo "#-> Counting number of alignments ($curtime)"
   set NUMALNS=`wc -l $BED_SUP150_FN | awk '{print $1}'`
   echo "#-> NUMALNS=$NUMALNS"

   set curtime=`date`
   echo "#-> Creating BigWig track of counts scaled to 10M alignments(=reads) ($curtime)"
   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $BED_SUP150_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 10000000 > $TEMP_BED_FN
   $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $SCALED10M_SUP150_BIGWIG_FN

   rm $TEMP_BED_FN
endif

## Step 5.5. Calc bw track of end positions (rather than full read coverage ala ChIP or RNAseq)
if (! -e $SCALED10M_SUP150_ENDS_BIGWIG_FN) then
   set TEMP_BED_FN="$scratchdir/temp.bed"
   set TEMP2_BED_FN="$scratchdir/temp2.bed"

   set curtime=`date`
   echo "# Calculating ends read coverage OF SUPER-150nt fragments ($curtime)"

   set curtime=`date`
   echo "#-> Counting number of alignments ($curtime)"
   set NUMALNS=`wc -l $BED_SUP150_FN | awk '{print $1}'`
   echo "#-> NUMALNS=$NUMALNS"

   echo "#-> Making ends BED file ($curtime)"
   awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,($2 + 1),".",0,"+"; print $1,($3 - 1),$3,".",0,"-"}' < $BED_SUP150_FN > $TEMP_BED_FN
   sort -k1,1 -k2,2n -k3,3n $TEMP_BED_FN > $TEMP2_BED_FN
   rm $TEMP_BED_FN

   set curtime=`date`
   echo "#-> Creating BigWig track of end counts scaled to 10M alignments(=reads) ($curtime)"
   ${BEDTOOLS_DIR}/genomeCoverageBed -bg -i $TEMP2_BED_FN -g ${CHROMSIZES_FN} | perl $SCALEBED_BIN $NUMALNS 5000000 > $TEMP_BED_FN
   $WIGTOBIGWIG_BIN $TEMP_BED_FN $CHROMSIZES_FN $SCALED10M_SUP150_ENDS_BIGWIG_FN

   rm $TEMP_BED_FN $TEMP2_BED_FN
endif


### Step 6. Count fragment coverage over TSS and genebodies
##if (! -e $TSS_COUNTS_SUB100NT_FN) then
##   set curtime=`date`
##   echo "${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUB100_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_tss1kbwin.bed > $TSS_COUNTS_SUB100NT_FN ($curtime)"
##   ${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUB100_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_tss1kbwin.bed > $TSS_COUNTS_SUB100NT_FN
##endif
##
##if (! -e $GENEBODIES_COUNTS_SUB100NT_FN) then
##   set curtime=`date`
##   echo "${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUB100_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_genebodies.bed > $GENEBODIES_COUNTS_SUB100NT_FN ($curtime)"
##   ${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUB100_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_genebodies.bed > $GENEBODIES_COUNTS_SUB100NT_FN
##endif
##
##if (! -e $TSS_COUNTS_SUP150NT_FN) then
##   set curtime=`date`
##   echo "${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUP150_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_tss1kbwin.bed > $TSS_COUNTS_SUP150NT_FN ($curtime)"
##   ${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUP150_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_tss1kbwin.bed > $TSS_COUNTS_SUP150NT_FN
##endif
##
##if (! -e $GENEBODIES_COUNTS_SUP150NT_FN) then
##   set curtime=`date`
##   echo "${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUP150_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_genebodies.bed > $GENEBODIES_COUNTS_SUP150NT_FN ($curtime)"
##   ${BEDTOOLS_DIR}/coverageBed -counts -a $BED_SUP150_FN -b ${FEATUREPOS_DIR}/mm10_igenomes_genebodies.bed > $GENEBODIES_COUNTS_SUP150NT_FN
##endif

if (! -e $FRAGLEN_PDF_FN) then
   set curtime=`date`
   echo "# Plot fragment length distribution  ($curtime)"
   awk '{print $3 - $2}' $BED_FN > $FRAGLEN_TMP_FN
   $RSCRIPT_BIN $RSCRIPT_LENDISTR_FN $FRAGLEN_TMP_FN $FRAGLEN_PDF_FN $FRAGLEN_LOG_PDF_FN
endif


rmdir $scratchdir

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
