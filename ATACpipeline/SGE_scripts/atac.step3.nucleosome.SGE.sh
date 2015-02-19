#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step3.nucleosome
#$ -e SGEOUT.atac.step3.nucleosome
#$ -r y
#$ -pe batch 16
#$ -t 1-6


set SAMPLEINFO_FN="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact/data/amo_samples_ATAC.txt"
set SAMPLE_OPTION="dataset run201409-merge"

set flowcells=( `perl txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
set rawfastq1_fns=( `perl txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r1 $SAMPLE_OPTION` )
set rawfastq2_fns=( `perl txt2tasklist.pl $SAMPLEINFO_FN raw_fastq_name_r2 $SAMPLE_OPTION` )
set names=( `perl txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )


set mouse_chr=( 'chr1' 'chr10' 'chr11' 'chr12' 'chr13' 'chr14' 'chr15' 'chr16' 'chr17' 'chr18' 'chr19' 'chr2' 'chr3' 'chr4' 'chr5' 'chr6' 'chr7' 'chr8' 'chr9' 'chrM' 'chrX' 'chrY' )

set FLOWCELL=$flowcells[$SGE_TASK_ID]
set RAWFASTQ1_FN=$rawfastq1_fns[$SGE_TASK_ID]
set RAWFASTQ2_FN=$rawfastq2_fns[$SGE_TASK_ID]
set NAME=$names[$SGE_TASK_ID]


set BASEDIR="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact"
set TRIM_FASTQF_DIR="$BASEDIR/results/ATACseq.20140915/trim_fastq/{$FLOWCELL}"
set TRIM_FASTQF1="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r1_trim.fastq.gz"
set TRIM_FASTQF2="$TRIM_FASTQF_DIR/{$FLOWCELL}_{$NAME}_r2_trim.fastq.gz"

set TRACK_DIR="$BASEDIR/results/ATACseq.20140915/tracks/${FLOWCELL}"
set VIZ_DIR="$BASEDIR/results/ATACseq.20140915/viz/${FLOWCELL}"
set SAMTOOLS_BIN="/groups/eddy/home/davisf/bin/samtools-0.1.19/samtools"
set BEDTOOLS_DIR="/groups/eddy/home/davisf/bin/bedtools-v2.15.0"
set WIGTOBIGWIG_BIN="/groups/eddy/home/davisf/software/bigWig/wigToBigWig"
set SCALEBED_BIN="/groups/eddy/home/davisf/work/genomics/fpdgenlib/fpdgenlib_work/trunk/src/scripts/scale_bed_rpm.pl"
set CHROMSIZES_FN="$BASEDIR/data/external/rsem/mm10_genome_size.tab"

set DANPOS_BIN="/groups/eddy/home/davisf/software/danpos/danpos-2.1.3/danpos.py"
set RSCRIPT_BIN="/groups/eddy/home/davisf/software/R/R-2.15.0/bin/Rscript"

set FRAGLENFIT_SCRIPT_FN="$BASEDIR/analysis/src/atacseq_fitLengthDistr.R"
set FRAGLENFIT_OUTBASE="$VIZ_DIR/${NAME}.fragLenFit"
set FRAGLENFIT_OUTFN="${FRAGLENFIT_OUTBASE}.components.txt"
set FRAGLENFIT_SUBSAMPLER_FN="$BASEDIR/analysis/src/getSubsampleBEDLengths.pl"
set FRAGCUTTER_FN="$BASEDIR/analysis/src/cutMonoNucFrags.pl"

## OUTPUT FILES
set BED_FN="$VIZ_DIR/${NAME}.bed"
set BED_SUB100_FN="$VIZ_DIR/${NAME}.sub100nt.bed"
set BED_SUP150_FN="$VIZ_DIR/${NAME}.sup150nt.bed"

set BED_MONONUC_UNSORT_FN="$VIZ_DIR/${NAME}.mononuc.unsort.bed"
set BED_MONONUC_FN="$VIZ_DIR/${NAME}.mononuc.bed"
set BAM_MONONUC_FN="$VIZ_DIR/${NAME}.mononuc.bam"
set BAI_MONONUC_FN="$VIZ_DIR/${NAME}.mononuc.bam.bai"

set GENOME_SIZE_FN="/groups/eddy/home/davisf/work/databases/iGenomes/Mus_musculus/UCSC/mm10/Sequence/Chromosomes/genome_size.tab"


set TWOBIT_FN="/groups/eddy/home/davisf/work/databases/ucsc_genomes/mm10/fasta/mm10.2bit"

set DANPOS_BASE_DIR="$BASEDIR/results/ATACseq.20140915/DANPOS_20140331/${FLOWCELL}/${NAME}"
set DANPOS_OUT_DIR="$DANPOS_BASE_DIR/result/bgsub"

set FORMAT_DANPOS_PEBED_BIN="perl /groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact/analysis/src/format_DANPOS_pairedend_BED.pl"

set BED_DANPOS_SUB100_FN="$DANPOS_BASE_DIR/${NAME}.sub100nt_danpos.bed"
set BED_DANPOS_MONONUC_FN="$DANPOS_BASE_DIR/${NAME}.mononuc_danpos.bed"
set BED_DANPOS_SUB100_BASE_FN="${NAME}.sub100nt_danpos.bed"
set BED_DANPOS_MONONUC_BASE_FN="${NAME}.mononuc_danpos.bed"
set BED_DANPOS_SUB100_BASE="${NAME}.sub100nt_danpos"
set BED_DANPOS_MONONUC_BASE="${NAME}.mononuc_danpos"

set BED_CHROMSPLIT_SCRIPT_FN="$BASEDIR/analysis/src/splitBEDbyChrom.pl"

##set DANPOS_DIFF_WIG_FN="$DANPOS_BASE_DIR/result/diff/${BED_DANPOS_MONONUC_BASE}-{$BED_DANPOS_SUB100_BASE}.pois_diff.wig"
##set DANPOS_DIFF_BIGWIG_FN="$DANPOS_BASE_DIR/result/diff/${BED_DANPOS_MONONUC_BASE}-{$BED_DANPOS_SUB100_BASE}.pois_diff.bw"

##set DANPOS_BGSUB_WIG_FN="$DANPOS_BASE_DIR/result/bgsub/${BED_DANPOS_MONONUC_BASE}.bgsub.wig"
##set DANPOS_BGSUB_BIGWIG_FN="$DANPOS_BASE_DIR/result/bgsub/${BED_DANPOS_MONONUC_BASE}.bgsub.bw"

set DANPOS_PEAKS_FN="$DANPOS_BASE_DIR/result/pooled/${BED_DANPOS_MONONUC_BASE}.bgsub.smooth.peaks.xls"
set DANPOS_PEAKS_BED_FN="$DANPOS_BASE_DIR/${BED_DANPOS_MONONUC_BASE}.bgsub.smooth.peaks.bed"
set DANPOS_POOLED_WIG_FN="$DANPOS_BASE_DIR/result/pooled/${BED_DANPOS_MONONUC_BASE}.bgsub.smooth.wig"
set DANPOS_POOLED_BIGWIG_FN="$DANPOS_BASE_DIR/result/pooled/${BED_DANPOS_MONONUC_BASE}.bgsub.smooth.bw"

###set DANPOS_DMESG_LOG_FN="$DANPOS_BASE_DIR/${BED_DANPOS_MONONUC_BASE}.DANPOS_RUN_DMESG.txt"


foreach T_OUTDIR ( $DANPOS_BASE_DIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"

if (! -e $BAM_MONONUC_FN) then #make 150-250nt mononuc BAM file
   set curtime=`date`
   echo "# Making mono-nucleosomal fragment BED file ($curtime)"

#1. Fit fragment length distribution
   set curtime=`date`
   echo "$RSCRIPT_BIN $FRAGLENFIT_SCRIPT_FN $BED_FN $FRAGLENFIT_OUTBASE $FRAGLENFIT_SUBSAMPLER_FN ($curtime)"
   $RSCRIPT_BIN $FRAGLENFIT_SCRIPT_FN $BED_FN $FRAGLENFIT_OUTBASE $FRAGLENFIT_SUBSAMPLER_FN

#2. Run perl script to make in silico cuts.
   set curtime=`date`
   echo "#perl $FRAGCUTTER_FN $BED_FN $FRAGLENFIT_OUTFN > $BED_MONONUC_UNSORT_FN ($curtime)"
   perl $FRAGCUTTER_FN $BED_FN $FRAGLENFIT_OUTFN > $BED_MONONUC_UNSORT_FN

   set curtime=`date`
   echo "#sort -k1,1 -k2,2n -k3,3n $BED_MONONUC_UNSORT_FN > $BED_MONONUC_FN ($curtime)"
   sort -k1,1 -k2,2n -k3,3n $BED_MONONUC_UNSORT_FN > $BED_MONONUC_FN

   set curtime=`date`
   echo "#rm $BED_MONONUC_UNSORT_FN ($curtime)"
   rm $BED_MONONUC_UNSORT_FN

#   awk '{if (($3 - $2) >= 150 && ($3 - $2 ) <= 250) print}' $BED_FN > $BED_MONONUC_FN
   set curtime=`date`
   echo "# ${BEDTOOLS_DIR}/bedtools bedtobam -i $BED_MONONUC_FN -g $GENOME_SIZE_FN > $BAM_MONONUC_FN ($curtime)"
   ${BEDTOOLS_DIR}/bedtools bedtobam -i $BED_MONONUC_FN -g $GENOME_SIZE_FN > $BAM_MONONUC_FN
endif

if (! -e $BAI_MONONUC_FN) then
   set curtime=`date`
   echo "# $SAMTOOLS_BIN index $BAM_MONONUC_FN ($curtime)"
   $SAMTOOLS_BIN index $BAM_MONONUC_FN
endif

if (! -e $BED_DANPOS_SUB100_FN) then
   set curtime=`date`
   echo "# $FORMAT_DANPOS_PEBED_BIN < $BED_SUB100_FN > $BED_DANPOS_SUB100_FN ($curtime)"
   $FORMAT_DANPOS_PEBED_BIN < $BED_SUB100_FN > $BED_DANPOS_SUB100_FN
endif

if (! -e $BED_DANPOS_MONONUC_FN) then
   set curtime=`date`
   echo "# $FORMAT_DANPOS_PEBED_BIN < $BED_MONONUC_FN > $BED_DANPOS_MONONUC_FN ($curtime)"
   $FORMAT_DANPOS_PEBED_BIN < $BED_MONONUC_FN > $BED_DANPOS_MONONUC_FN
endif

set curdir=`pwd`
if (! -e $DANPOS_OUT_DIR) then 
   set curtime=`date`
   cd $DANPOS_BASE_DIR

   echo "# Splitting DANPOS input BED files by chromosome ($curtime)"

   set curtime=`date`
   echo "# perl $BED_CHROMSPLIT_SCRIPT_FN ${BED_DANPOS_MONONUC_BASE} < ${BED_DANPOS_MONONUC_BASE_FN} ($curtime)"
   perl $BED_CHROMSPLIT_SCRIPT_FN ${BED_DANPOS_MONONUC_BASE} < ${BED_DANPOS_MONONUC_BASE_FN}

   set curtime=`date`
   echo "# perl $BED_CHROMSPLIT_SCRIPT_FN ${BED_DANPOS_SUB100_BASE} < ${BED_DANPOS_SUB100_BASE_FN} ($curtime)"
   perl $BED_CHROMSPLIT_SCRIPT_FN ${BED_DANPOS_SUB100_BASE} < ${BED_DANPOS_SUB100_BASE_FN}

   foreach CUR_CHR ( $mouse_chr )
      set CUR_MONONUC_BASE_FN = "${BED_DANPOS_MONONUC_BASE}.${CUR_CHR}.bed"
      set CUR_SUB100_BASE_FN = "${BED_DANPOS_SUB100_BASE}.${CUR_CHR}.bed"

      set curtime=`date`
      echo "# python ${DANPOS_BIN} ${CUR_MONONUC_BASE_FN} -b ${CUR_MONONUC_BASE_FN}:${CUR_SUB100_BASE_FN} -x 1 -k 1 -p 1 -a 1 -d 20 --clonalcut 0 -o ./ ($curtime)"
      python ${DANPOS_BIN} ${CUR_MONONUC_BASE_FN} -b ${CUR_MONONUC_BASE_FN}:${CUR_SUB100_BASE_FN} -x 1 -k 1 -p 1 -a 1 -d 20 --clonalcut 0 -o ./
   end

   set curtime=`date`
   echo "# Merging per-chromosome DANPOS outputs into whole-genome outputs ($curtime)"
   head -1 $DANPOS_BASE_DIR/result/pooled/${BED_DANPOS_MONONUC_BASE}.chr1.bgsub.smooth.peaks.xls > $DANPOS_PEAKS_FN
   foreach CUR_CHR ( $mouse_chr )
      sed 1d $DANPOS_BASE_DIR/result/pooled/${BED_DANPOS_MONONUC_BASE}.${CUR_CHR}.bgsub.smooth.peaks.xls >> $DANPOS_PEAKS_FN
      cat $DANPOS_BASE_DIR/result/pooled/${BED_DANPOS_MONONUC_BASE}.${CUR_CHR}.bgsub.smooth.wig >> $DANPOS_POOLED_WIG_FN
   end


## RUN ALL DATA TOGETHER
##   echo "# python ${DANPOS_BIN} ${BED_DANPOS_MONONUC_BASE_FN} -b ${BED_DANPOS_MONONUC_BASE_FN}:${BED_DANPOS_SUB100_BASE_FN} -x 1 -k 1 -p 1 -a 5 -d 20 --clonalcut 0 -o ./ ($curtime)"
##   python ${DANPOS_BIN} ${BED_DANPOS_MONONUC_BASE_FN} -b ${BED_DANPOS_MONONUC_BASE_FN}:${BED_DANPOS_SUB100_BASE_FN} -x 1 -k 1 -p 1 -a 5 -d 20 --clonalcut 0 -o ./

endif

if (! -e $DANPOS_POOLED_BIGWIG_FN) then
   set curtime=`date`
   cd $DANPOS_BASE_DIR
   echo "# $WIGTOBIGWIG_BIN $DANPOS_POOLED_WIG_FN $CHROMSIZES_FN $DANPOS_POOLED_BIGWIG_FN -clip ($curtime)"
   $WIGTOBIGWIG_BIN $DANPOS_POOLED_WIG_FN $CHROMSIZES_FN $DANPOS_POOLED_BIGWIG_FN -clip
endif


if (! -e $DANPOS_PEAKS_BED_FN) then
   cat $DANPOS_PEAKS_FN | sed 1d | awk -F"\t" 'OFS="\t" {print $1,$2,$3,"",$5}' > $DANPOS_PEAKS_BED_FN
endif

cd $curdir

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
