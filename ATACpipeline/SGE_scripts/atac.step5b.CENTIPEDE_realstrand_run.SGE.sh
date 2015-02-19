#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.step5b.CENTIPEDE_realstrand_run
#$ -e SGEOUT.atac.step5b.CENTIPEDE_realstrand_run
#$ -r y
#$ -pe batch 8
##$ -t 1-6

if ( ! ( $?SGE_TASK_ID ) ) then
   echo "ERROR: Have to specify -t on the command line when qsub'ing this script"
   exit
endif

source /sge/current/default/common/settings.csh
   

## Purpose: run PIQ footprint

set SAMPLEINFO_FN="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact/data/amo_samples_ATAC.txt"
set SAMPLE_OPTION="dataset run201409-merge"

set flowcells=( `perl txt2tasklist.pl $SAMPLEINFO_FN flowcell $SAMPLE_OPTION` )
set names=( `perl txt2tasklist.pl $SAMPLEINFO_FN sample_name $SAMPLE_OPTION` )

set MOTIFS_FN="/groups/eddy/home/davisf/work/databases/motifs/meme_20140123/mouse_pwms_withgenenames.meme"
set motifnames=( `grep ^MOTIF $MOTIFS_FN | cut -d" " -f2` )

set FLOWCELL=$flowcells[$SGE_TASK_ID]
set NAME=$names[$SGE_TASK_ID]



## 97 separate batches of motifs: 20 motifs each.
set NUM_MOTIF_BATCHES=97
set STARTMOTIFS=( '1' '21' '41' '61' '81' '101' '121' '141' '161' '181' '201' '221' '241' '261' '281' '301' '321' '341' '361' '381' '401' '421' '441' '461' '481' '501' '521' '541' '561' '581' '601' '621' '641' '661' '681' '701' '721' '741' '761' '781' '801' '821' '841' '861' '881' '901' '921' '941' '961' '981' '1001' '1021' '1041' '1061' '1081' '1101' '1121' '1141' '1161' '1181' '1201' '1221' '1241' '1261' '1281' '1301' '1321' '1341' '1361' '1381' '1401' '1421' '1441' '1461' '1481' '1501' '1521' '1541' '1561' '1581' '1601' '1621' '1641' '1661' '1681' '1701' '1721' '1741' '1761' '1781' '1801' '1821' '1841' '1861' '1881' '1901' '1921' )
set ENDMOTIFS=( '20' '40' '60' '80' '100' '120' '140' '160' '180' '200' '220' '240' '260' '280' '300' '320' '340' '360' '380' '400' '420' '440' '460' '480' '500' '520' '540' '560' '580' '600' '620' '640' '660' '680' '700' '720' '740' '760' '780' '800' '820' '840' '860' '880' '900' '920' '940' '960' '980' '1000' '1020' '1040' '1060' '1080' '1100' '1120' '1140' '1160' '1180' '1200' '1220' '1240' '1260' '1280' '1300' '1320' '1340' '1360' '1380' '1400' '1420' '1440' '1460' '1480' '1500' '1520' '1540' '1560' '1580' '1600' '1620' '1640' '1660' '1680' '1700' '1720' '1740' '1760' '1780' '1800' '1820' '1840' '1860' '1880' '1900' '1920' '1925' )


set BASEDIR="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact"
set VIZ_DIR="$BASEDIR/results/ATACseq.20140915/viz/${FLOWCELL}"

set RSCRIPT_BIN="/groups/eddy/home/davisf/software/R/R-2.15.0/bin/Rscript"
set CENTIPEDE_STRANDFIX_R_SCRIPT="$BASEDIR/analysis/run.20140915/atac.step5b.CENTIPEDE_strandfix_run.R"

set ALL_ENDS_BIGWIG_FN="$VIZ_DIR/${NAME}.rawcounts.ends.bw"
set ALL_ENDS_REALSTRAND_FW_BIGWIG_FN="$VIZ_DIR/${NAME}.rawcounts.fw_ends.bw"
set ALL_ENDS_REALSTRAND_RV_BIGWIG_FN="$VIZ_DIR/${NAME}.rawcounts.rv_ends.bw"

set CENTIPEDE_FP_OUTBASEDIR="$BASEDIR/results/ATACseq.20140915/CENTIPEDE/${FLOWCELL}"
set FIMO_SCAN_DIR="$BASEDIR/results/ATACseq.20140915/FIMO.scan"
set BWTOOL_BIN="/groups/eddy/home/davisf/software/bwtools/bwtool-1.0-gamma/bwtool/bwtool"
set PHYLOP_BIWGWIG_FN="/groups/eddy/home/davisf/work/databases/ucsc_genomes/mm10/phyloP60way/mm10.60way.phyloP60wayPlacental.bw"

foreach T_OUTDIR ( $CENTIPEDE_FP_OUTBASEDIR )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end


set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"


if ( ! ( $?start_motif )  ) then
   echo "# NOTE: controller node, job $JOB_ID task $SGE_TASK_ID"

   foreach motifbatch ( `seq 1 $NUM_MOTIF_BATCHES` )

      set start_motif=$STARTMOTIFS[${motifbatch}]
      set end_motif=$ENDMOTIFS[${motifbatch}]

      set curtime=`date`

      echo "Submitting for batch "$motifbatch": qsub -t "$SGE_TASK_ID" -v start_motif="$start_motif",end_motif="$end_motif" "$JOB_NAME" ("$curtime")"
      qsub -t $SGE_TASK_ID -v "start_motif=$start_motif,end_motif=$end_motif" $JOB_NAME
   end

else

   if (! ( $?end_motif ) ) then
      exit "ERROR: end_motif must be specified"
   endif

   echo "# NOTE: processing node, job $JOB_ID task $SGE_TASK_ID start_motif $start_motif end_motif $end_motif"

   set scratchdir="/scratch/davisf/nodejobtmp-$$-$JOB_ID-$SGE_TASK_ID"
   if (! -e $scratchdir) then
      mkdir -p $scratchdir
   endif

   foreach curmotif ( `seq $start_motif $end_motif` )

      set MOTIFNAME=$motifnames[$curmotif]
      set CUR_MOTIF_FIRST3=`echo $MOTIFNAME | awk '{print substr($0,1,3)}'`
      set FP_OUTDIR="$CENTIPEDE_FP_OUTBASEDIR/$NAME/$CUR_MOTIF_FIRST3/$MOTIFNAME"

      set FIMO_FN="$FIMO_SCAN_DIR/$CUR_MOTIF_FIRST3/$MOTIFNAME.fimo_pval1e-5.bed.gz"

      set BWTOOLOUT_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.bwtool.txt"
      set CUTMATRIX_REALSTRANDED_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.cutmatrix_stranded_real.txt"
      set PHYLOPMATRIX_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.phylop_score.txt"

      set FPOUTBED_WITHCONS_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.withcons_stranded_real_fix.bed"
      set FPOUTPDF_WITHCONS_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.withcons_stranded_real_fix.pdf"
      set FPOUTBED_NOCONS_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.nocons_stranded_real_fix.bed"
      set FPOUTPDF_NOCONS_FN="$FP_OUTDIR/$NAME.$MOTIFNAME.nocons_stranded_real_fix.pdf"

      set TMPBED_FN="$scratchdir/tmp.$MOTIFNAME.bed"
      set TMPBWTOOLOUT_FW_FN="$scratchdir/tmp.$MOTIFNAME.$NAME.fw.bwtools.out"
      set TMPBWTOOLOUT_RV_FN="$scratchdir/tmp.$MOTIFNAME.$NAME.rv.bwtools.out"

      set TMP_PHYLOP_BWTOOLOUT_FN="$scratchdir/tmp.$MOTIFNAME.$NAME.phylop.bwtools.out"
   

      if (! -s $FPOUTPDF_WITHCONS_FN) then
         if (! -s $FIMO_FN) then
            echo "WARNING: skipping $MOTIFNAME for $NAME because $FIMO_FN not found or empty"
         else
   
      
            if (! -e $FP_OUTDIR) then
               mkdir -p $FP_OUTDIR
            endif
      
      
            zcat $FIMO_FN | awk -F"\t" 'BEGIN{OFS="\t"} {print $1,$2,$3,$4,0,$6}' > $TMPBED_FN
      
      
            set curtime=`date`
            echo "$BWTOOL_BIN summary $TMPBED_FN $PHYLOP_BIWGWIG_FN $TMP_PHYLOP_BWTOOLOUT_FN ($curtime)"
            $BWTOOL_BIN summary $TMPBED_FN $PHYLOP_BIWGWIG_FN $TMP_PHYLOP_BWTOOLOUT_FN
            cut -f8 $TMP_PHYLOP_BWTOOLOUT_FN | sed 's/na/0/' > $PHYLOPMATRIX_FN
            if (-e $PHYLOPMATRIX_FN) then
               rm $PHYLOPMATRIX_FN.gz
            endif
            gzip $PHYLOPMATRIX_FN

            set curtime=`date`
            echo "# $BWTOOL_BIN matrix 100:100 $TMPBED_FN $ALL_ENDS_REALSTRAND_FW_BIGWIG_FN $TMPBWTOOLOUT_FW_FN ($curtime)"
            $BWTOOL_BIN matrix 100:100 $TMPBED_FN $ALL_ENDS_REALSTRAND_FW_BIGWIG_FN $TMPBWTOOLOUT_FW_FN

            set curtime=`date`
            echo "# $BWTOOL_BIN matrix 100:100 $TMPBED_FN $ALL_ENDS_REALSTRAND_RV_BIGWIG_FN $TMPBWTOOLOUT_RV_FN ($curtime)"
            $BWTOOL_BIN matrix 100:100 $TMPBED_FN $ALL_ENDS_REALSTRAND_RV_BIGWIG_FN $TMPBWTOOLOUT_RV_FN

            paste $TMPBWTOOLOUT_FW_FN $TMPBWTOOLOUT_RV_FN | sed 's/-nan/0/g' > $CUTMATRIX_REALSTRANDED_FN
            gzip $CUTMATRIX_REALSTRANDED_FN
      
            rm $TMP_PHYLOP_BWTOOLOUT_FN; rm $TMPBWTOOLOUT_FW_FN; rm $TMPBWTOOLOUT_RV_FN; rm $TMPBED_FN

      
            set curtime=`date`
            echo "# $RSCRIPT_BIN $CENTIPEDE_STRANDFIX_R_SCRIPT withcons $CUTMATRIX_REALSTRANDED_FN.gz $PHYLOPMATRIX_FN.gz $FIMO_FN $FPOUTBED_WITHCONS_FN $FPOUTPDF_WITHCONS_FN ($curtime)"
            $RSCRIPT_BIN $CENTIPEDE_STRANDFIX_R_SCRIPT withcons $CUTMATRIX_REALSTRANDED_FN.gz $PHYLOPMATRIX_FN.gz $FIMO_FN $FPOUTBED_WITHCONS_FN $FPOUTPDF_WITHCONS_FN

            set curtime=`date`
            echo "# $RSCRIPT_BIN $CENTIPEDE_STRANDFIX_R_SCRIPT nocons $CUTMATRIX_REALSTRANDED_FN.gz $FIMO_FN $FPOUTBED_NOCONS_FN $FPOUTPDF_NOCONS_FN ($curtime)"
            $RSCRIPT_BIN $CENTIPEDE_STRANDFIX_R_SCRIPT nocons $CUTMATRIX_REALSTRANDED_FN.gz $FIMO_FN $FPOUTBED_NOCONS_FN $FPOUTPDF_NOCONS_FN

         endif
      else
         echo "NOTE: skipping $MOTIFNAME for $NAME because already done: non-zero $FPOUTPDF_WITHCONS_FN"
      endif

   end

endif

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
