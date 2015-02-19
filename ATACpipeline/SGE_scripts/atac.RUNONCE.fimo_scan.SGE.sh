#!/bin/csh
#$ -S /bin/csh
#$ -cwd
#$ -o SGEOUT.atac.fimo_scan
#$ -e SGEOUT.atac.fimo_can
#$ -r y
#$ -pe batch 4
#$ -t 1-193

set NUMCPU=4

## Purpose: 1. scan for known motifs across whole genome
##          2. run CENTIPEDE to footprint motifs

set STARTTASKS=( '1' '11' '21' '31' '41' '51' '61' '71' '81' '91' '101' '111' '121' '131' '141' '151' '161' '171' '181' '191' '201' '211' '221' '231' '241' '251' '261' '271' '281' '291' '301' '311' '321' '331' '341' '351' '361' '371' '381' '391' '401' '411' '421' '431' '441' '451' '461' '471' '481' '491' '501' '511' '521' '531' '541' '551' '561' '571' '581' '591' '601' '611' '621' '631' '641' '651' '661' '671' '681' '691' '701' '711' '721' '731' '741' '751' '761' '771' '781' '791' '801' '811' '821' '831' '841' '851' '861' '871' '881' '891' '901' '911' '921' '931' '941' '951' '961' '971' '981' '991' '1001' '1011' '1021' '1031' '1041' '1051' '1061' '1071' '1081' '1091' '1101' '1111' '1121' '1131' '1141' '1151' '1161' '1171' '1181' '1191' '1201' '1211' '1221' '1231' '1241' '1251' '1261' '1271' '1281' '1291' '1301' '1311' '1321' '1331' '1341' '1351' '1361' '1371' '1381' '1391' '1401' '1411' '1421' '1431' '1441' '1451' '1461' '1471' '1481' '1491' '1501' '1511' '1521' '1531' '1541' '1551' '1561' '1571' '1581' '1591' '1601' '1611' '1621' '1631' '1641' '1651' '1661' '1671' '1681' '1691' '1701' '1711' '1721' '1731' '1741' '1751' '1761' '1771' '1781' '1791' '1801' '1811' '1821' '1831' '1841' '1851' '1861' '1871' '1881' '1891' '1901' '1911' '1921' )
set ENDTASKS=( '10' '20' '30' '40' '50' '60' '70' '80' '90' '100' '110' '120' '130' '140' '150' '160' '170' '180' '190' '200' '210' '220' '230' '240' '250' '260' '270' '280' '290' '300' '310' '320' '330' '340' '350' '360' '370' '380' '390' '400' '410' '420' '430' '440' '450' '460' '470' '480' '490' '500' '510' '520' '530' '540' '550' '560' '570' '580' '590' '600' '610' '620' '630' '640' '650' '660' '670' '680' '690' '700' '710' '720' '730' '740' '750' '760' '770' '780' '790' '800' '810' '820' '830' '840' '850' '860' '870' '880' '890' '900' '910' '920' '930' '940' '950' '960' '970' '980' '990' '1000' '1010' '1020' '1030' '1040' '1050' '1060' '1070' '1080' '1090' '1100' '1110' '1120' '1130' '1140' '1150' '1160' '1170' '1180' '1190' '1200' '1210' '1220' '1230' '1240' '1250' '1260' '1270' '1280' '1290' '1300' '1310' '1320' '1330' '1340' '1350' '1360' '1370' '1380' '1390' '1400' '1410' '1420' '1430' '1440' '1450' '1460' '1470' '1480' '1490' '1500' '1510' '1520' '1530' '1540' '1550' '1560' '1570' '1580' '1590' '1600' '1610' '1620' '1630' '1640' '1650' '1660' '1670' '1680' '1690' '1700' '1710' '1720' '1730' '1740' '1750' '1760' '1770' '1780' '1790' '1800' '1810' '1820' '1830' '1840' '1850' '1860' '1870' '1880' '1890' '1900' '1910' '1920' '1925' )

set starttask=$STARTTASKS[$SGE_TASK_ID]
set endtask=$ENDTASKS[$SGE_TASK_ID]


## OUTPUT FILES
set BASEDIR="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact"
set FIMO_BIN="/groups/eddy/home/davisf/bin/meme/bin/fimo"
set GENOMEFA="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact/data/external/rsem/mm10.fa"
set FIMO_SCAN_DIR="$BASEDIR/results/ATACseq.20140718/FIMO.scan"
set FIMO_PVAL_CUTOFF="1e-5"

set MOTIFS_FN="/groups/eddy/home/davisf/work/databases/motifs/meme/mouse_pwms_withgenenames.meme"

set scratchdir="/scratch/davisf/nodejobtmp-$$-$JOB_ID-$SGE_TASK_ID"
foreach T_OUTDIR ( $FIMO_SCAN_DIR $scratchdir )
   if (! -e $T_OUTDIR) then
      mkdir -p $T_OUTDIR
   endif
end


set MOTIF_FN="/groups/eddy/home/davisf/work/databases/motifs/meme_20140123/mouse_pwms_withgenenames.meme"
set motifnames=( `grep ^MOTIF ~/work/databases/motifs/meme/*meme | cut -d" " -f2` )

set EXTRACT_MEME_MOTIF_SCRIPT="/groups/eddy/home/davisf/work/genomics/seq_data_analysis/mouse_intact/analysis/src/extract_meme_motif.pl"

set curdir=`pwd`
set curhost=`hostname`
set curtime=`date`
echo "#sgejob run started on $curhost at $curtime"


foreach curmotif ( `seq $starttask $endtask` )
   set curtime=`date`

   set MOTIFNAME=$motifnames[$curmotif]

   set CUR_MOTIF_FN="$scratchdir/motif$curmotif.meme"
   echo "perl $EXTRACT_MEME_MOTIF_SCRIPT $MOTIFS_FN $curmotif > $CUR_MOTIF_FN"
   perl $EXTRACT_MEME_MOTIF_SCRIPT $MOTIFS_FN $curmotif > $CUR_MOTIF_FN


   set CUR_MOTIF_FIRST3=`echo $MOTIFNAME | awk '{print substr($0,1,3)}'`
   set CUR_FIMO_OUTDIR="$FIMO_SCAN_DIR/$CUR_MOTIF_FIRST3"
   if (! -e $CUR_FIMO_OUTDIR) then
      mkdir -p $CUR_FIMO_OUTDIR
   endif


   set CUR_FIMO_RAW_OUTFN="$CUR_FIMO_OUTDIR/$MOTIFNAME.fimo_pval$FIMO_PVAL_CUTOFF.txt"
   set CUR_FIMO_BED_OUTFN="$CUR_FIMO_OUTDIR/$MOTIFNAME.fimo_pval$FIMO_PVAL_CUTOFF.bed"

   if (! -s $CUR_FIMO_RAW_OUTFN) then
      echo "# $FIMO_BIN --text --output-pthresh $FIMO_PVAL_CUTOFF --max-stored-scores 500000 $CUR_MOTIF_FN $GENOMEFA > $CUR_FIMO_RAW_OUTFN ($curtime)"
      $FIMO_BIN --text --output-pthresh $FIMO_PVAL_CUTOFF --max-stored-scores 500000 $CUR_MOTIF_FN $GENOMEFA > $CUR_FIMO_RAW_OUTFN

      rm $CUR_MOTIF_FN
   endif

   if (! -s $CUR_FIMO_BED_OUTFN) then
      sed 1d $CUR_FIMO_RAW_OUTFN | awk -F"\t" 'BEGIN{OFS="\t"} {print $2,$3 - 1,$4,$1,$6,$5}' > $CUR_FIMO_BED_OUTFN
   endif

   if (-s $CUR_FIMO_BED_OUTFN) then
      gzip $CUR_FIMO_BED_OUTFN
   endif

   if (-s $CUR_FIMO_RAW_OUTFN) then
      gzip $CUR_FIMO_RAW_OUTFN
   endif


#Motif   Seq     Start   Stop   Strand  Log-odds        p-value Site
#Barx1__UP00181_1_Barx1_2877.1__homeodomain      chr10   3446066 3446081 +       16.295  1.21677e-07     AAAGTAATTGGTACAT

end
rmdir $scratchdir

set curtime=`date`
echo "#sgejob run finished on $curhost at $curtime"
