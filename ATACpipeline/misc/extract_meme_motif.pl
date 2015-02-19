#!/usr/local/bin/perl
#fpd140718_0936 
# Purpose: given MEME motif databsae, extract motif # specified as ARGV

use strict;
use warnings;

main() ;

sub main {

   my $usage = __FILE__." meme_motif_file motif_number" ;
   my $motif_fn = $ARGV[0] || die $usage;
   my $motif_num = $ARGV[1] || die $usage;

   my $motif_fh ;
   open($motif_fh, $motif_fn) ;

   my $headers = 1 ;
   my $cur_motif_num = 0 ; my $motif_found = 0 ;

   my $outlines = '';
   while (my $line = readline($motif_fh)) {
      chomp $line;
      if ($line =~ /^MOTIF/) {
         $cur_motif_num++ ;
         $headers = 0;

         if ($motif_num == $cur_motif_num) {
            $motif_found++ ;
            $outlines .= $line."\n";
            while (my $motifline = readline($motif_fh)) {
               chomp $motifline;
               if ($motifline =~ /^MOTIF/) {last;}
               $outlines .= $motifline."\n";
            }
            last;
         }
      } elsif ($headers) {
         $outlines .= $line."\n";
      }
   }
   close($motif_fh) ;

   if ($motif_found) {
      print $outlines ;
   } else {
      print STDERR "ERROR: didn't find motif $motif_num\n";
   }

}
