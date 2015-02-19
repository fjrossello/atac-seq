#!/usr/local/bin/perl
#fpd140120_2110 
# purpose: given a BED file of fragments (each line = 1 pair of reads spanning from start of one to end of the other)
#          create a BED file wth 2 entries each representing each end


use strict;
use warnings;
main() ;

sub main {

   my $read_len = "15" ;

   my $n = 1;
   while (my $line = <STDIN>) {
#chr1    3000147 3000215 .       0       +
      chomp $line;
      my @t = split(/\t/, $line) ;
      my ($chr, $start, $end) = ($t[0], $t[1], $t[2]) ;

      my $strand = '+';
      if ($#t >= 5) {$strand = $t[5];}

      my $rev_strand = '-' ;
      if ($strand eq '-') {$rev_strand = '+';}

      print join("\t", $chr, $start, $start + $read_len, "read$n/1", 0, $strand)."\n" ;
      print join("\t", $chr, $end - $read_len, $end, "read$n/2", 0, $rev_strand)."\n" ;
      $n++ ;
   }

}
