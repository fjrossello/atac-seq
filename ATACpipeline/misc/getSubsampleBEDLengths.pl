#!/usr/local/bin/perl
#fpd140328_1610 
# goal: subsample BED file and compute length of elements

use strict;
use warnings;

while (my $line = <STDIN>) {
   if (rand() > $ARGV[0]) {next;}
   chomp $line;
   my @t = split(/\t/, $line);
   my $len = $t[2]-$t[1] ;
   print $len."\n";
}
