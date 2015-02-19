#!/usr/local/bin/perl

=head1 NAME

splitBEDbyChrom.pl

=head1 VERSION

fpd140331_2103

=head1 AUTHOR

Fred P. Davis, JFRC (davisf@janelia.hhmi.org)

=head1 SYNOPSIS

splitBEDbyChrom.pl splits BED entries by chromosome

=head1 DESCRIPTION

=cut


use strict;
use warnings;

main() ;

sub main {

   my $usage = __FILE__." OUT_BASE < INPUT_BED file
will split the <INPUT_BED> file into <OUT_BASE>.chrXX.bed" ;

   my $out_base = $ARGV[0] || die $usage ;
   my $chr2fh = {};
   my $chr2fn = {};
   while (my $line = <STDIN>) {
      chomp $line;
      my @t = split(/\t/, $line) ;
      if (!exists $chr2fh->{$t[0]}) {
         $chr2fn->{$t[0]} = $out_base.".".$t[0].".bed" ;
         open($chr2fh->{$t[0]}, ">".$chr2fn->{$t[0]}) ;
      }
      print {$chr2fh->{$t[0]}} $line."\n";
   }

   foreach my $chr (sort keys %{$chr2fn}) {
      close($chr2fh->{$chr}) ;
      print $chr2fn->{$chr}."\n";
   }

}
