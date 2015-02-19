#!/usr/local/bin/perl

=head1 NAME

txt2tasklist.pl - create SGE tasklist from tab-delimited text file

=head1 VERSION

140620_0911
orig: fpd140516_1804

=head1 AUTHOR

Fred P. Davis, JFRC (davisf@janelia.hhmi.org)

=head1 SYNOPSIS

txt2tasklist.pl

=head1 DESCRIPTION

designed to be called from an SGE script itself.

=cut

use strict;
use warnings;

main() ;

sub main {

   my $usage = __FILE__." TXT_FILE FIELD_NAME|samplecount" ;
   my $txt_fn = $ARGV[0] || die $usage; #assumes header line at top
   my $field_name = $ARGV[1] || die $usage;

   my $sample_props = {}; # if specified, read in desired sample properties
   my $num_sample_props = 0 ;
   {
      my $j = 2 ;
      while ($j <= $#ARGV) {
         my $prop = $ARGV[$j] ; $prop =~ s/^-// ;
         $sample_props->{$prop} = $ARGV[($j + 1)] ;
         $j += 2 ;
         $num_sample_props++ ;
      }
   }
#   die ;


   open(INF, $txt_fn) ;
   my $f2i = {} ;
   {
      my $line = <INF> ; chomp $line;
      $line =~ s/\^#// ;
      my @t = split(/\t/, $line) ;
      map {$f2i->{$t[$_]} = $_;} (0 .. $#t) ;
   }

   $f2i->{samplecount} = 0 ;

   if (!defined $f2i->{$field_name}) {
      die "problem, can't find field $field_name in headers\n";
   }

   my @out ;
   while (my $line = <INF>) {
      chomp $line;
      if ($line =~ /^\#/) {next;}
      my @t = split(/\t/, $line) ;

      my $skipit = 0 ;
      if ($num_sample_props)  {
         foreach my $prop (keys %{$sample_props}) {
            if ($t[$f2i->{$prop}] ne $sample_props->{$prop}) {$skipit++; last;}
         }
      }

      if ($skipit) {next;}
      if ($t[$f2i->{$field_name}] eq '') {
         $t[$f2i->{$field_name}] = "NA";}
      push @out, $t[$f2i->{$field_name}] ;
   }

   if ($field_name eq 'samplecount') {
      if ($#out >= 0)   { print "".($#out + 1)."\n";}
      else              { print "0\n"; }
   } else {
      print join(" ", @out)."\n" ;
   }

}
