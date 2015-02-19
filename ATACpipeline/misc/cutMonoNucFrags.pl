=head1 NAME

cutMonoNucFrags.pl

=head1 VERSION

fpd140330_0050

=head1 AUTHOR

Fred P. Davis, JFRC (davisf@janelia.hhmi.org)

=head1 SYNOPSIS

../src/cutMonoNucFrags.pl

=head1 DESCRIPTION

Given fragment length boundaries inferred by gaussian fits to ATAC-seq fragment
length distribution, read original BED file and split to simulate 
mono-nucleosomal cuts.

=cut

use strict;
use warnings;
use POSIX qw/floor/ ;

main() ;


sub main {

   my $usage = __FILE__." bed_fn fragment_boundaries" ;

   my $specs = {} ;
   $specs->{bed_fn} = $ARGV[0] || die $usage;
   $specs->{bound_fn} = $ARGV[1] || die $usage;

   my $len2comp = read_bound_fn({fn => $specs->{bound_fn}}) ;

   open(BEDF, $specs->{bed_fn});
   while (my $line = <BEDF>) {
      chomp $line;
      my @t = split(/\t/, $line) ;

      my $start = $t[1] ;
      my $end = $t[2] ;
      my $length = $end - $start ;

      if (exists $len2comp->{$length}) {
         my $comp = $len2comp->{$length} ;
         if ($comp == 1) {
            print $line."\n" ;
         } elsif ($comp > 1) {
            my $split_length = POSIX::floor($length / $comp) ;
            my (@split_starts, @split_ends) ;
            foreach my $split (1 .. $comp) {
               my $split_start = $start + $split_length * ($split - 1) ;
               my $split_end   = $split_start + $split_length ;

               if ($split == $comp) {$split_end = $end; }

               $t[1] = $split_start ;
               $t[2] = $split_end ;
#               print STDERR "split #$split of $start-$end is $split_start - $split_end\n";
               print join("\t", @t)."\n";
            }
         }
      }

   }

}

sub read_bound_fn {

   my $in = shift ;
   my $comp2cuts = [] ;

   open(BOUNDF, $in->{fn}) ;
   my $header = <BOUNDF> ; chomp $header;
   my $f2i = {} ; my @fields ;
   {
      my @t = split(/\t/, $header) ;
      @fields = @t ;
      map {$f2i->{$t[$_]} = $_} (0 .. $#t) ;
   }
   while (my $line = <BOUNDF>) {
      chomp $line;
      my @t = split(/\t/, $line) ;
      $t[$f2i->{'n'}]-- ; #switch R numbering 1.. to perl 0..
      $comp2cuts->[$t[$f2i->{'n'}]] = {
         start => $t[$f2i->{'start_pr90'}],
         end   => $t[$f2i->{'end_pr90'}]
      } ;
   }
   close(BOUNDF) ;

   my $len2comp = {};
   foreach my $i ( 0 .. $#{$comp2cuts}) {
      if ($comp2cuts->[$i]->{start} =~ /Inf/) {next;}
      foreach my $j ( $comp2cuts->[$i]->{start} .. $comp2cuts->[$i]->{end}) {
         $len2comp->{$j} = $i ; }
   }

   return($len2comp) ;
}
