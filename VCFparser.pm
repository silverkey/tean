#!/usr/bin/perl

package VCFparser;

use strict;
use warnings;

=head1 get_vcf
 Title   : get_vcf
 Usage   : my $vcf_aref = VCFparser->get_vcf($vcf_file);
 Function: Read and parse a VCF and returns an hashref
 Returns : arrayref
 Args    : VCF input file to read and parse
 Note    : 
=cut

sub read_vcf {
  my $class = shift;
  my $file = shift;
  my $data = [];
  open(IN,$file);
  while(my $row = <IN>) {
    next if $row =~ /^\#/;
    chomp($row);
    my $aref = {};
    my @cell = split(/\t/,$row);
    $aref->{chrom} = $cell[0];
    $aref->{pos} = $cell[1];
    $aref->{id} = $cell[2];
    $aref->{ref} = $cell[3];
    $aref->{alt} = $cell[4];
    $aref->{qual} = $cell[5];
    $aref->{filter} = $cell[6];
    $aref->{info} = _get_info($cell[7]);
    $aref->{genotype} = _get_genotype($cell[8],$cell[9]);
    push(@$data,$aref);
  }
  return $data;
}

sub _get_info {
  my $info = shift;
  my $href = {};
  my @info = split(/\;/,$info);
  foreach my $ele(@info) {
    my($key,$val) = split(/\=/,$ele);
    $href->{$key} = $val;
  }
  return $href;
}

sub _get_genotype {
  my $gtkey = shift;
  my $gtval = shift;
  my $href = {};
  my @gtkey = split(/\:/,$gtkey);
  my @gtval = split(/\:/,$gtval);

# print "\nWARNINGS UNPAIRED KEYS/VALUES FOR GENOTYPE INFO: $gtkey --- $gtval\n"
#       unless scalar(@gtkey) == scalar(@gtval);

  foreach(0..scalar(@gtkey)-1) {
    $href->{$gtkey[$_]} = $gtval[$_] || 'NA';
  }
  return $href;
}


=head1 get_melt_vcf_with_meinfo
 Title   : get_melt_vcf_with_meinfo
 Usage   : my $vcf_href = VCFparser->get_melt_vcf_with_meinfo($vcf_file);
 Function: Read and parse a VCF produced using MELT and returns an hashref
           containing also the information MEINFO in the first nesting
           of the structure (in addition to the 'info' nested hashref).
           The hasref is indexed on chromosome -> position.
 Returns : hashref
 Args    : VCF input file to read and parse
 Note    :
=cut

sub get_melt_vcf_with_meinfo {
  my $class = shift;
  my $file = shift;
  my $vcf = VCFparser->read_vcf($file);
  my $vcfchrom = {};
  foreach my $mei(@$vcf) {
    ($mei->{meid},$mei->{mestart},$mei->{meend},$mei->{mestrand}) = split(/\,/,$mei->{info}->{MEINFO});
    $vcfchrom->{$mei->{chrom}}->{$mei->{pos}} = $mei;
  }
  return $vcfchrom;
}

1;

__END__
