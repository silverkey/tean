#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;
use lib '.';
use VCFparser;

my $ob1sub = 'melt_bimaculoides/IIT_OCT_OB1_SUB_S7/RTE8BF6ALB220127.final_comp.vcf';
my $ob1gil = 'melt_bimaculoides/IIT_OCT_OB1_GILL_S2/RTE8BF6ALB220127.final_comp.vcf';
my $ob2sub = 'melt_bimaculoides/IIT_OCT_OB2_SUB_S5/RTE8BF6ALB220127.final_comp.vcf';
my $ob2gil = 'melt_bimaculoides/IIT_OCT_OB2_GILL_S1/RTE8BF6ALB220127.final_comp.vcf';

my $vcf1 = VCFparser->get_melt_vcf_with_meinfo($ob1sub);
#print Dumper $vcf1;

my $vcf2 = VCFparser->get_melt_vcf_with_meinfo($ob1gil);
#print Dumper $vcf2;

my $vcf3 = VCFparser->get_melt_vcf_with_meinfo($ob2sub);
#print Dumper $vcf3;

my $vcf4 = VCFparser->get_melt_vcf_with_meinfo($ob2gil);
#print Dumper $vcf4;

my $counter = 1;
my $range = 10;

foreach my $chrom1(keys %$vcf1) {
  foreach my $pos1(keys %{$vcf1->{$chrom1}}) {
    $vcf1->{$chrom1}->{$pos1}->{CLUSTER} = 'NA' unless exists $vcf1->{$chrom1}->{$pos1}->{CLUSTER};
    $vcf1->{$chrom1}->{$pos1}->{SAMPLE} = 'OB1_SUB' unless exists $vcf1->{$chrom1}->{$pos1}->{SAMPLE};
    search_in_vcf($chrom1,$pos1,$vcf1,$vcf2,'OB1_GIL');
    search_in_vcf($chrom1,$pos1,$vcf1,$vcf3,'OB2_SUB');
    search_in_vcf($chrom1,$pos1,$vcf1,$vcf4,'OB2_GIL');
  }
}

foreach my $chrom2(keys %$vcf2) {
  foreach my $pos2(keys %{$vcf2->{$chrom2}}) {
    $vcf2->{$chrom2}->{$pos2}->{CLUSTER} = 'NA' unless exists $vcf2->{$chrom2}->{$pos2}->{CLUSTER};
    $vcf2->{$chrom2}->{$pos2}->{SAMPLE} = 'OB1_GIL' unless exists $vcf2->{$chrom2}->{$pos2}->{SAMPLE};
    search_in_vcf($chrom2,$pos2,$vcf2,$vcf3,'OB2_SUB');
    search_in_vcf($chrom2,$pos2,$vcf2,$vcf4,'OB2_GIL');
  }
}

foreach my $chrom3(keys %$vcf3) {
  foreach my $pos3(keys %{$vcf3->{$chrom3}}) {
    $vcf3->{$chrom3}->{$pos3}->{CLUSTER} = 'NA' unless exists $vcf3->{$chrom3}->{$pos3}->{CLUSTER};
    $vcf3->{$chrom3}->{$pos3}->{SAMPLE} = 'OB2_SUB' unless exists $vcf3->{$chrom3}->{$pos3}->{SAMPLE};
    search_in_vcf($chrom3,$pos3,$vcf3,$vcf4,'OB2_GIL');
  }
}

foreach my $chrom4(keys %$vcf4) {
  foreach my $pos4(keys %{$vcf4->{$chrom4}}) {
    $vcf4->{$chrom4}->{$pos4}->{CLUSTER} = 'NA' unless exists $vcf4->{$chrom4}->{$pos4}->{CLUSTER};
    $vcf4->{$chrom4}->{$pos4}->{SAMPLE} = 'OB2_GIL' unless exists $vcf4->{$chrom4}->{$pos4}->{SAMPLE};
  }
}

my $header = join("\t",'sample','samples','cluster','filter','pos','chrom','gt','dp','mestart','mend','mestrand','svleng','tsd','lp','rp','sr','diff');
print $header."\n";
print_data($vcf1,'OB1_SUB');
print_data($vcf2,'OB1_GIL');
print_data($vcf3,'OB2_SUB');
print_data($vcf4,'OB2_GIL');

sub search_in_vcf {
  my $qchrom = shift;
  my $qpos = shift;
  my $qvcf = shift;
  my $tvcf = shift;
  my $tsample = shift;
  my $cluster;
  if(exists $tvcf->{$qchrom}) {
    foreach my $tpos(keys %{$tvcf->{$qchrom}}) {
      if($tpos <= ($qpos+$range) and $tpos >= ($qpos-$range)) {
        if($qvcf->{$qchrom}->{$qpos}->{CLUSTER} ne 'NA') {
          $cluster = $qvcf->{$qchrom}->{$qpos}->{CLUSTER};
        }
        else {
          $cluster = "$counter";
          $counter ++;
        }
        $qvcf->{$qchrom}->{$qpos}->{CLUSTER} = $cluster;
        $tvcf->{$qchrom}->{$tpos}->{CLUSTER} = $cluster;
        $qvcf->{$qchrom}->{$qpos}->{SAMPLE} = $qvcf->{$qchrom}->{$qpos}->{SAMPLE}.";$tsample";
        $tvcf->{$qchrom}->{$tpos}->{SAMPLE} = $qvcf->{$qchrom}->{$qpos}->{SAMPLE};
      }
    }
  }
}

sub print_data {
  my $vcf = shift;
  my $sample = shift;
  foreach my $chrom (keys %$vcf) {
    foreach my $pos (keys %{$vcf->{$chrom}}) {
      print join("\t",
        $sample,
        $vcf->{$chrom}->{$pos}->{SAMPLE},
        $vcf->{$chrom}->{$pos}->{CLUSTER},
        $vcf->{$chrom}->{$pos}->{filter},
        $vcf->{$chrom}->{$pos}->{pos},
        $vcf->{$chrom}->{$pos}->{chrom},
        $vcf->{$chrom}->{$pos}->{genotype}->{GT},
        $vcf->{$chrom}->{$pos}->{genotype}->{DP},
        $vcf->{$chrom}->{$pos}->{mestart},
        $vcf->{$chrom}->{$pos}->{meend},
        $vcf->{$chrom}->{$pos}->{mestrand},
        $vcf->{$chrom}->{$pos}->{info}->{SVLEN},
        $vcf->{$chrom}->{$pos}->{info}->{TSD},
        $vcf->{$chrom}->{$pos}->{info}->{LP},
        $vcf->{$chrom}->{$pos}->{info}->{RP},
        $vcf->{$chrom}->{$pos}->{info}->{SR},
        $vcf->{$chrom}->{$pos}->{info}->{DIFF},
      );
      print "\n";
    }
  }
}

__END__



          'Scaffold15862' => {
                               '233463' => {
                                             'ref' => 'C',
                                             'meend' => '3206',
                                             'pos' => '233463',
                                             'info' => {
                                                         'SR' => '2',
                                                         'SVTYPE' => 'RTE8BF6ALB220127',
                                                         'RP' => '38',
                                                         'LP' => '1',
                                                         'SVLEN' => '1872',
                                                         'MEINFO' => 'RTE8BF6ALB220127,1334,3206,-',
                                                         'ASSESS' => '3',
                                                         'PRIOR' => 'false',
                                                         'TSD' => 'null',
                                                         'DIFF' => '0.07:n1-1334,a1347g,a1390g,a1394g,t1431c,t1443g,g1459a,t1460g,c1471a,c1479t,a1515g,c1522t,c1532a,t1548c,t1563c,a1566t,g1569a,c1575t,t1605c,g1611t,a1612t,g1618a,c1625t,n1631-3206',
                                                         'RA' => '-5.248',
                                                         'INTERNAL' => 'Ocbimv22008358m.g,TERMINATOR'
                                                       },
                                             'mestrand' => '-',
                                             'id' => '.',
                                             'genotype' => {
                                                             'DP' => '58',
                                                             'GT' => '0/1',
                                                             'GL' => '-280,-34.92,-360',
                                                             'AD' => '28'
                                                           },
                                             'meid' => 'RTE8BF6ALB220127',
                                             'filter' => 'lc',
                                             'mestart' => '1334',
                                             'chrom' => 'Scaffold15862',
                                             'CLUSTER' => 54,
                                             'qual' => '.',
                                             'alt' => '<INS:ME:RTE8BF6ALB220127>'
                                           }
                             },


