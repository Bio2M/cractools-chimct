#! /usr/bin/perl
#
use strict;
use warnings;

use Statistics::Descriptive;

=head1 NAME

chimeraStats.pl - Generate statistics from ChimCT CSV output.

=head1 SYNOPSIS

chimeraStats.pl chimeras.csv

=head1 DESCRIPTION

The chimeraStats script generate a complete report from a chimeraPipeline CSV output.

Number of chimera

Répartition par classes

Répartition par score

Croiser les classes et les score

Moyenne de span junc par class

Moyenne de span junc par score

=cut

if(@ARGV < 1) {
  die "Missing CSV file in argument";
}

my $chimera_file = shift;

open(IN,$chimera_file) or die("Cannot open $chimera_file");

my %causes = (ANCHORED        => 'chimera_anchored',
              PSEUDOGENE      => 'pseudogene',
              LOW_SUPPORT     => 'low_support',
              PE_SUPPORT      => 'strange_paired_end_support',
              HIGH_SUPPORT    => 'high_support',
            );

my %chimeras = ();

my %stats = (NB_CHIMERAS      => 0,
             RANKS            => [],
             SPAN_JUNC        => [],
             SPAN_PE          => [],
             SPAN_RATIO       => [],
            );

# INIT Stats
foreach my $class (1,2,3,4) {
  foreach my $cause (keys %causes) {
    $stats{CLASS}{$class}{$cause} = 0;
  }
  $stats{CLASS}{$class}{NB_CHIMERAS} = 0;
  $stats{CLASS}{$class}{SAME_GENE} = 0;
  $stats{CLASS}{$class}{SPAN_JUNC} = [];
  $stats{CLASS}{$class}{SPAN_PE} = [];
  $stats{CLASS}{$class}{RANKS} = [];
}

while(<IN>) {
  # skip headers
  next if $_ =~ /^#/;
  my ($id,$name,$chr1,$pos1,$strand1,$chr2,$pos2,$strand2,$rank,$span_junc,$span_PE,$class,$comments,@others) = split('\t',$_);
  $chimeras{$id} = {name => $name, rank => $rank, span_junc => $span_junc, span_PE => $span_PE, class => $class};
  $stats{NB_CHIMERAS}++; 
  push($stats{RANKS},$rank);
  push($stats{SPAN_JUNC},$span_junc);
  push($stats{SPAN_PE},$span_PE);
  if($span_PE > 0) {
    push($stats{SPAN_RATIO},$span_junc/$span_PE);
  } else {
    push($stats{SPAN_RATIO},0);
  }

  # By class
  $stats{CLASS}{$class}{NB_CHIMERAS}++;
  foreach my $cause (keys %causes) {
    $stats{CLASS}{$class}{$cause}++ if $comments =~ /$causes{$cause}/;
  }
  push($stats{CLASS}{$class}{RANKS},$rank);
  push($stats{CLASS}{$class}{SPAN_JUNC},$span_junc);
  push($stats{CLASS}{$class}{SPAN_PE},$span_PE);
  
  my ($geneA,$geneB) = split '---', $name;
  $stats{CLASS}{$class}{SAME_GENE}++ if $geneA eq $geneB && $geneA ne 'N/A';
}

my $rank_stats = Statistics::Descriptive::Full->new();
$rank_stats->add_data($stats{RANKS});

my $span_junc_stats = Statistics::Descriptive::Full->new();
$span_junc_stats->add_data($stats{SPAN_JUNC});

my $span_PE_stats = Statistics::Descriptive::Full->new();
$span_PE_stats->add_data($stats{SPAN_PE});

my $span_ratio_stats = Statistics::Descriptive::Full->new();
$span_ratio_stats->add_data($stats{SPAN_RATIO});

print "---------> GENRAL STATISTICS <-----------
Found $stats{NB_CHIMERAS} chimeras :
> Class 1: $stats{CLASS}{1}{NB_CHIMERAS}
> Class 2: $stats{CLASS}{2}{NB_CHIMERAS}
> Class 3: $stats{CLASS}{3}{NB_CHIMERAS}
> Class 4: $stats{CLASS}{4}{NB_CHIMERAS}

Ranks stats:
> mean  : ".int($rank_stats->mean+0.5)."

Ranks distribution: 
> 100   : ".scalar(grep {$_ == 100} @{$stats{RANKS}})."
> 99-80 : ".scalar(grep {$_ < 100 && $_ >= 80} @{$stats{RANKS}})."
> 89-60 : ".scalar(grep {$_ < 80 && $_ >= 60} @{$stats{RANKS}})."
> 59-40 : ".scalar(grep {$_ < 60 && $_ >= 40} @{$stats{RANKS}})."
> 39-0  : ".scalar(grep {$_ < 40} @{$stats{RANKS}})."

Spanning Reads :
> Spanning junction : 
  > Sum ".$span_junc_stats->sum."
  > Mean ".$span_junc_stats->mean."
> Spanning PE
  > Sum ".$span_PE_stats->sum."
  > Mean ".$span_PE_stats->mean."
> Mean(sp_jun/sp_PE); ".$span_ratio_stats->mean."

---------> CLASS STATISTICS <-----------
";
foreach my $class (1,2,3,4) {
  my $class_rank_stats = Statistics::Descriptive::Full->new();
  $class_rank_stats->add_data($stats{CLASS}{$class}{RANKS});

  my $class_span_junc_stats = Statistics::Descriptive::Full->new();
  $class_span_junc_stats->add_data($stats{CLASS}{$class}{SPAN_JUNC});

  my $class_span_PE_stats = Statistics::Descriptive::Full->new();
  $class_span_PE_stats->add_data($stats{CLASS}{$class}{SPAN_PE});

  my ($anchored_percent) = split '\.', $stats{CLASS}{$class}{ANCHORED} / $stats{CLASS}{$class}{NB_CHIMERAS} * 100;
  my ($pseudo_percent) = split '\.', $stats{CLASS}{$class}{PSEUDOGENE} / $stats{CLASS}{$class}{NB_CHIMERAS} * 100;

  print "Class $class : $stats{CLASS}{$class}{NB_CHIMERAS} chimeras
- Ranks
  > mean : ".$class_rank_stats->mean."
- Ranks distribution: 
  > 100   : ".scalar(grep {$_ == 100} @{$stats{CLASS}{$class}{RANKS}})."
  > 99-80 : ".scalar(grep {$_ < 100 && $_ >= 80} @{$stats{CLASS}{$class}{RANKS}})."
  > 89-60 : ".scalar(grep {$_ < 80 && $_ >= 60} @{$stats{CLASS}{$class}{RANKS}})."
  > 59-40 : ".scalar(grep {$_ < 60 && $_ >= 40} @{$stats{CLASS}{$class}{RANKS}})."
  > 39-0  : ".scalar(grep {$_ < 40} @{$stats{CLASS}{$class}{RANKS}})."
- Causes
";

foreach my $cause (keys %causes) {
  my ($cause_percent) = split '\.', $stats{CLASS}{$class}{$cause} / $stats{CLASS}{$class}{NB_CHIMERAS} * 100;
  print "  > $cause   : ".$stats{CLASS}{$class}{$cause}." ($cause_percent%)\n";
}

my ($same_gene_percent) = split '\.', $stats{CLASS}{$class}{SAME_GENE} / $stats{CLASS}{$class}{NB_CHIMERAS} * 100;
print "- Other
  > Same gene : $stats{CLASS}{$class}{SAME_GENE} ($same_gene_percent%)
";

print "- Spanning reads
  > Spanning junction : ".$class_span_junc_stats->sum." (mean : ".$class_span_junc_stats->mean.")
  > Spanning PE       : ".$class_span_PE_stats->sum." (mean : ".$class_span_PE_stats->mean.")

";
}


