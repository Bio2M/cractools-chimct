#! /usr/bin/perl
#
use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;

=head1 NAME

extractChimeras.pl - Create chimeras from a chimCT with specific criterias

=head1 SYNOPSIS

  extractChimeras.pl chimeras.tsv [--type STR] [--print STR] [--min-chim-value INT] [--min-spanning-reads INT]

  Options:

    --type=STR                Parametrized himeric type to extract. (default : "all")
                              Possible values: read-through, tandem-repeat, class2, class3, class4

    --print=STR               Print mode (default : "full")
                              Possible values: full, transcript-ids

    --min-chim-value=INT      Threshold to filter chimeras based on their chim-value
    --min-spanning-reads=INT  Threshold to filter chimeras based on the number of spanning reads

  Read-through options: (Only used if type is set to "read-through")

    --max-read-through-dist=INT Maximum distance used to detect read-through. 
    --max-exon-end-dist=INT     Maximum distance from annotated exons ends.

=head1 AUTHOR

Jérôme Audoux C<jerome.audoux@inserm.fr>

=cut

use CracTools::Utils;

my $sep_field = "---";
my $min_chim_value = 75;
my $min_spanning_reads = 3;
my $max_read_through_dist = 300000;
my $max_exon_end_dist = 20;
my $mode = "all";
my $print_mode = "full";

GetOptions("type=s"                   => \$mode,
           "print=s"                  => \$print_mode,
           "min-chim-value=i"         => \$min_chim_value,
           "min-spanning-reads=i"     => \$min_spanning_reads,
           "max-read-through-dist=i"  => \$max_read_through_dist,
           "max-exon-end-dist=i"      => \$max_exon_end_dist,
) or pod2usage(0);

my $chim_file = shift @ARGV;

pod2usage(0) if !defined $chim_file;

my $chim_it = CracTools::Utils::chimCTFileIterator($chim_file);

while (my $chim = $chim_it->()) {
  # Test chimeras for read-through
  my ($trans_left,$trans_right) = split $sep_field, $chim->{extended_fields}->{mRNA_IDs};
  my ($exon_left,$exon_right) = split $sep_field, $chim->{extended_fields}->{Exon_IDs};

  # Skip chimeras with low counts
  next if $chim->{extended_fields}->{Nb_spanning_reads} < $min_spanning_reads;

  # Skip chimeras with low chim_value
  next if $chim->{chim_value} < $min_chim_value;

  my ($annot_left,$annot_right) = split $sep_field, $chim->{extended_fields}->{Annotation_type};
  my ($dist_left,$dist_right)   = split $sep_field, $chim->{extended_fields}->{Exon_end_distances};

  if($mode eq "read-through") {
    if($chim->{class} == 2 &&
      abs($chim->{pos1} - $chim->{pos2}) <= $max_read_through_dist &&
      #$annot_left eq $annot_right && $annot_left eq "EXON" &&
      $dist_left <= $max_exon_end_dist && $dist_right <= $max_exon_end_dist) {
      #print "$trans_left\n$trans_right\n";
    } else {
      next;
    }
  } elsif($mode eq "tandem-repeat") {
    if($chim->{class} == 3 &&
      $annot_left eq $annot_right && $annot_left eq "EXON" &&
      $exon_left eq $exon_right &&
      defined $chim->{comments}->{FusionDistance} && $chim->{comments}->{FusionDistance} eq "overlap") {
      #print "$trans_left\n";
      #print $chim->{_original_line}."\n";
    } else {
      next;
    }
  } elsif($mode eq "class3") {
    if($chim->{class} == 3) {
      #print "$trans_left\n";
    } else {
      next;
    }
  } elsif($mode eq "class2") {
    if($chim->{class} == 2) {
      #print "$trans_left\n$trans_right\n";
    } else {
      next;
    }
  } elsif($mode eq "class4") {
    if($chim->{class} == 4) {
      #print "$trans_left\n$trans_right\n";
    } else {
      next;
    }
  } elsif($mode eq "class1") {
    if($chim->{class} == 1) {
      #print "$trans_left\n$trans_right\n";
    } else {
      next;
    }
  } else {
    #print "$trans_left\n$trans_right\n";
  }
  if($print_mode eq "full") {
    print $chim->{_original_line},"\n";
  } elsif($print_mode eq "transcript-id") {
    if($trans_left ne $trans_right) {
      print "$trans_left\n$trans_right\n";
    } else {
      print "$trans_left\n";
    }
  }
}
