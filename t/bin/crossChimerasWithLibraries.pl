#! /usr/bin/perl
#
use strict;
use warnings;

=head1 NAME

crossChimerasWithLibraries - Add extended fields to chimCT output

=head1 SYNOPSIS

crossChimerasWithLibraries --library lib1.csv --library lib2.csv ... chimeras.csv

=head1 DESCRIPTION



=head1 AUTHOR

Jérôme Audoux <jerome.audoux@gmail.com>

=cut

use Pod::Usage;
use Getopt::Long; # Get options
use File::Basename;
use CracTools::ChimCT::OverlapStructure;

use constant MAX_OVERLAPPING_DISTANCE => 10;

my @lib_files;
my %libs_overlap_structures;
my %libs_chimeras_samples;
my $max_overlapping_distance = MAX_OVERLAPPING_DISTANCE;

GetOptions("l|library=s" => \@lib_files,
  "m|max-overlapping-distance=i" => \$max_overlapping_distance,
);

if (@ARGV == 0) {
  print STDERR "Missing arguments\n";
  pod2usage(-verbose => 1);
}

# Loading libraries
foreach my $file (@lib_files) {
  # Get libname from file's basename
  my $lib_name = basename($file);
  # Create a new overlap structure for this library
  $libs_overlap_structures{$lib_name} = CracTools::ChimCT::OverlapStructure->new(max_overlapping_distance => $max_overlapping_distance);
  open(my $fh, $file) or die ("Cannot open $file");
  while (<$fh>) {
    # skip header lines
    next if $_ =~ /^#/;
    # Extract, sample and chimera key for the first field
    my ($sample,$chim_key) = extractSampleAndChimKey($_);
    # Add the chimera to overlap structure
    $libs_overlap_structures{$lib_name}->addChimera(CracTools::ChimCT::Chimera->newFromKey($chim_key));
    # Add the sample to the chimera list
    if(!defined $libs_chimeras_samples{$lib_name}{$chim_key}) {
      $libs_chimeras_samples{$lib_name}{$chim_key} = ();
    }
    push(@{$libs_chimeras_samples{$lib_name}{$chim_key}},$sample);
  }
}

# Crossing file with libraries
print "# crossChimerasWithLibraries.pl (libraries: ".join(",",@lib_files).")\n";
foreach my $file (@ARGV){
  open(my $fh,$file) or die("Cannot open $file");
  while(<$fh>) {
    if ($_ =~ /^#/) {
      print $_;
      next;
    }
    chomp;
    my ($sample,$chim_key) = extractSampleAndChimKey($_);
    # Look in each lib to find a hits for this chimera
    foreach my $lib (keys %libs_overlap_structures) {
      $_ .= "\tlib_$lib=";
      # Query the overlap structure of this lib to find matching chimeras
      my @overlapping_chimeras = @{$libs_overlap_structures{$lib}->getOverlappingChimeras(CracTools::ChimCT::Chimera->newFromKey($chim_key))};
      my %sample_hits;
      # if there is at least one it
      if(@overlapping_chimeras) {
        foreach my $chim_hit (@overlapping_chimeras) {
          my @chim_samples = @{$libs_chimeras_samples{$lib}{$chim_hit->getKey}};
          foreach my $sample (@chim_samples) {
            $sample_hits{$sample} = () if !defined $sample_hits{$sample};
            push(@{$sample_hits{$sample}},$chim_hit->getKey);
          }
        }
        my @sample_list = ();
        foreach my $sample (keys %sample_hits) {
          push(@sample_list,"$sample(".join(",",@{$sample_hits{$sample}}).")");
        }
        $_ .= join(",",@sample_list);
      } else {
        $_ .= "Not_found";
      }
    }
    print $_,"\n";
  }
}

sub extractSampleAndChimKey {
  my $line = shift;
  my ($id) = split("\t",$line,2);
  return split(/:/,$id);
}















