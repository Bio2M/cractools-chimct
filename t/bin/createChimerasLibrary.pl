#! /usr/bin/perl
#
use strict;
use warnings;

=head1 NAME

createChimerasLibrary.pl - Create a library from multiple chimCT (csv) output

=head1 SYNOPSIS

createChimerasLibrary chimCT-file1.csv chimCT-file2.csv > library.csv

=cut

use Pod::Usage;
use Getopt::Long; # Get options
use CracTools::Output;

my @ARGV_copy = @ARGV;

if (@ARGV == 0) {
  print STDERR "Missing arguments\n";
  pod2usage(-verbose => 1);
}

my $output = CracTools::Output->new();

# Print headers in output
$output->printHeaders(version => $CracTools::ChimCT::VERSION, 
                      args => \@ARGV_copy);

# Print Header
foreach my $file (@ARGV) {
  open(my $fh,$file) or die("Cannot open $file");
  while(<$fh>) {
    next if $_ =~ /^#/;
    print $_;
  }
}
