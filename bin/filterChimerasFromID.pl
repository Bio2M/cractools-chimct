#! /usr/bin/perl
## author: nicolas PHILIPPE
## email: nicolas.philippe@inserm.fr
## goal: filter chimeras from a list of chimeras_id

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use CracTools::Output;

=head1 NAME

    filterChimeras.pl - Filter chimeras from a list of chimeras_id

=head1 SYNOPSIS

    filterChimeras.pl [-h] [--output <output_file>] <chimCT_file> <id_file>
    The mandatory arguments are an output file generated by chimCT and a file contaning chimeras_id (one id by line).
    --output: put the filename to write the results (STDOUT by default)
 
=head1 OPTIONS
    
    -h, --help

=head1 OUTPUT FORMAT

    A chimCT output file filtered

=back

=head1 DESCRIPTION
    
    This script permets to filter chimeras from a list of chimeras_id

=head1 REQUIRES

    Perl5.
    Getopt::Long
    Pod::Usage
    
=head1 AUTHOR

    Nicolas PHILIPPE <nicolas.philippe@inserm.fr>

=cut   


my ($help, $output_file) = (0, undef);

GetOptions( "output=s" => \$output_file,
	    "help" => \$help,
        ) or pod2usage(-verbose => 0);

pod2usage(-verbose => 0)  if ($help);


if (scalar @ARGV < 2){
    pod2usage(-verbose => 0);
    exit 1;
}


open(IN,$ARGV[1]) or die("enable to open $ARGV[1]");

my %hash;
while(<IN>){
    chomp;
    if ($_ =~ /^\S+:\S+/){
	my ($id) = $_ =~ /^\S+:(\S+)/;
	$hash{$id}=1;
    }
}
close(IN);

open(IN,$ARGV[0]) or die("enable to open $ARGV[1]");

my $output = CracTools::Output->new(file => $output_file);
$output->printHeaders(args => \@ARGV);
$output->printHeaderLine("Goal: filter chimeras in from a list of chimeras_id");
$output->printHeaderLine();

my $nb_filtered = 0;
while(<IN>){
    chomp;
    if ($_ !~ /^#/){
	my ($id) = $_ =~ /^\S+:(\S+)/;
	if (defined $hash{$id}){ 
	    $nb_filtered++;
	    print STDERR "$_\n";
	}else{
	    $output->printLine($_);
	}
    }else{
	$output->printLine($_);
    }
}
close(IN);

print STDERR "$nb_filtered chimeras are filtered\n";
