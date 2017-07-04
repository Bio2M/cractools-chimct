#!/usr/bin/perl
#
# ABSTRACT: extracts Tags from chimCT file / list of files 

=head1 SYNOPSIS

    chimCTtoTag.pl -f file_chimCT.tsv > output.fa

=head1 DESCRIPTION

    This tools extracts tags from chimCT output.

=cut

use strict;
use warnings;

use Getopt::Long qw(:config auto_version); # Get options
use Pod::Usage;   # Printing pod documentation in terminal
use Data::Dumper;
use CracTools::Utils;

=head1 OPTIONS

=head2 General

    --help      Print this usefull help
    --man       Open man page 

=cut

=head2 Mandatory arguments

    -f, --file           ChimCT file
    -l, --tag-length     Tag length

=cut

=head2 Syntax fasta header
    
    >id:chimKey:mean_crac_score:chimScore:chimValue:chimName

=cut

my @ARGV_copy = @ARGV;

# Variables
my ($help,$man);
my($chimct_files,$tag);

# DEFAULT
$tag = 11;

GetOptions("man"        => \$man,
    "h|help"            => \$help,
    "l|tag-length=i"    => \$tag,

# Input
    "f|file=s" => \$chimct_files,

) or pod2usage(-verbose => 1);

pod2usage(-verbose => 1)  if ($help);
pod2usage(-verbose => 2)  if ($man);

pod2usage(

    -message => "Mandatory argument '--file' is missing",
    -verbose => 1,

) unless defined $chimct_files;

open IN,$chimct_files;
while(<IN>) {

    next if $_=~/^#/;
    chomp $_ ;
    my($id,$chimname,$ch1,$p1,$s1,$ch2,$p2,$s2,$cv,$sj,$sp,$class,$comments,@others) = split(/\s+/, $_);
    my %extend_fields;
    my %com_fields;
    my ($chimScore, $cracScore);
    foreach my $e (@others) {

        my ($key, $value) = split("=",$e);
        $extend_fields{$key} = $value;
    }

    if(defined $comments) {

        my @com = split(",",$comments);
        foreach(@com) {
            my($key,$value) = split("=",$_);
            $com_fields{$key} = $value;
        } 
    }

    # Print Fasta Header
    print ">".$id.":";
    print "$com_fields{mean_crac_score}" unless !defined $com_fields{mean_crac_score};
    print ":";
    print "$com_fields{ChimScore}" unless!defined $com_fields{ChimScore};
    print ":";
    print "$cv:";
    print "$chimname\n";

    my @seqA = split('',substr($extend_fields{Read_seq},0,$extend_fields{Pos_junction}));
    my @seqB = split('',substr($extend_fields{Read_seq},$extend_fields{Pos_junction}));
    my $posA = $extend_fields{Pos_junction} - $tag;
    my $seqLeft = join('',@seqA[$posA..$extend_fields{Pos_junction}-1]);   
    my $seqRight = join('',@seqB[0..$tag-1]);

    # Print tag
    print join('',$seqLeft,$seqRight)."\n";
}



