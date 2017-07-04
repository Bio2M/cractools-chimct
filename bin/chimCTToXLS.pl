#! /usr/bin/perl
#

=pod NAME

chimCTToXLS - Convert chimCT output file to XLS Workbook

=pod SYNOPSIS

chimCTToXLS file-chimeras.csv file-chimeras.xls

=head1 OPTIONS

  --field=<INT>   Choose a field to split output XLS file in multiple sheets. Default: field 12 (Class).

=pod AUTHOR

Jérôme Audoux <jerome.audoux@gmail.com>

=cut

use strict;
use warnings;

use Pod::Usage;
use Getopt::Long;
use Excel::Writer::XLSX;

my $split_field = 12;

GetOptions( "f|field=i" => \$split_field);

my $chimCT_file = shift;
my $xls_file = shift;

pod2usage(-verbose => 1) unless defined $chimCT_file and defined $xls_file;


open(my $fh, $chimCT_file) or die("Cannot open $chimCT_file");
#my @header_raw = qw(Id Name Chr1 Pos1 Strand1 Chr2 Pos2 Strand2 Rank Spanning_junction_normalized Spanning_PE_normalized Class Comments);

my $workbook  = Excel::Writer::XLSX->new($xls_file);
my %class_worksheets;
my %nb_chimeras;
my $last_header_line;
my @header_raw;

while(<$fh>) {
  chomp;
  # Skip header lines
  if ($_ =~ /^#/) {
    ($last_header_line) = $_ =~ /^#\s+(.*)/;
    next;
  }

  @header_raw = split(/\s+/,$last_header_line) if @header_raw == 0;
   
  my @fields = split("\t",$_);
  my $split_field_value = $fields[$split_field-1];

  if(!defined $class_worksheets{$split_field_value}) {
    $class_worksheets{$split_field_value} = $workbook->add_worksheet("$header_raw[$split_field-1] $split_field_value");
    $class_worksheets{$split_field_value}->write_row(0,0,\@header_raw);
    $nb_chimeras{$split_field_value} = 1;
  }

  my $nb_fields = 0;
  foreach my $field (@fields) {
    $class_worksheets{$split_field_value}->write($nb_chimeras{$split_field_value}, $nb_fields, $field);
    $nb_fields++;
  }
  $nb_chimeras{$split_field_value}++;
}


