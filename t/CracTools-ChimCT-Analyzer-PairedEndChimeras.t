#! /usr/bin/perl

use strict;
use warnings;

use Test::More tests => 1;
use File::Temp;
use CracTools::ChimCT::Analyzer::PairedEndChimeras;
use CracTools::ChimCT::Chimera;
use CracTools::ChimCT::Structure;

my	$line	=	'15	read1	1|1,24	pos_location=2	GATC	12,13,12,11	1,1,1,1	read2	1|1,20	pos_location=3	GGCA 11,15,12,13 1,1,1,1';

ok(1);
#my $paired_analyzer = CracTools::ChimCT::Analyzer::PairedEndChimeras->new(
#  chimera_struct => $chimera_struct,
#  paired_end_chimeras_file => $paired_end_chimeras_file,
#);

#my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = CracTools::ChimCT::Analyzer::PairedEndChimeras::_extractChimericCoordinates($line);
#
#is($chr1,1,'_extractChimericCoordinates (1)');
#is($pos1,24,'_extractChimericCoordinates (2)');
#is($strand1,1,'_extractChimericCoordinates (3)');
#is($chr2,1,'_extractChimericCoordinates (4)');
#is($pos2,20,'_extractChimericCoordinates (5)');
#is($strand1,1,'_extractChimericCoordinates (6)');
#
#my $paired_end_chimeras_file = new File::Temp( SUFFIX => '.paired-chimera', UNLINK => 1);
#print $paired_end_chimeras_file $line;
#close $paired_end_chimeras_file;

#$chimera_struct->addChimera($chim);
#my $p_chim1 = CracTools::ChimCT::Chimera->new(1,30,1,1,10,1);
#my @p_chimeras = ($p_chim1);

#my $chimera_struct = CracTools::ChimCT::Structure->new();
#my $paired_analyzer = CracTools::ChimCT::Analyzer::PairedEndChimeras->new(
#  chimera_struct => $chimera_struct,
#  paired_end_chimeras_file => $paired_end_chimeras_file,
#  paired_end_chimera_max_distance => 2,
#);
#
#my $chim = CracTools::ChimCT::Chimera->new(1,24,1,1,20,1);
#ok($paired_analyzer->_foundPairedEndChimera($chim,\@p_chimeras),'_foundPairedEndChimera');
