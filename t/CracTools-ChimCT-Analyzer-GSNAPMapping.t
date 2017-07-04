#! /usr/bin/perl

use strict;
use warnings;

#use Test::More tests => 2;
#use CracTools::ChimCT::Analyzer::GSNAPMapping;
#use CracTools::ChimCT::Chimera;
#use CracTools::ChimCT::Read;
#use CracTools::ChimCT::Structure;
#
#
#my $GSNAP_EXE = '/usr/bin/gsnap';
#my $GSNAP_GENOME_DIRECTORY = '/data/indexes/gsnap';
#my $GSNAP_GENOME_NAME = 'GRCh37';
#
#my $read1 = CracTools::ChimCT::Read->new(id => 1, pos_junction => 1,seq => 'GGCAACAGGCCCAGTTCCTCTTCAGGAGAACTTCAGTGAATAAGCAAGAAATCAGGGTGAGGAAGAAAAAGGGGAGCCATGGCCTGAGCCTGGCTGGCTGG');
#
#my $chim1 = CracTools::ChimCT::Chimera->new(1,20,1,1,24,1);
#$chim1->addRead($read1);
#
#my $chim_struct = CracTools::ChimCT::Structure->new();
#$chim_struct->addChimera($chim1);
#
#my $gsnap_analyzer = CracTools::ChimCT::Analyzer::GSNAPMapping->new(chimera_struct => $chim_struct, gsnap_exe => $GSNAP_EXE, gsnap_genome_dir => $GSNAP_GENOME_DIRECTORY, gsnap_genome_name => $GSNAP_GENOME_NAME);
#
#is($gsnap_analyzer->nbReadMapped($chim1),1,'nbReadMapped');
#my ($score,$comment) = $gsnap_analyzer->getScore($chim1);
#is($score,0,'getScore');



use Test::More;
plan( skip_all => "Need GSNAP configuration" );
