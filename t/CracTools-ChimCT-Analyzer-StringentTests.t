#! /usr/bin/perl

use strict;
use warnings;

use Test::More tests => 2;
use CracTools::ChimCT::Analyzer::StringentTests;
use CracTools::ChimCT::Chimera;
use CracTools::ChimCT::Read;
use CracTools::ChimCT::Structure;

my $read1 = CracTools::ChimCT::Read->new(id => 1, pos_junction => 1,seq => 'ATGC');
my $chim1 = CracTools::ChimCT::Chimera->new(1,20,1,1,24,1);
$chim1->addRead($read1);
my $chim_struct = CracTools::ChimCT::Structure->new();
$chim_struct->addChimera($chim1);

my $stringent_analyzer = CracTools::ChimCT::Analyzer::StringentTests->new(chimera_struct => $chim_struct, min_support => 2, max_support => 3);

is($stringent_analyzer->minSupport,2,'minSupport');
is($stringent_analyzer->maxSupport,3,'maxSupport');
#ok(!$stringent_analyzer->_hasMinSupport($chim1),'_hasMinSupport (1)');
#$chim1->addRead($read1);
#ok($stringent_analyzer->_hasMinSupport($chim1),'_hasMinSupport (2)');
#ok(!$stringent_analyzer->_isTooSupported($chim1),'_isTooSupported');
#$chim1->addRead($read1);
#$chim1->addRead($read1);
#ok($stringent_analyzer->_isTooSupported($chim1),'_isTooSupported');
