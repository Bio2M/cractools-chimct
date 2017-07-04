#! /usr/bin/perl

use Test::More tests => 11;
use CracTools::ChimCT::OverlapStructure qw(getChimeraPosMinMax);
#use CracTools::ChimCT::Chimera;

sub getChimeraPosMinMax { return CracTools::ChimCT::OverlapStructure::getChimeraPosMinMax(@_); }
sub overlapDistance { return CracTools::ChimCT::OverlapStructure::overlapDistance(@_); }
sub newChimera { return CracTools::ChimCT::Chimera->new(@_); }

my $struct = CracTools::ChimCT::OverlapStructure->new();

# Test basic methods
my ($min,$max) = getChimeraPosMinMax(newChimera('X',4,-1,'Y',2,1));
ok($min < $max);
($min,$max) = getChimeraPosMinMax(newChimera('X',2,-1,'Y',4,1));
ok($min < $max);
#is(CracTools::ChimCT::OverlapStructure::getOverlappingChimeras(

# Test class 1
my $chim1 = newChimera('X',2,-1,'Y',4,1);
$struct->addChimera($chim1);

is(overlapDistance($chim1,$chim1),0);
is(@{$struct->getOverlappingChimeras($chim1)}[0],$chim1);
is($struct->nbOverlappingChimeras($chim1),1,'nbOverlappingChimeras');

$struct->addChimera(newChimera('Y',12,1,'3',2,1));
is($struct->nbOverlappingChimeras($chim1),1,'nbOverlappingChimeras');
is($struct->nbOverlappingChimeras(newChimera('X',1,-1,'Y',1,1)),1,'nbOverlappingChimeras');

# Test class 2,3,4
my $chim2_a = newChimera('X',2,1,'X',4,1);
my $chim2_b = newChimera('X',3,1,'X',5,1);
$struct->addChimera($chim2_a);
$struct->addChimera($chim2_b);
is($struct->nbOverlappingChimeras(newChimera('X',1,1,'X',5,1)),2,'nbOverlappingChimeras');

# test overlapping distance
$struct->maxOverlappingDistance(1);
is($struct->nbOverlappingChimeras(newChimera('X',3,1,'X',6,1)),1,'nbOverlappingChimeras');

# test stranded protocol
$struct->isStranded(0);
is($struct->nbOverlappingChimeras(newChimera('X',6,-1,'X',3,-1)),1,'nbOverlappingChimeras');
is($struct->nbOverlappingChimeras(newChimera('X',6,-1,'X',3,-1)),1,'nbOverlappingChimeras');

