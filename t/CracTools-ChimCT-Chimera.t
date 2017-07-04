#! /usr/bin/perl
##
use Test::More tests => 23;
use CracTools::ChimCT::Chimera;
use CracTools::ChimCT::Read;
#use CracTools::SAMReader;

my $chimera = CracTools::ChimCT::Chimera->new('X',2,-1,'Y',4,1);

# Testing accessors
is($chimera->chr1,'X','chr1');
is($chimera->chr2,'Y','chr2');
is($chimera->pos1,2,'pos1');
is($chimera->pos2,4,'pos2');
is($chimera->strand1,-1,'strand1');
is($chimera->strand2,1,'strand2');

# testing some utility functions
ok($chimera->isSameChimera($chimera),'isSameChimera');

# Testint specific constructor that use a key
my $key = $chimera->getKey();
my $chimera_bis = CracTools::ChimCT::Chimera->newFromKey($key);
ok($chimera->isSameChimera($chimera_bis), 'newFromKey');

# lets add some read that cover the chimera
#my $sam_line1 = "read1\t0\tX\t2\t254\t10M\t*\t0\t0\tAAATTTTTTT\t*\tXE:Z:1:1:chimera:3:X|-1,2:Y|1,4";
#my $read1 = CracTools::SAMReader::SAMline->new($sam_line1);
my $read1 = CracTools::ChimCT::Read->new(id => 'read1', seq => 'AAATTTTTTT', pos_junction => 3);
$chimera->addRead($read1);

is($chimera->nbReads,1,'nbReads (1)');

#my $sam_line2 = "read2\t0\tX\t4\t254\t10M\t*\t0\t0\tAAAAATTTTT\t*\tXE:Z:1:1:chimera:5:X|-1,2:Y|1,4";
#my $read2 = CracTools::SAMReader::SAMline->new($sam_line2);
my $read2 = CracTools::ChimCT::Read->new(id => 'read2', seq => 'AAAAATTTTT', pos_junction => 5);
$chimera->addRead($read2);

is($chimera->nbReads,2,'nbReads (2)');
#my ($chr2,$pos2,$strand2) = @{$chimera->getChimeraEventFromRead($read1)->{loc2}}{'chr','pos','strand'};
#is("$chr2,$pos2,$strand2","Y,4,1",'getChimeraEventFromRead');
is($chimera->getBestRead,$read2,'getBestRead');

# Test if chimera is anchored (all the read's sequences are the same)
ok(!$chimera->isAnchored,'isAnchored (1)');
# $read2->seq($read1->seq); # Now we set the same sequences for both reads

my $chimera_2 = CracTools::ChimCT::Chimera->new('X',2,-1,'Y',4,1);
$chimera_2->addRead($read2);
$chimera_2->addRead($read2);
ok($chimera_2->isAnchored,'isAnchored (2)');

$chimera->reverseChimera();
is($chimera->chr1,'Y','chr1 reversed');
is($chimera->chr2,'X','chr2 reversed');
is($chimera->pos1,4,'pos1 reversed');
is($chimera->pos2,2,'pos2 reversed');
is($chimera->strand1,-1,'strand1 reversed');
is($chimera->strand2,1,'strand2 reversed');

# Testing chimera class
$chimera = CracTools::ChimCT::Chimera->new('X',2,-1,'Y',4,1);
is($chimera->getClass,1,'getClass (1)');
$chimera = CracTools::ChimCT::Chimera->new('X',2,1,'X',400000,1);
is($chimera->getClass,2,'getClass (2)');
$chimera = CracTools::ChimCT::Chimera->new('X',2,1,'X',1,1);
is($chimera->getClass,3,'getClass (3)');
$chimera = CracTools::ChimCT::Chimera->new('X',2,1,'X',3,-1);
is($chimera->getClass,4,'getClass (4)');
