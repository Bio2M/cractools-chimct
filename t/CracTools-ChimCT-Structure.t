#! /usr/bin/perl

use Test::More tests => 5;
use CracTools::ChimCT::Structure;
use File::Temp;
use Scalar::Util;

my $struct = CracTools::ChimCT::Structure->new();

my $chim1 = CracTools::ChimCT::Chimera->new('X',2,-1,'Y',4,1);
$struct->addChimera($chim1);
is($struct->nbChimeras,1,'nbChimeras (1)');
my $chim2 = CracTools::ChimCT::Chimera->new('X',2,1,'X',4,1);
$struct->addChimera($chim2);
is($struct->nbChimeras,2,'nbChimeras (2)');
my $chim3 = CracTools::ChimCT::Chimera->new('X',2,-1,'Y',4,1);
$struct->addChimera($chim3);
is($struct->nbChimeras,2,'nbChimeras (3)');
my $iterator = $struct->chimerasIterator();
my $it_chim = $iterator->();
is(Scalar::Util::blessed($it_chim),"CracTools::ChimCT::Chimera",'chimerasIterator');

# Let's try to populate a structure from a SAM file
my $sam = "\@HD\tVN:1.5\tSO:coordinate\n".
          "\@SQ\tSN:ref\tLN:45\n".
          "r001\t163\tref\t7\t30\t8M2I4M1D3M\t=\t37\t39\tTTAGATAAAGGATACTG\t*\tXE:Z:1:1:chimera:7:X|-1,2:Y|1,4\n".
          "r002\t0\tref\t9\t30\t3S6M1P1I4M\t*\t0\t0\tAAAAGATAAGGATA\t*\tXE:Z:1:1:chimera:8:X|-1,2:Y|1,4\n".
          "r003\t0\tref\t9\t30\t5S6M\t*\t0\t0\tGCCTAAGCTAA\t*\tXE:Z:1:1:chimera:8:Y|-1,4:X|1,2\n".
          "r004\t0\tref\t16\t30\t6M14N5M\t*\t0\t0\tATAGCTTCAGC\t*\tXE:Z:1:1:chimera:8:X|-1,2:1|1,4\n".
          "r003\t2064\tref\t29\t17\t6H5M\t*\t0\t0\tTAGGC\t*\tSA:Z:ref,9,+,5S6M,30,1;\n".
          "r001\t83\tref\t37\t30\t9M\t=\t7\t-39\tCAGCGGCAT\t*\tNM:i:1\n";

my $sam_file = new File::Temp( SUFFIX => '.sam', UNLINK => 1);
print $sam_file $sam;
close $sam_file;

$struct = CracTools::ChimCT::Structure->new(SAM => $sam_file);
is($struct->nbChimeras,2,'loading SAM file');
