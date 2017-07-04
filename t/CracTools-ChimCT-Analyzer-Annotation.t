#! /usr/bin/perl

use strict;
use warnings;

use Test::More tests => 28;
use CracTools::ChimCT::Analyzer::Annotation;
use CracTools::ChimCT::Structure;
use File::Temp 0.23;
use Inline::Files 0.68;
#use Data::Dumper;

# Create a temp file with the GFF lines described below
my $gff_file = new File::Temp( SUFFIX => '.gff', UNLINK => 1);
while(<GFF>) {print $gff_file $_;}
close $gff_file;

# Create a temp file with the GFF lines described below
my $est_file = new File::Temp( SUFFIX => '.gff', UNLINK => 1);
while(<EST>) {print $est_file $_;}
close $est_file;

#chim(chr1,pos1,strand1,chr2,pos2,strand2)
## basic tests
my $chim = CracTools::ChimCT::Chimera->new(1,64,1,1,60,1);
my $chim2 = CracTools::ChimCT::Chimera->new(1,58,-1,1,3,-1);
my $chim3 = CracTools::ChimCT::Chimera->new(2,46,1,2,21,1);
my $chim4 = CracTools::ChimCT::Chimera->new(2,33,1,2,69,1);
my $chim5 = CracTools::ChimCT::Chimera->new(2,110,1,2,140,1);
my $chim6 = CracTools::ChimCT::Chimera->new(2,1,1,2,41000000,1);
## bugs
my $chim_17618 = CracTools::ChimCT::Chimera->new(7,98984412,1,7,98985657,1);
my $chim_17623 = CracTools::ChimCT::Chimera->new(9,117160718,-1,9,117156692,-1);
my $chim_17647 = CracTools::ChimCT::Chimera->new(10,67818421,1,10,77818421,1);  
my $chim_17676 = CracTools::ChimCT::Chimera->new(17,62781339,-1,17,62760666,-1);  
my $chim_17994 = CracTools::ChimCT::Chimera->new(6,32551357,-1,6,32525206,-1);

#global test for classification
is($chim->getClass(),3,'getClass (1)');
is($chim2->getClass(),2,'getClass (2)');
is($chim3->getClass(),3,'getClass (3)');
is($chim4->getClass(),2,'getClass (4)');
is($chim5->getClass(),2,'getClass (5)');
is($chim6->getClass(),2,'getClass (6)');

#specific test for non stranded protocol
my $chim_struct_nonstranded = CracTools::ChimCT::Structure->new();
$chim_struct_nonstranded->addChimera($chim);
$chim_struct_nonstranded->addChimera($chim2);
my $annot_analyzer_nonstranded = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_nonstranded,gff_file => $gff_file);
is($annot_analyzer_nonstranded->isStranded(),0, 'stranded (1)');
is($chim_struct_nonstranded->nbChimeras(),2,'nbChimeras (1)');

#specific test for stranded protocol
my $chim_struct = CracTools::ChimCT::Structure->new();
$chim_struct->addChimera($chim); #class 3
$chim_struct->addChimera($chim2); #class 2 sense_(protein_coding)---sense_(non_coding)
is($chim_struct->nbChimeras(),2,'nbChimeras (2)');
my $annot_analyzer = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct,gff_file => $gff_file, is_stranded => 1);
is($annot_analyzer->isStranded(),1, 'stranded (2)');
is($chim_struct->nbChimeras(),2,'nbChimeras (3)');
my @headers = $annot_analyzer->getHeaders;
my @output = $annot_analyzer->getOutput($chim);
is($annot_analyzer->getFusionGenes($chim), 'TOTO---TOTO','getFusionGenes');
is($annot_analyzer->getFusionGenes($chim2), 'TOTO2---TITI','getFusionGenes (2)');
my ($score) = $annot_analyzer->getScore($chim);
is($score, 100,'getScore');
is($annot_analyzer->getExonsDistance($chim), "37---1", 'getExonsDistance (1)');
is($annot_analyzer->getExonsDistance($chim2), "4---16", 'getExonsDistance (2)');

#specific test for EST annotation
my $chim_struct_est = CracTools::ChimCT::Structure->new();
$chim_struct_est->addChimera($chim3);
$chim_struct_est->addChimera($chim4); #2/ a new splice (sense---NA) with a distance < SPLICE_EVENT_THRESHOLD 
$chim_struct_est->addChimera($chim5); #1/ a splice (Class2, sense, same gene)
$chim_struct_est->addChimera($chim6); 
is($chim_struct_est->nbChimeras(),4,'nbChimeras (4)');
my $annot_analyzer_est = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_est,gff_file => $gff_file,est_file => $est_file);
is($annot_analyzer_est->getFusionGenes($chim6), 'N/A---N/A','getOutput (3)');
is($chim_struct_est->nbChimeras(),2,'nbChimeras (5)');
## TODO: be careful, this following test is dangerous :) (that's remove $annot_analyzer_est but i don't know why)
#is($annot_analyzer->getExonsDistance($chim3), "7---20", 'getExonsDistance (3)');
is($annot_analyzer_est->getExonsDistance($chim3), "7---20", 'getExonsDistance (3)');
is($annot_analyzer_est->getExonsDistance($chim6), "N/A---58999999", 'getExonsDistance (5)');


# bug 17618 (submitted by T. Guignard) 
my $chim_struct_17618 = CracTools::ChimCT::Structure->new();
$chim_struct_17618->addChimera($chim_17618);
is($chim_struct_17618->nbChimeras(),1,'nbChimeras (6)');
my $annot_17618 = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_17618,gff_file => $gff_file, is_stranded => 1);
is($chim_struct_17618->nbChimeras(),0,'nbChimeras (7)');

# bug 17623 (submitted by T. Guignard) 
my $chim_struct_17623 = CracTools::ChimCT::Structure->new();
$chim_struct_17623->addChimera($chim_17623);
my $annot_17623 = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_17623,gff_file => $gff_file,est_file => $est_file, is_stranded => 1);
is($annot_17623->getExonsDistance($chim_17623),'1---2','getExonsDistance (6)');

# bug 17647 (submitted by T. Guignard) 
my $chim_struct_17647 = CracTools::ChimCT::Structure->new();
$chim_struct_17647->addChimera($chim_17647);
my $annot_17647 = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_17647,gff_file => $gff_file,est_file => $est_file, is_stranded => 1);
my @output_17647 = $annot_17647->getOutput($chim_17647);
is($annot_17647->getExonsDistance($chim_17647),'N/A---2','getExonsDistance (7)');

# bug 17676 (submitted by T. Guignard) 
my $chim_struct_17676 = CracTools::ChimCT::Structure->new();
$chim_struct_17676->addChimera($chim_17676);
my $annot_17676 = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_17676,gff_file => $gff_file,est_file => $est_file, is_stranded => 1);
is($chim_struct_17676->nbChimeras(),1,'nbChimeras (8)');
my @output_17676 = $annot_17676->getOutput($chim_17676);
ok($output_17676[1] =~ /INTRON---EXON/,'getOutput (4)');

# bug17994 (submitted by J. Audoux)
my $chim_struct_17994 = CracTools::ChimCT::Structure->new();
$chim_struct_17994->addChimera($chim_17994);
my $annot_17994 = CracTools::ChimCT::Analyzer::Annotation->new(chimera_struct => $chim_struct_17994,gff_file => $gff_file, is_stranded => 1);
my ($score_17994,$comments_17994) = $annot_17994->getScore($chim_17994);
is($score_17994,0, 'getScore() with pseudogenes');

__GFF__
1	Ensembl_CORE	exon	55	102	.	+	.	ID=ENSE00002706393;Parent=ENST00000578939
1	Ensembl_CORE	cds	55	102	.	+	.	ID=ENST00000578939.cds;Parent=ENST00000578939
1	Ensembl_CORE	mRNA	55	102	.	+	.	ID=ENST00000578939;Parent=ENSG00000266142;Exons_NB=1;type=protein_coding
1	Ensembl_CORE	exon	60	102	.	+	.	ID=ENSE00002706394;Parent=ENST00000578940
1	Ensembl_CORE	mRNA	60	102	.	+	.	ID=ENST00000578940;Parent=ENSG00000266142;Exons_NB=1;type=protein_coding
1	Ensembl_CORE	gene	55	102	.	+	.	ID=ENSG00000266142;Name=TOTO;Transcripts_NB=1
1	Ensembl_CORE	exon	55	102	.	-	.	ID=ENSE00002706493;Parent=ENST00000579939
1	Ensembl_CORE	mRNA	55	102	.	-	.	ID=ENST00000579939;Parent=ENSG00000269142;Exons_NB=1;type=protein_coding
1	Ensembl_CORE	gene	55	102	.	-	.	ID=ENSG00000269142;Name=TOTO2;Transcripts_NB=1
1	Ensembl_CORE	exon	1	20	.	-	.	ID=ENSE0000270622;Parent=ENST000005788
1	Ensembl_CORE	mRNA	1	20	.	-	.	ID=ENST000005788;Parent=ENSG00000266141;Exons_NB=1;type=miRNA
1	Ensembl_CORE	gene	1	20	.	-	.	ID=ENSG00000266141;Name=TITI;Transcripts_NB=1
7	Ensembl_CORE	exon	98957168	98957361	.	+	.	ID=ENSE00003632159;Parent=ENST00000262942,ENST00000432884;exon_rank=6
7	Ensembl_CORE	exon	98951532	98951744	.	+	.	ID=ENSE00003627403;Parent=ENST00000262942,ENST00000432786,ENST00000432884;exon_rank=7
7	Ensembl_CORE	exon	98930948	98931040	.	+	.	ID=ENSE00003557955;Parent=ENST00000432884;exon_rank=9
7	Ensembl_CORE	exon	98933076	98933107	.	+	.	ID=ENSE00003489116;Parent=ENST00000432884;exon_rank=10
7	Ensembl_CORE	exon	98984308	98984412	.	+	.	ID=ENSE00003466137;Parent=ENST00000432884;exon_rank=11
7	Ensembl_CORE	exon	98961166	98961256	.	+	.	ID=ENSE00003628680;Parent=ENST00000262942,ENST00000432884;exon_rank=12
7	Ensembl_CORE	exon	98937481	98937632	.	+	.	ID=ENSE00003668136;Parent=ENST00000432884;exon_rank=15
7	Ensembl_CORE	exon	98935804	98935908	.	+	.	ID=ENSE00003521456;Parent=ENST00000432884;exon_rank=19
7	Ensembl_CORE	exon	98985662	98985787	.	+	.	ID=ENSE00001664745;Parent=ENST00000441989,ENST00000432884;exon_rank=20
7	Ensembl_CORE	exon	98955963	98956038	.	+	.	ID=ENSE00003546812;Parent=ENST00000262942,ENST00000432884;exon_rank=22
7	Ensembl_CORE	exon	98923521	98923627	.	+	.	ID=ENSE00001683310;Parent=ENST00000441989,ENST00000432884;exon_rank=23
7	Ensembl_CORE	exon	98983325	98983401	.	+	.	ID=ENSE00003496371;Parent=ENST00000432884;exon_rank=27
7	Ensembl_CORE	exon	98941916	98942138	.	+	.	ID=ENSE00003481841;Parent=ENST00000262942,ENST00000432786,ENST00000432884;exon_rank=30
7	Ensembl_CORE	exon	98946475	98946582	.	+	.	ID=ENSE00003668403;Parent=ENST00000262942,ENST00000432786,ENST00000432884;exon_rank=37
7	Ensembl_CORE	mRNA	98923521	98985787	.	+	.	ID=ENST00000441989;Parent=ENSG00000241685;exons_nb=14;type=protein_coding:protein_coding
7	Ensembl_CORE	mRNA	98923533	98963880	.	+	.	ID=ENST00000262942;Parent=ENSG00000241685;exons_nb=10;type=protein_coding:protein_coding
7	Ensembl_CORE	mRNA	98923550	98963849	.	+	.	ID=ENST00000432786;Parent=ENSG00000241685;exons_nb=10;type=protein_coding:protein_coding
7	Ensembl_CORE	mRNA	98951558	98957865	.	+	.	ID=ENST00000471960;Parent=ENSG00000241685;exons_nb=4;type=protein_coding:protein_coding
7	Ensembl_CORE	mRNA	98956289	98963837	.	+	.	ID=ENST00000463009;Parent=ENSG00000241685;exons_nb=4;type=protein_coding:protein_coding
7	Ensembl_CORE	mRNA	98961088	98963885	.	+	.	ID=ENST00000477240;Parent=ENSG00000241685;exons_nb=2;type=protein_coding:protein_coding
7	Ensembl_CORE	mRNA	98923521	98985787	.	+	.	ID=ENST00000432884;Parent=ENSG00000241685;exons_nb=14;type=protein_coding:protein_coding
7	Ensembl_CORE	gene	98923521	98985787	.	+	.	ID=ENSG00000241685;Name=ARPC1A;transcripts_nb=7;exons_nb=39
9	Ensembl_CORE	exon	117156637	117156695	.	-	.	ID=ENSE00001357841;Parent=ENST00000307564;exon_rank=33
9	Ensembl_CORE	mRNA	117096436	117156695	.	-	.	ID=ENST00000307564;Parent=ENSG00000106948;exons_nb=22;type=protein_coding:protein_coding
9	Ensembl_CORE	gene	117096436	117156695	.	-	.	ID=ENSG00000106948;Name=AKNA;transcripts_nb=9;exons_nb=44
10	Ensembl_CORE	exon	77818424	77818541	.	+	.	ID=ENSE00003524454;Parent=ENST00000593699;exon_rank=32
10	Ensembl_CORE	mRNA	77407267	78317135	.	+	.	ID=ENST00000593699;Parent=ENSG00000148655;exons_nb=8;type=protein_coding:protein_coding
10	Ensembl_CORE	gene	77360998	78319925	.	+	.	ID=ENSG00000148655;Name=C10orf11;transcripts_nb=14;exons_nb=39
17	Ensembl_CORE	exon	62777552	62777734	.	-	.	ID=ENSE00002709454;Parent=ENST00000578036;exon_rank=21
17	Ensembl_CORE	exon	62781340	62781485	.	-	.	ID=ENSE00002720759;Parent=ENST00000578036;exon_rank=16
17	Ensembl_CORE	mRNA	62775377	62833257	.	-	.	ID=ENST00000578036;Parent=ENSG00000214176;exons_nb=18;type=protein_coding:pseudogene
17	Ensembl_CORE	gene	62775377	62833272	.	-	.	ID=ENSG00000214176;Name=PLEKHM1P;transcripts_nb=7;exons_nb=35
6	Ensembl_CORE	gene	32546546	32557625	.	-	.	ID=ENSG00000196126;Name=HLA-DRB1;transcripts_nb=1;exons_nb=6
6	Ensembl_CORE	mRNA	32546546	32557625	.	-	.	ID=ENST00000360004;Parent=ENSG00000196126;exons_nb=6;type=protein_coding:protein_coding
6	Ensembl_CORE	gene	32520490	32527799	.	-	.	ID=ENSG00000229391;Name=HLA-DRB6;transcripts_nb=4;exons_nb=12
6	Ensembl_CORE	mRNA	32522507	32527799	.	-	.	ID=ENST00000440876;Parent=ENSG00000229391;exons_nb=3;type=protein_coding:pseudogene
__EST__
2	Ensembl_EST	exon	2	40	.	+	.	ID=ENSESTE0000270622;Parent=ENSESTT000005788
2	Ensembl_EST	mRNA	2	40	.	+	.	ID=ENSESTT000005788;Parent=ENSESTG00000266141;Exons_NB=1;type=miRNA
2	Ensembl_EST	gene	2	40	.	+	.	ID=ENSESTG00000266141;Name=NONCODING_EST1;Transcripts_NB=1
2	Ensembl_EST	exon	100	200	.	+	.	ID=ENSESTE0000270623;Parent=ENSESTT000005789
2	Ensembl_EST	mRNA	100	200	.	+	.	ID=ENSESTT000005789;Parent=ENSESTG00000266142;Exons_NB=1;type=protein_coding
2	Ensembl_EST	gene	100	200	.	+	.	ID=ENSESTG00000266142;Name=PROTEIN_EST2;Transcripts_NB=1
2	Ensembl_EST	exon	100000000	200000000	.	+	.	ID=ENSESTE00002706230;Parent=ENSESTT0000057890
2	Ensembl_EST	mRNA	100000000	200000000	.	+	.	ID=ENSESTT0000057890;Parent=ENSESTG000002661420;Exons_NB=1;type=protein_coding:proteing_coding
2	Ensembl_EST	gene	100000000	200000000	.	+	.	ID=ENSESTG000002661420;Name=PROTEIN_EST_LONG;Transcripts_NB=1
9	Ensembl_EST	exon	117160718	117160783	.	-	.	ID=ENSESTE00000266792;Parent=ENSESTT00000078915;exon_rank=9
9	Ensembl_EST	mRNA	117160718	117160783	.	-	.	ID=ENSESTT00000078915;Parent=ENSESTG00000031284;exons_nb=5;type=protein_coding:protein_coding
9	Ensembl_EST	gene	117160718	117160783	.	-	.	ID=ENSESTG00000031284;Name=TONTON;transcripts_nb=3;exons_nb=9
17	Ensembl_EST	five	62760581	62760865	.	-	.	ID=ENSESTT00000061176.five;Parent=ENSESTT00000061176
17	Ensembl_EST	exon	62760569	62760865	.	-	.	ID=ENSESTE00000206155;Parent=ENSESTT00000061176;exon_rank=5
17	Ensembl_EST	mRNA	62754360	62760865	.	-	.	ID=ENSESTT00000061176;Parent=ENSESTG00000024173;exons_nb=5;type=protein_coding:protein_coding
17	Ensembl_EST	gene	62754360	62760865	.	-	.	ID=ENSESTG00000024173;Name=N/A;transcripts_nb=1;exons_nb=5
