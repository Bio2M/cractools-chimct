#! /usr/bin/perl

use strict;
use warnings;

use Test::More;
plan( skip_all => "Not yet integrated" );

# TODO Integrate test when parsing SAM files will be integrated as methods
# of CracTools/ChimCT.pm 

#use Test::More tests => 1;
#use CracTools::ChimCT::Analyzer::Annotation;
#use CracTools::ChimCT::Structure;
#use File::Temp 0.23;
#use Inline::Files 0.68;

# Create a temp file with the GFF lines described below
#my $sam_file = new File::Temp( SUFFIX => '.sam', UNLINK => 1);
#my $chim_struct = 
#while(<SAM>) {print $gff_file $_;}
#close $gff_file;

#__SAM__
#GA-C_0001:5:1:6619:21247#0/1    83      11      294172  254     26M23S  =       293903  -295    TGGTCTTCGGGTGCACGGGGTTCAGGGTCAGGCGCTGTGAACTTCCTGA       bbb`aa``bbb`bb`bbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbbb       IH:i:2  MC:Z:23M173N24M2S MQ:i:254        NH:i:1  NM:i:0  R2:Z:CCAACATGGCTGAACCCTTCAAGGTGTGGACGGAGAATGCAGACGGGGC  SA:Z:11,294196,-,26S23M,254,0   XB:Z:0:is_duplicated=0;genome_indels=-76;score_intra=-0.424704;score_inter=-0.738405;deviation=-0.666683;falling_left=0;falling_right=0;inside_first_quartile=2;inside_last_quartile=2;inside_score=2;outside_score=3.36364;average_low_inside=2;average_high_inside=2.33333;has_no_start_break=0;has_no_end_break=0;is_deviated=0;is_nice_break=0;is_very_nice_break=0;pos_start_break=7;pos_end_break=25      XC:i:1  XD:i:0  XE:Z:0:0:chimera:25:11|-1,294121:11|-1,294196:3:0.8     XM:i:0  XN:i:0  XO:Z:11|-1,294174       XP:Z:loc:1:0:0  XQ:i:26 XR:Z:p_support=5,5,5,5,5,5,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,3,3,3;p_loc=1,1,1,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,1,1,1,1    XU:i:1
