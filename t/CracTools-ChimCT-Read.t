#! /usr/bin/perl

use Test::More tests => 5;
use CracTools::ChimCT::Read;

my $read = CracTools::ChimCT::Read->new(seq => 'ATGc', id => 'read1', pos_junction => 2);

# Testing Getters & Setters
is($read->seq,'ATGc','seq');
is($read->id,'read1','id');
is($read->posJunction,2,'posJunction');

# Revesing read and testing values
$read->reverseRead;
is($read->seq,'gCAT','reverseRead (1)');
# TODO reversing read, is it the right value for pos_junction???
is($read->posJunction,1,'reverseRead (2)');
