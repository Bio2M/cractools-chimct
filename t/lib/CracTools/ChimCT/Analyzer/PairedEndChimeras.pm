package CracTools::ChimCT::Analyzer::PairedEndChimeras;
# ABSTRACT: Find consistent PE chimeras

use strict;
use warnings;

use Carp;
use Fcntl qw( SEEK_SET );
use CracTools::SAMReader::SAMline;
use CracTools::ChimCT::Const;
use parent 'CracTools::ChimCT::Analyzer';

sub new {
  my $class = shift;

  # Creating Annotation Analyzer using the generic analyzer
  my $self  = $class->SUPER::new(@_);

  my %args = @_;

  croak "Missing paired_chimera_overlap_struct file in arguments" unless defined $args{paired_chimera_overlap_struct};
  $self->{paired_chimera_overlap_struct} = $args{paired_chimera_overlap_struct};

  return $self;
}

sub getHeaders {
  my $self = shift;
  return ('Spanning_PE');
}

sub getOutput {
  my $self = shift;
  my $chim = shift;
  return ($self->{nb_paired_end_chimera}{$chim->getKey});
}

sub getNbSpanningReads {
  my $self = shift;
  my $chim = shift;
  return $self->{paired_chimera_overlap_struct}->nbOverlappingChimeras($chim);
}

sub getNbPairedEndReads {
  my $self = shift;
  return $self->{paired_chimera_overlap_struct}->nbChimeras();
}

sub getScore {
  my $self = shift;
  my $chimera = shift;
  my ($score,$comments);
  croak("Missing chimera argument") unless defined $chimera;
  my $spanning_junction = $chimera->nbReads;
  my $spanning_paired = $self->getNbSpanningReads($chimera);
  if($spanning_paired >= $spanning_junction) {
    $score = 100;
  } else {
    $score = $spanning_paired / $spanning_junction * 100;
    $comments = "PairedEndChimeras=strange_paired_end_support";
  }
  return ($score,$comments);
}

1;
