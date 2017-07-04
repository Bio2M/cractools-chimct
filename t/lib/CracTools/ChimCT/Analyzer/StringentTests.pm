package CracTools::ChimCT::Analyzer::StringentTests;
# ABSTRACT: Stringent tests on chimeras

use strict;
use warnings;

use Carp;
use CracTools::ChimCT::Const;
use parent 'CracTools::ChimCT::Analyzer';

sub new {
  my $class = shift;

  # Creating Annotation Analyzer using the generic analyzer
  my $self  = $class->SUPER::new(@_);

  my %args = @_;
  if(defined $args{min_support}) {
    $self->minSupport($args{min_support});
  } else {
    $self->minSupport($CracTools::ChimCT::Const::MIN_SUPPORT);
  }
  if(defined $args{max_support}) {
    $self->maxSupport($args{max_support});
  } else {
    $self->maxSupport($CracTools::ChimCT::Const::MAX_SUPPORT);
  }

  return $self;
}

sub getHeaders {
  return ();
}

sub getOutput {
  return ();
}

sub getScore {
  my $self = shift;
  my $chimera = shift;
  my $score = 100;
  my $comments;
  my $nb_reads_normalized = int($chimera->nbReads / $self->chimeraStruct->nbReadsImported * $CracTools::ChimCT::Const::NORMALIZATION_CONST + 0.5);

  if($nb_reads_normalized > $self->maxSupport) {
    $score = 0;
    $comments .= 'high_support;';
  }elsif($nb_reads_normalized < $self->minSupport) {
    $score = 0;
    $comments .= 'low_support;';
  }
  if ($chimera->isAnchored()) {
    #$score = 0;
    $comments .= 'chimera_anchored;';
  }
  if(defined $comments) {
    chop($comments);
    $comments = "StringentTests=$comments";
  }
  return ($score,$comments);
}

sub minSupport {
  my $self = shift;
  my $new_support = shift;
  if(defined $new_support) {
    $self->{min_support} = $new_support;
  }
  return $self->{min_support};
}

sub maxSupport {
  my $self = shift;
  my $new_support = shift;
  if(defined $new_support) {
    $self->{max_support} = $new_support;
  }
  return $self->{max_support};
}

1;
