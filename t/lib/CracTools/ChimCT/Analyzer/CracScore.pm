package CracTools::ChimCT::Analyzer::CracScore;
# ABSTRACT: Integrate the score gave to the chimera by CRAC as an Analyzer

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
  my $threshold = shift;
  my $best_read = $chimera->getBestRead();
  my ($score,$comments) = (0,undef);

  # That means we have at least one read with a score available
  # Class 2 chimeras that were recovered splice does not have a chimeraEvent score
  # so we give them a 1.
  if(defined $best_read && defined $best_read->chimeraEvent() && defined $best_read->chimeraEvent->{score}) {
    foreach my $read (@{$chimera->reads}) {
      $score+=$read->chimeraEvent->{score};  
    }
    $score = $score / @{$chimera->reads};
    $comments = "mean_crac_score=$score";
  } else {
    $score = 1;
  }

  return ($score,$comments);
  #return (100,undef);
}


1;
