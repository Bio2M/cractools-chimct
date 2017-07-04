package CracTools::ChimCT::Analyzer::ChimScore;
# ABSTRACT: Integrate the score gave to the chimera by Random Forest regression

use strict;
use warnings;

use Carp;
use Statistics::Basic qw(:all);
use CracTools::ChimCT::Const;
use CracTools::ChimCT::Analyzer::Annotation;

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

sub getChimScore {
  my $self = shift;
  my $chimera = shift;
  my ($chimScore,$comments) = (1,undef);
  my $best_read = $chimera->getBestRead();
  my @p_support = split(",",$best_read->pSupport);
  my $spanning_junction_count=$chimera->nbReads;
  my $support_profile_stddev = stddev(@p_support);
  my $cv = $support_profile_stddev / $spanning_junction_count unless $spanning_junction_count==0 && $support_profile_stddev==0;
  $chimScore = 1-(1/(1+exp(-(0.00984317624465517 - 0.00318469422109822*$cv + 0.00984317624465517*$spanning_junction_count))));
  $comments = "ChimScore=$chimScore"; 
  return ($comments,$chimScore);
}

1;
