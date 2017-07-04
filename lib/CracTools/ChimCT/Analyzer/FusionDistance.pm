package CracTools::ChimCT::Analyzer::FusionDistance;
# ABSTRACT: Compute a fusion distance "score" for class 3 chimeras

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

  $self->fusionDistanceThreshold($args{fusion_distance_threshold});

  return $self;
}

sub fusionDistanceThreshold {
  my $self = shift;
  my $new_dist = shift;
  if(defined $new_dist) {
    $self->{fusion_distance_threshold} = $new_dist;
  } elsif(!defined $self->{fusion_distance_threshold}) {
    $self->{fusion_distance_threshold} = $CracTools::ChimCT::Const::FUSION_DISTANCE_THRESHOLD;
  }
  return $self->{fusion_distance_threshold};
}

sub getScore {
  my $self = shift;
  my $chimera = shift;
  croak("Missing chimera argument") unless defined $chimera; 

  # this analyzer is only available for class 3 chimeras
  # If chimera is not a class 3 we return MAX_SCORE (100).
  if($chimera->getClass ne 3) {
    return 100;
  } else {
    my ($score,$comments);
    my $read = $chimera->getBestRead;
    my ($d1,$d2) = ($chimera->pos1,$chimera->pos2);

    # try to update positions to search for chimeric overlap
    if(defined $read) {
      if($chimera->strand1 eq 1) {
        $d1 = $d1 - $read->posJunction;
        $d2 = $d2 + (length($read->seq) - $read->posJunction)
      } else {
        $d1 = $d1 + $read->posJunction;
        $d2 = $d2 - (length($read->seq) - $read->posJunction)
      }
    }

    # Compute distance between both side of the chimera
    my $fusion_distance = abs($d1 - $d2);
    # In this case there is an overlap
    if(($chimera->strand1 eq 1 && $d2 > $d1) || ($chimera->strand1 eq -1 && $d1 > $d2)) {
      $score = 0;
      $comments = 'FusionDistance=overlap';
    } elsif($fusion_distance > $self->fusionDistanceThreshold) {
      $score = 100;
    } else {
      $score = $fusion_distance / $self->fusionDistanceThreshold * 100;
      $comments = "FusionDistance=short_distance($fusion_distance)";
    }
    return ($score,$comments);
  }
}

1;
