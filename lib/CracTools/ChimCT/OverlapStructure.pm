package CracTools::ChimCT::OverlapStructure;
# ABSTRACT: Structure to store chimeras and search for overlaps

=head1 SYNOPSIS

This structure is designed to store chimeras with indexes in order to search for overlaps.
This aim of OverlapStructure is to store paired-end chimera and cross them with spanning
chimeras where we know the exact position of the junction.

=cut

use strict;
use warnings;

use Set::IntervalTree;
use Carp;

use CracTools::ChimCT::Const;
use CracTools::ChimCT::Chimera;

=head2 new

  Exemple     : $chimeras = CracTools::ChimCT::Structure->new();
  Description : Create a new CracTools::ChimCT::Structure object
  ReturnType  : CracTools::ChimCT::Structure
  Exceptions  : none

=cut

sub new {
  my $class = shift;
  my %args = @_;

  my $self = bless{
    indexes => {}, # In this hash we store data structure for quick access to chimeras
    is_stranded => $args{is_stranded},
    nb_chimeras => 0,
  }, $class;

  if(defined $args{max_overlapping_distance}) {
    $self->{max_overlapping_distance} = $args{max_overlapping_distance};
  } else {
    $self->{max_overlapping_distance} = $CracTools::ChimCT::Const::PAIRED_MAX_OVERLAPPING_DISTANCE;
  }

  return $self;
}

sub isStranded {
  my $self = shift;
  my $is_stranded = shift;
  if(defined $is_stranded) {
    $self->{is_stranded} = $is_stranded;
  }
  if(defined $self->{is_stranded}) {
    return $self->{is_stranded};
  } else {
    return 0;
  }
}
  
sub nbChimeras {
  my $self = shift;
  return $self->{nb_chimeras};
}

sub maxOverlappingDistance {
  my $self = shift;
  my $dist = shift;
  if(defined $dist) {
    $self->{max_overlapping_distance} = $dist;
  }
  return $self->{max_overlapping_distance};
}

sub isValidOverlappingDistance {
  my $self = shift;
  my ($chim1,$chim2) = @_;
  return overlapDistance($chim1,$chim2) <= $self->maxOverlappingDistance();
}

sub addChimera {
  my $self = shift;
  my $chim = shift;
  my $class = $chim->getClass();

  $self->{nb_chimeras}++;

  # CLASS 1 CHIMERAS, we use arrays as index
  if($class == 1) {
    my $key1 = $chim->chr1."@".$chim->strand1;
    my $key2 = $chim->chr2."@".$chim->strand2;
    if(!defined $self->{indexes}{$class}{$key1}{$key2}) {
      $self->{indexes}{$class}{$key1}{$key2} = [];
    }
    push(@{$self->{indexes}{$class}{$key1}{$key2}},$chim);
  # CLASS 2,3,4 CHIMERAS, we use intervalTrees as index
  } else {
    my ($pos_min,$pos_max) = getChimeraPosMinMax($chim);
    
    if(!defined $self->{indexes}{$class}{$chim->chr1}{$chim->strand1}) {
      $self->{indexes}{$class}{$chim->chr1}{$chim->strand1} = Set::IntervalTree->new;
    }
    $pos_max++ if $pos_min == $pos_max; # Avoid null width intervals
    $self->{indexes}{$class}{$chim->chr1}{$chim->strand1}->insert($chim,$pos_min,$pos_max);
  }
}

sub getOverlappingChimeras {
  my $self = shift;
  my $chim = shift;
  my $class = $chim->getClass();
  my @overlapping_chimeras = ();

  # CLASS 1 CHIMERAS
  if($class == 1) {
    my $key1 = $chim->chr1."@".$chim->strand1;
    my $key2 = $chim->chr2."@".$chim->strand2;
    foreach my $candidate_chimera (@{$self->{indexes}{$class}{$key1}{$key2}}) {
      if($self->isValidOverlappingDistance($chim,$candidate_chimera)) {
        push(@overlapping_chimeras,$candidate_chimera);
      }
    }

    if(!$self->isStranded) {
      $key1 = $chim->chr2."@".-$chim->strand2;
      $key2 = $chim->chr1."@".-$chim->strand1;
      $chim->reverseChimera();
      foreach my $candidate_chimera (@{$self->{indexes}{$class}{$key1}{$key2}}) {
        #$candidate_chimera->reverseChimera(); # Reverse chimera coordinates to compute overlapping distance
        if($self->isValidOverlappingDistance($chim,$candidate_chimera)) {
          #$candidate_chimera->reverseChimera(); # Reverse back chimera coordinates
          push(@overlapping_chimeras,$candidate_chimera);
        }
      }
      $chim->reverseChimera();
    }
  # CLASS 2,3,4 CHIMERAS
  } else {
    my ($pos_min,$pos_max) = getChimeraPosMinMax($chim);
    $pos_max++ if $pos_min == $pos_max; # Avoid null width intervals

    my $hits;
    if(defined $self->{indexes}{$class}{$chim->chr1}{$chim->strand1}) {
      $hits = $self->{indexes}{$class}{$chim->chr1}{$chim->strand1}->fetch($pos_min,$pos_max);

      foreach my $candidate_chimera (@{$hits}) {
        if($self->isValidOverlappingDistance($chim,$candidate_chimera)) {
          push(@overlapping_chimeras,$candidate_chimera);
        }
      }
    }

    if(!$self->isStranded && defined $self->{indexes}{$class}{$chim->chr2}{-$chim->strand2}) {
      $hits = $self->{indexes}{$class}{$chim->chr2}{-$chim->strand2}->fetch($pos_min,$pos_max);
      $chim->reverseChimera();
      foreach my $candidate_chimera (@{$hits}) {
        #$candidate_chimera->reverseChimera(); # Reverse chimera coordinates to compute overlapping distance
        if($self->isValidOverlappingDistance($chim,$candidate_chimera)) {
          #$candidate_chimera->reverseChimera(); # Reverse back chimera coordinates
          push(@overlapping_chimeras,$candidate_chimera);
        }
      }
      $chim->reverseChimera();
    }
  }

  return \@overlapping_chimeras;
}

sub nbOverlappingChimeras {
  my $self = shift;
  my $chim = shift;
  return @{$self->getOverlappingChimeras($chim)};
}

# STATIC METHODS
#
sub getChimeraPosMinMax {
  my $chim = shift;
  my ($pos_min,$pos_max) = ($chim->pos1,$chim->pos2);
  ($pos_min,$pos_max) = ($pos_max,$pos_min) if $pos_max < $pos_min;
  return ($pos_min,$pos_max);
}

sub overlapDistance {
  my ($chim1,$chim2) = @_;
  return abs($chim1->pos1 - $chim2->pos1) + abs($chim1->pos2 - $chim2->pos2);
}

1;
