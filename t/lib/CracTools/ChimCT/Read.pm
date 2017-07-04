package CracTools::ChimCT::Read;
# ABSTRACT: Representation of a chimeric read

use strict;
use warnings;

use POSIX;

use CracTools::Utils;

=head2 new

  Arg [id]            : Integer - read id
  Arg [pos_junction]  : Integer - Position in the read where the break start
  Arg [seq]           : String - Sequence of the read

  Description : Create a new read object
  ReturnType  : CracTools::ChimCT::Read
  Exceptions  : none

=cut

sub new {
  my $class = shift;

  my %args = @_;

  my $self = bless{ 
      id => $args{id},
      pos_junction => $args{pos_junction},
      seq => $args{seq},
  }, $class;

  $self->samLine($args{sam_line});
  $self->chimeraEvent($args{chimera_event});

  return $self;
}

sub seq {
  my $self = shift;
  my $new_seq = shift;
  if(defined $new_seq) {
    $self->{seq} = $new_seq;
  }
  return $self->{seq};
}

sub seqLength {
  my $self = shift;
  return length $self->seq;
}

sub id {
  my $self = shift;
  my $new_id = shift;
  if(defined $new_id) {
    $self->{id} = $new_id;
  }
  return $self->{id};
}

sub posJunction {
  my $self = shift;
  my $new_pos_junction = shift;
  if(defined $new_pos_junction) {
    $self->{pos_junction} = $new_pos_junction;
  }
  return $self->{pos_junction};
}

sub samLine {
  my $self = shift;
  my $sam_line = shift;
  if(defined $sam_line) {
    $self->{sam_line} = $sam_line;
  }
  return $self->{sam_line};
}

sub chimeraEvent {
  my $self = shift;
  my $chimera_event = shift;
  if(defined $chimera_event) {
    $self->{chimera_event} = $chimera_event;
  }
  return $self->{chimera_event};
}

sub isReversed {
  my $self = shift;
  my $reverse = shift;
  if(defined $reverse) {
    $self->{reversed} = $reverse;
  } elsif (!defined $self->{reversed}) {
    $self->{reversed} = 0;
  }
  return $self->{reversed};
}

sub pSupport {
  my $self = shift;
  if(defined $self->samLine) {
    my $p_support = $self->samLine->pSupport;
    if(defined $p_support && $self->isReversed) {
      $p_support = CracTools::Utils::reverse_tab($p_support);
    }
    return $p_support;
  } else {
    return undef;
  }
}

sub pLoc {
  my $self = shift;
  if(defined $self->samLine) {
    my $p_loc = $self->samLine->pLoc;
    if(defined $p_loc && $self->isReversed) {
      $p_loc = CracTools::Utils::reverse_tab($p_loc);
    }
    return $p_loc;
  } else {
    return undef;
  }
}

=head2 reverseRead

  Example : $read->reverseRead();
  Description : Reverse complemente the read sequence as well as the support profil
                and the location profil
  ReturnType  : none
  Exceptions  : none

=cut

sub reverseRead {
  my $self = shift;
  $self->posJunction(length($self->seq) - $self->posJunction - 1);
  $self->seq(CracTools::Utils::reverseComplement($self->seq));
  if($self->isReversed) {
    $self->isReversed(0);
  } else {
    $self->isReversed(1);
  }
}

1;
