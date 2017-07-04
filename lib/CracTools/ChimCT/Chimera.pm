=head1 SYNOPSIS

  my $chimera = CracTools::ChimCT::ChimeraStructure::Chimera->new('2',-1,12345,'2',1,12346);

  # Reverse chimera
  $chimera->reverseChimera();

  # Add a read that cover the chimera
  $chimera->addRead($read);

  # Print read data
  $chimera->print();

  # Check intersiting stuff
  $chimera->isAnchored();

=cut

package CracTools::ChimCT::Chimera;
# ABSTRACT: Object representing a chimera found by CRAC.

use strict;
use warnings;
use Carp;
use CracTools::ChimCT::Const;

use constant KEY_SEPARATOR => '@';

=head1 METHODS

=cut

=head2 new

  Arg [1]     : String - Chr1
  Arg [2]     : Integer - Pos1
  Arg [3]     : Integer - Strand1
  Arg [4]     : String - Chr2
  Arg [5]     : Integer - Pos2
  Arg [6]     : Integer - Strand2

  Exemple     : $read = CracTools::ChimCT::ChimeraStructure::Chimera->new(...);
  Description : Create a new chimera object
  ReturnType  : CracTools::ChimCT::ChimeraStructure::Chimera
  Exceptions  : none

=cut

sub new {
  my $class = shift;
  if (@_ < 6) {
    croak "Invalid number of arguments";
  }
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = @_;

  my $self = bless { 
      chr1 => $chr1,
      strand1 => $strand1,
      pos1 => $pos1,
      chr2 => $chr2,
      strand2 => $strand2,
      pos2 => $pos2,
      nb_reads => 0, # Total number of read added for the chimera
      is_anchored => 0, # True if nb_read > 1 and all read have the same "posJunction"
      reads => [], # An array that store the N first read added, whit N = $CracTools::ChimCT::Const::MAX_READS_BY_CHIMERA
  }, $class;

  return $self;
}

=head2 newFromKey

  Arg [1]     : String - Chr1
  Arg [2]     : Integer - Pos1
  Arg [3]     : Integer - Strand1
  Arg [4]     : String - Chr2
  Arg [5]     : Integer - Pos2
  Arg [6]     : Integer - Strand2

  Exemple     : $read = CracTools::ChimCT::ChimeraStructure::Chimera->new(...);
  Description : Create a new chimera object
  ReturnType  : CracTools::ChimCT::ChimeraStructure::Chimera
  Exceptions  : none

=cut

sub newFromKey {
  my $class = shift;
  my $key = shift;
  my ($chr1,$strand1,$pos1,$chr2,$strand2,$pos2) = split(KEY_SEPARATOR,$key);
  return $class->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
}

=head2 isSameChimera

  Arg [1]     : CracTools::ChimCT::Chimera

  Description : Return true if the chimera in argument is the same
                the chimera object.

=cut

sub isSameChimera {
  my $self = shift;
  my $ch = shift;
  return ($self->chr1 eq $ch->chr1 && $self->pos1 eq $ch->pos1 && $self->strand1 eq $ch->strand1 &&
          $self->chr2 eq $ch->chr2 && $self->pos2 eq $ch->pos2 && $self->strand2 eq $ch->strand2);
}

sub getChimericCoordinates {
  my $self = shift;
  return ($self->chr1,$self->pos1,$self->strand1,$self->chr2,$self->pos2,$self->strand2);
}

sub getReversedChimericCoordinates {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = $self->getChimericCoordinates;
  my $tmp = $chr1;
  $chr1 = $chr2;
  $chr2 = $tmp;
  $tmp = $pos1;
  $pos1 = $pos2;
  $pos2 = $tmp;
  $tmp = $strand2*-1;
  $strand2 = $strand1*-1;
  $strand1 = $tmp;
  return ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
}

=head2 reverseChimera

  Example     : $chimera->reverseChimera();
  Description : Reverse complement the chimera positions (before and after the junction).
  ReturnType  : none
  Exceptions  : none

=cut

sub reverseChimera {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = $self->getReversedChimericCoordinates;
  $self->chr1($chr1);
  $self->pos1($pos1);
  $self->strand1($strand1);
  $self->chr2($chr2);
  $self->pos2($pos2);
  $self->strand2($strand2);
}

=head2 addRead

  Arg [1]     : CracTools::SAMReader::SAMline $read
  Description : Add a read that cover the chimera. There is a maximum of $CracTools::ChimCT::Const::MAX_READS_BY_CHIMERA
                read saved. After this threshold only the count is increased.

=cut

sub addRead {
  my ($self,$new_read) = @_;
  if(!defined $new_read) {
    croak "Missing argument";
  }

  # 1/ increment the counter of reads
  $self->{nb_reads}++;

  # 2/ Check if the read has the same sequence as the previous one, in this case
  #    the chimera is considered as "Anchored"
  #    TODO Anchored should not be that simple! We do not want just sequences not to be equal
  #    but sequence that do not completly overlap (for ex: is there is a sequence error in one
  #    of the read covering the chimera, that does not mean it is not anchored...);
  $self->isAnchored(1) if $self->nbReads == 2;
  if($self->nbReads > 1 && $self->isAnchored() &&
     ($self->{last_added_read}->posJunction != $new_read->posJunction || $self->{last_added_read}->seqLength != $new_read->seqLength)) {
    $self->isAnchored(0);
    # We do not need to store that information anymore
    # this will save some memory...
    delete $self->{last_added_read};
  }

  # 3/ We add this new read to the list of reads if we have not reach the maximum number of reads
  # by chimeras
  if(@{$self->reads} < $CracTools::ChimCT::Const::MAX_READS_BY_CHIMERA) {
    push(@{$self->{reads}},$new_read);
  }

  # TODO Determine the Best_read dynamically in order to keep only one "sam_line" per chimera
  # This could be done right here ----------> .

  # We save this read as the last one in order to check if the chimera is anchored even if we have outpass the MAX_READS_BY_CHIMERA
  $self->{last_added_read} = $new_read if $self->isAnchored() || $self->nbReads < 2;
}

=head2 nbReads

  Example     : $chimera->nbReads();
  Description : Return the number of reads that cover the chimera.
                This method does a simple count of all the reads that have been
                added to the chimera using the addRead() method.
  ReturnType  : Integer
  Exceptions  : none

=cut

sub nbReads {
  my $self = shift;
  ##return scalar @{$self->{reads}};
  return $self->{nb_reads};
}

=head2 getBestRead

  Example     : $chimera->getBestRead();
  Description : Return a read that cover the chimera with the most centered
                junction position. 
                The read should not contain sequence error wherever possible.
  ReturnType  : CracTools::ChimCT::ChimeraStructure::Read or undef

=cut

sub getBestRead {
  my $self = shift;
  my $best_read;
  # Min distance betwen the chimera junction and the left/right size of the read.
  my $best_read_min_dist = -1;

  # If there is some reads that cover the chimera
  if(defined $self->reads) {
    foreach my $read (@{$self->reads}) {
      # If there is no seqError inside the read
      my $line = $read->samLine();
      if (!defined $line || scalar @{ $line->events('Error')} == 0){
        my $left_dist = $read->posJunction;
        if(defined $read->posJunction) {
          my $right_dist = (length($read->seq) - $left_dist);
          if ($left_dist > $best_read_min_dist && $right_dist > $best_read_min_dist) {
            $best_read = $read;
            if($left_dist > $right_dist) {
              $best_read_min_dist = $right_dist;
            } else {
              $best_read_min_dist = $left_dist;
            }
          }
        } elsif(!defined $best_read) {
          $best_read = $read;
        }
      }
    }
    # if all of chimeric reads contain at least one Error we get the first one
    if (!defined $best_read){
      $best_read = $self->reads->[0];
    }
  }
  return $best_read;
}

#sub getChimeraEventFromRead {
#  my $self = shift;
#  my $read = shift;
#  foreach my $chim_event (@{$read->events('chimera')}) {
#    my ($chr1,$pos1,$strand1) = @{$chim_event->{loc1}}{'chr','pos','strand'};
#    my ($chr2,$pos2,$strand2) = @{$chim_event->{loc2}}{'chr','pos','strand'};
#    if($self->isSameChimera($chr1,$pos1,$strand1,$chr2,$pos2,$strand2)) {
#      return $chim_event;
#    }
#  }
#  croak "The read ".$read->qname." doesn't contain the chimera that he was added to";
#}

=head2 isAnchored

  Example     : if($chimera->isAnchored()) {
                  doSomething();
                }
  Description : Return true if the chimera is Anchored (ie : every read that covers
                the chimera have the same sequence).
                However if there is only one read that covers the chimera "isAnchored()"
                will return false.
  ReturnType  : Integer
  Exceptions  : none

=cut

sub isAnchored {
  my $self = shift;
  my $is_anchored_value = shift;
  if (defined $is_anchored_value) {
    $self->{is_anchored} = $is_anchored_value;
  }
  return $self->{is_anchored};
}

=head2 getClass

  Example     : $chimera->getClass();
  Description : Return the class of the chimera (from 1 to 5) :
                - 1: If the chimera has a tipical profil, i.e locations on two different chr.
                - 2: If the two locations are both on same chr and strand but the distance is too large for a splice.
                     This class is inffered by CRAC classification. The threshold used to classify
                     this chimera was declared in CRAC configuration.
                - 3: If the two locations are both on same chr and strand but not in the good order. 
                - 4: If the two locations are on same chr but on different strand.  
  ReturnType  : Integer
  Exceptions  : none

=cut

sub getClass {
  my $self = shift;
  if ($self->chr1 eq $self->chr2){
    if ($self->strand1 == $self->strand2){
      if (($self->strand1 == 1 && $self->pos1 > $self->pos2) || ($self->strand1 == -1 && $self->pos2 > $self->pos1)){
        return 3;
      } else {
        return 2;
      }
    } else {
      return 4;
    }
  } else {
    return 1;
  }
}

=head2 getUniqId

  Return a Uniq identifier for the chimera.
  This identifier is the same for the reversed chimera.

=cut

sub getUniqId {
  my $self = shift;
  my $id = '';
  my $key1 = $self->getKey;
  my $key2 = $self->getReverseKey;
  if($key1 gt $key2) {
    return $key1;
  } else {
    return $key2;
  }
}

=head2 getKey

  Example     : $chimera->getKey();
  Description : Return the key to use for storing the chimera in a hash
  ReturnType  : String
  Exceptions  : none

=cut

sub getKey {
  my $self = shift;
  return $self->chr1.KEY_SEPARATOR.$self->strand1.KEY_SEPARATOR.$self->pos1.KEY_SEPARATOR.$self->chr2.KEY_SEPARATOR.$self->strand2.KEY_SEPARATOR.$self->pos2;
}

=head2 getReverseKey

  Example     : $key = $chimera->getReverseKey();
                $hash{$key} = $chimera;
  Description : Return the key to use for storing the reversed chimera in a hash.
                This method does not reverse the chimera, it just generate a key that
                correspond to the reversed chimera. As an illusatration :

                > $key1 = $chimera->getReverseKey();
                > $chimera->reverseChimera();
                > $key2 = $chimera->getKey();
                > ($key1 == $key2) # Will return true

  ReturnType  : String
  Exceptions  : none

=cut

sub getReverseKey {
  my $self = shift;
  my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = $self->getReversedChimericCoordinates;
  return $chr1.KEY_SEPARATOR.$strand1.KEY_SEPARATOR.$pos1.KEY_SEPARATOR.$chr2.KEY_SEPARATOR.$strand2.KEY_SEPARATOR.$pos2;
}

=head1 GETTERS AND SETTERS

=head2 score

=cut

sub score {
  my $self = shift;
  my $score = shift;
  if(defined $score) {
    $self->{score} = $score;
  }

  if(defined $self->{score}) {
    return $self->{score};
  } else {
    return 0;
  }
}

=head2 chr1

  Arg [1] : (optional) String - the chr of the chimera before the junction
  Example : $read->chr1(X);
  Description : Getter/setter for attribute chr1
  ReturnType  : String or undef
  Exceptions  : none

=cut

sub chr1 {
  my $self = shift;
  if(@_) {
    $self->{chr1} = shift;
  }
  return $self->{chr1};
}

=head2 strand1

  Arg [1] : (optional) Integer (1 or -1) - the strand of the chimera before the junction
  Example : $read->strand1(1);
  Description : Getter/setter for attribute strand1
  ReturnType  : Integer or undef
  Exceptions  : none

=cut

sub strand1 {
  my $self = shift;
  if(@_) {
    $self->{strand1} = shift;
  }
  return $self->{strand1};
}

=head2 pos1

  Arg [1] : (optional) Integer - the position of the chimera before the junction
  Example : $read->pos1(1234);
  Description : Getter/setter for attribute pos1
  ReturnType  : Integer or undef
  Exceptions  : none

=cut

sub pos1 {
  my $self = shift;
  if(@_) {
    $self->{pos1} = shift;
  }
  return $self->{pos1};
}

=head2 chr2

  Arg [1] : (optional) String - the chr of the chimera after the junction
  Example : $read->chr2(X);
  Description : Getter/setter for attribute chr2
  ReturnType  : String or undef
  Exceptions  : none

=cut

sub chr2 {
  my $self = shift;
  if(@_) {
    $self->{chr2} = shift;
  }
  return $self->{chr2};
}

=head2 strand2

  Arg [1] : (optional) Integer (1 or -1) - the strand of the chimera after the junction
  Example : $read->strand2(-1);
  Description : Getter/setter for attribute strand2
  ReturnType  : Integer or undef
  Exceptions  : none

=cut

sub strand2 {
  my $self = shift;
  if(@_) {
    $self->{strand2} = shift;
  }
  return $self->{strand2};
}

=head2 pos2

  Arg [1] : (optional) Integer - the position of the chimera after the junction
  Example : $read->pos2(1234);
  Description : Getter/setter for attribute pos2
  ReturnType  : Integer or undef
  Exceptions  : none

=cut

sub pos2 {
  my $self = shift;
  if(@_) {
    $self->{pos2} = shift;
  }
  return $self->{pos2};
}

=head2 reads

  Description : Getter for attribute reads
  ReturnType  : Array of CracTools::SAMReader::SAMline

=cut

sub reads {
  my $self = shift;
  return $self->{reads};
}

#=head2 read
#
#  Description : Getter for attribute read
#  ReturnType  : CracTools::SAMReader::SAMline
#
#=cut 
#
#sub read {
#  my $self = shift;
#  my $read = shift;
#  if (defined $read) {
#    $self->{read} = $read;
#  }
#  return $self->{read};
#}

1;
