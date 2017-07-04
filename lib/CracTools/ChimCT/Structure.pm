package CracTools::ChimCT::Structure;
# ABSTRACT: Structure to store CracTools::ChimCT::Chimera objects

use strict;
use warnings;

use CracTools::ChimCT::Chimera;
use CracTools::ChimCT::Read;
use CracTools::SAMReader;
use CracTools::SAMReader::SAMline;

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
    chimeras => {},
    chimera_events => {},
    sam_lines => [],
  }, $class;

  # If there is a SAM file in args we load it
  if(defined $args{SAM}) {
    $self->importFromSAM($args{SAM},$args{stranded});
  }

  return $self;
}

=head2 importFromSAM

  Arg [1]     : String - The SAM/BAM file (this file can be gzipped).
  Arg [2]     : (Optional) Boolean - Stranded protocol

  Description : Populate the structure with chimeras extracted from the SAM file.
  ReturnType  : (Number of chimeras added to the structure, Number of reads in the file)
  Exceptions  : none

=cut

sub importFromSAM {

  # Get arguments
  my ($self,$file,$is_stranded) = @_; 

  if(!defined $is_stranded) {
    $is_stranded = 0;
  }

  # Create a SAMReader to parse the file looking for chimeras
  my $sam_reader = CracTools::SAMReader->new($file,'CRAC');
  #my $sam_it = $sam_reader->iterator();
  my $sam_it = $sam_reader->iteratorFile('IGNORE_HEADERS');
  my $read_id = 0;
  my $nb_chimeras_added_to_the_structure = 0;

  # Loop over each line of the SAM file
  while (my ($line,$line_number) = $sam_it->()) {

    # If there is a chimera in this line
    #if(defined $line->events('chimera')) {
    #  my @chimeras = @{$line->events('chimera')};
    if(CracTools::SAMReader::SAMline::hasEvent($line,'chimera')) {
      # Create SAMline object from the unparsed sam line
      my $sam_line = CracTools::SAMReader::SAMline->new($line);
      my @chimeras = @{$sam_line->events('chimera')};

      # Loop over each chimera of the current SAM line
      foreach my $chim_event (@chimeras) {

        # Get chimeras positions
        my ($chr1,$pos1,$strand1) = @{$chim_event->{loc1}}{'chr','pos','strand'};
        my ($chr2,$pos2,$strand2) = @{$chim_event->{loc2}}{'chr','pos','strand'};
        my $chimera = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
        my $seq = $sam_line->seq;

        my $read = CracTools::ChimCT::Read->new(id => $sam_line->qname,
                                                         seq => $sam_line->getOriginalSeq, 
                                                         pos_junction => $chim_event->{pos}, 
                                                         sam_line => $sam_line, 
                                                         chimera_event => $chim_event);

        if($is_stranded && $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{FIRST_SEGMENT})) {
          $read->reverseRead();  
          $chimera->reverseChimera();
        } elsif (!$is_stranded && defined $self->getChimera($chimera->getReverseKey)) {
          $read->reverseRead();  
          $chimera->reverseChimera();
        }

        my $existing_chimera = $self->getChimera($chimera->getKey());
        if(defined $existing_chimera) {
          $existing_chimera->addRead($read);
        } else {
          $chimera->addRead($read);
          $self->addChimera($chimera);
        }
      } 
      $sam_line->genericInfo('line_number',$line_number);
      #$self->addSamLine($sam_line);
    }
    $read_id++;
  }
  $self->nbReadsImported($read_id);
  return ($nb_chimeras_added_to_the_structure,$read_id);
}

sub addChimera {
  my $self = shift;
  my $new_chimera = shift;
  my $new_chimera_key = $new_chimera->getKey;

  # We use chimera position as the key of the chimeras hash
  if(defined $self->getChimera($new_chimera_key)) {
    print STDERR "Cannot add chimera $new_chimera_key to the structure, another chimera with the same coordinates is already register\n";
  } else {
    $self->{chimeras}{$new_chimera_key} = $new_chimera;
  }
}

#sub addSamLine {
#  my $self = shift;
#  my $sam_line = shift;
#  if(defined $sam_line) {
#    push($self->{sam_lines},$sam_line);
#  } else {
#    croak("Missing argument");
#  }
#}

=head2 removeChimera

  Arg [1] : CracTools::ChimCT::Chimera object
  Arg [2] : (optional) Cause :
            - 'splice' : remove the chimera and classified it as splice

=cut

sub removeChimera {
  my $self = shift,
  my $chimera = shift;
  my ($cause,$info) = @_;
  croak('Missing chimera argument') unless defined $chimera;
  if(defined $self->getChimera($chimera->getKey)) {
    # TODO create a regex that would re-classify this chimera in the SAM file
    #if(defined $cause && $cause =~ /splice/) {
    #  foreach my $read (@{$chimera->reads}) {
    #    my $chim_event = $read->chimeraEvent;
    #    if(defined $chim_event) {
    #      # TODO 
    #      # - Use original 'type'
    #      # - Add $info to event
    #      $read->samLine->updateEvent($chim_event,'Junction',(type => 'normal',
    #                       pos  =>  $chim_event->{pos},
    #                       loc  =>  $chim_event->{loc1},
    #                       gap  =>  abs($chim_event->{loc1}{pos} - $chim_event->{loc2}{pos})));
    #    }
    #  }
    #}
    #print STDERR "Delete chimera (cause : $cause)\n";
    delete $self->{chimeras}{$chimera->getKey};
  } else {
    carp("Try do delete unknown chimera");
  }
}

sub getChimerasOrdered {
  my $self = shift;
  my @sorted_by_score = sort {
    $self->getChimera($b)->score <=> $self->getChimera($a)->score ||
    $self->getChimera($a)->getClass <=> $self->getChimera($b)->getClass ||
    $self->getChimera($b)->nbReads <=> $self->getChimera($a)->nbReads
  } keys %{$self->chimeras};
  #my @sorted_by_score = sort {$self->chimeras->{$b}->score <=> $self->chimeras->{$a}->score} keys $self->chimeras;
  my @chimeras = map {$self->chimeras->{$_}} @sorted_by_score;
  return \@chimeras;
}

sub getChimera {
  my $self = shift;
  my $key = shift;
  croak('Missing key argument') unless defined $key;
  return $self->{chimeras}{$key};
}

sub nbChimeras {
  my $self = shift;
  my %args = @_;
  my $chim_it = $self->chimerasIterator;
  my $nb_chimeras = 0;
  if(defined $args{class}) {
    my $class = $args{class};
    while (my $chim = $chim_it->()) {
      if($chim->getClass eq $class) {
        $nb_chimeras++;
      }
    }
  } elsif(defined $args{score}) {
    my $score = $args{score};
    while (my $chim = $chim_it->()) {
      if($chim->score eq $score) {
        $nb_chimeras++;
      }
    }
  } else {
    $nb_chimeras = scalar keys %{$self->chimeras};
  }
  return $nb_chimeras;
}

sub nbReads {
  my $self = shift;
  my $nb_reads = 0;
  my $chim_it = $self->chimerasIterator;
  while (my $chim = $chim_it->()) {
    $nb_reads += $chim->nbReads;
  }
  return $nb_reads;
}

# Number total of reads RNA-seq run used for normalization
sub nbReadsImported {
  my $self = shift;
  my $nb_reads = shift;
  if(defined $nb_reads) {
    $self->{nb_reads_imported} = $nb_reads;
  }
  return $self->{nb_reads_imported};
}

sub chimerasIterator {
  my $self = shift;
  my @keys = keys %{$self->chimeras};
  my $i = 0;
  return sub {
    return if $i == @keys;
    return $self->getChimera($keys[$i++]);
  }
}

sub getStats {
  my $self = shift;
  my $stats = "> Nb chimeras      : ".$self->nbChimeras."\n".
              "> Nb reads         : ".$self->nbReads."\n".
              "> Nb chimeras 100% : ".$self->nbChimeras(score => 100)."\n".
              "> Classification   :\n";
  for(my $class = 1; $class < 5; $class++) {
    $stats .= "\t> Class $class : ".$self->nbChimeras(class => $class)."\n";
  }
  return $stats;
}

#sub getPatch {
#  my $self = shift;
#  my $patch = '';
#  foreach my $sam_line (@{$self->{sam_lines}}) {
#    my $line_patch = $sam_line->getPatch($sam_line->genericInfo('line_number'));
#    if($line_patch) {
#      $patch .= $line_patch."\n";
#    }
#  }
#  return $patch;
#}

sub chimeras {
  my $self = shift;
  return $self->{chimeras};
}

1;
