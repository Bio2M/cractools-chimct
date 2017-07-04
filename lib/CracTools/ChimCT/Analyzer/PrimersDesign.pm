package CracTools::ChimCT::Analyzer::PrimersDesign;
# ABSTRACT: Create a meta read from PE reads in order to easily design primers

use strict;
use warnings;

use Carp;
use CracTools::ChimCT::Const;
use CracTools::SAMReader;
use CracTools::SAMReader::SAMline;
use parent 'CracTools::ChimCT::Analyzer';

sub new {
  my $class = shift;

  # Creating Annotation Analyzer using the generic analyzer
  my $self  = $class->SUPER::new(@_);

  my %args = @_;

  croak("Missing sam_file in argument") unless defined $args{sam_file};

  $self->_init($args{sam_file});

  return $self;
}

sub getOutput {
  my $self = shift;
  my $chimera = shift;
  my $primer = $CracTools::ChimCT::Const::NOT_AVAILABLE;
  if(defined $chimera && defined $self->{primers}{$chimera->getKey()}) {
    $primer = $self->{primers}{$chimera->getKey()};
  }
  return "Primers=".$primer;
}

sub _init {
  my $self = shift;
  my $sam_file = shift;

  my %chim_reads = ();

  # First we record all chimeras read_id that belong to paired-end reads
  foreach my $chim (values %{$self->chimeraStruct->chimeras}) {
    my $read = $chim->getBestRead();
    if(defined $read) {
      # If this is a paired-end read we register the paired read_id
      my $paired_type;
      my $paired_id = $read->id;
      if($read->samLine->isFlagged($CracTools::SAMReader::SAMline::flags{FIRST_SEGMENT})) {
        $paired_type = 1;
      } else {
        $paired_type = 2;
      }

      # if PE reads have and extra /1 or /2 we remove it
      $paired_id = _removePairedIdentifier($paired_id);

      my $pos_junc = $read->posJunction;
      my ($seqStartBreak,$seqEndBreak) = $read->seq =~ /^([CATG]{$pos_junc})([CATG]+)$/;
      
      $chim_reads{$paired_id}{read} = $read;
      $chim_reads{$paired_id}{seq} = "$seqStartBreak\*$seqEndBreak";
      $chim_reads{$paired_id}{paired_type} = $paired_type;
      $chim_reads{$paired_id}{chim_key} = $chim->getKey();
      # seq => $seq, paired_type => $paired_type, chim_key => $chim->getKey());
    }
  }

  my $sam_reader = CracTools::SAMReader->new($sam_file);
  my $sam_it = $sam_reader->iteratorFile('IGNORE_HEADERS');
  while (my ($line,$line_number) = $sam_it->()) {
    my ($read_id) = $line =~ /^(\S+)/;
    # if PE reads have and extra /1 or /2 we remove it
    $read_id = _removePairedIdentifier($read_id);

    # TODO we could get the mate sequence using the "R2" field
    if(defined $chim_reads{$read_id}) {
      my $sam_line = CracTools::SAMReader::SAMline->new($line);
      
      # We want the primary alignement
      next if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{SECONDARY_ALIGNMENT}) || $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{CHIMERIC_ALIGNMENT});

      # We want the mate read not the one we already have
      next if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{FIRST_SEGMENT}) && $chim_reads{$read_id}{paired_type} eq 1;
      next if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{LAST_SEGMENT}) && $chim_reads{$read_id}{paired_type} eq 2;


      my($read1_seq,$read2_seq);
      my $primer = "";

      if($chim_reads{$read_id}{paired_type} eq 1) {
        if($chim_reads{$read_id}{read}->isReversed()) {
          $primer = $sam_line->getOriginalSeq."#".$chim_reads{$read_id}{seq};
        } else {
          $primer = $chim_reads{$read_id}{seq}."#".CracTools::Utils::reverseComplement($sam_line->getOriginalSeq);
        }
      } else {
        if($chim_reads{$read_id}{read}->isReversed()) {
          $primer = $sam_line->getOriginalSeq."#".$chim_reads{$read_id}{seq};
        } else {
          $primer = $chim_reads{$read_id}{seq}."#".CracTools::Utils::reverseComplement($sam_line->getOriginalSeq);
        }
      }
      $self->{primers}{$chim_reads{$read_id}{chim_key}} = $primer;
    }
  }
}

sub _removePairedIdentifier {
  my $read_id = shift;
  if($read_id =~ /\/1$|\/2$/) {
    $read_id =~ s/\/.$//;
  }
  return $read_id;
}

1;
