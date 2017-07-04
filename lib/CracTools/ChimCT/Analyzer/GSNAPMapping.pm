package CracTools::ChimCT::Analyzer::GSNAPMapping;
# ABSTRACT: Use GSNAP to perform a second mapping on chimeric reads

use strict;
use warnings;

use Carp;
use CracTools::ChimCT::Const;
use CracTools::SAMReader::SAMline;
use File::Temp; #To open temporary files in the safe way
use parent 'CracTools::ChimCT::Analyzer';

sub new {
  my $class = shift;

  # Creating Annotation Analyzer using the generic analyzer
  my $self  = $class->SUPER::new(@_);

  my %args = @_;

  croak "Missing gsnap_exe argument" unless defined $args{gsnap_exe};
  croak "Missing gsnap_genome_directory argument" unless defined $args{gsnap_genome_dir};
  croak "Missing gsnap_genome_name argument" unless defined $args{gsnap_genome_name};

  $self->keepGsnapOutput($args{keep_gsnap_output});
  $self->gsnapExe($args{gsnap_exe});
  $self->gsnapGenomeDir($args{gsnap_genome_dir});
  $self->gsnapGenomeName($args{gsnap_genome_name});
  $self->gsnapNbThreads($args{gsnap_nb_threads});
  $self->gsnapSoftclipThreshold($args{gsnap_softclip_threshold});
  $self->tmpDir($args{tmp_dir});

  $self->gsnapNbThreads($CracTools::ChimCT::Const::GSNAP_NB_THREADS) unless defined $self->gsnapNbThreads;
  $self->gsnapSoftclipThreshold($CracTools::ChimCT::Const::GSNAP_SOFTCLIP_THRESHOLD) unless defined $self->gsnapSoftclipThreshold;

  $self->_init;

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
  my $nb_read_mapped = $self->nbReadsMapped($chimera);
  my $score = 100 * (@{$chimera->reads} - $nb_read_mapped) / @{$chimera->reads};
  my $comments;
  if($score < 100) {
    $comments = "GSNAPMapping='$nb_read_mapped read mapped out of ".@{$chimera->reads}."'";
  }
  return ($score,$comments);
}

sub keepGsnapOutput {
  my $self = shift;
  my $keep_gsnap_output = shift;
  if(defined $keep_gsnap_output) {
    $self->{keep_gsnap_output} = $keep_gsnap_output;
  }
  return $self->{keep_gsnap_output};
}

sub gsnapExe {
  my $self = shift;
  my $gsnap_exe = shift;
  if(defined $gsnap_exe) {
    $self->{gsnap_exe} = $gsnap_exe;
  }
  return $self->{gsnap_exe};
}

sub gsnapGenomeDir {
  my $self = shift;
  my $gsnap_genome_dir = shift;
  if(defined $gsnap_genome_dir) {
    $self->{gsnap_genome_dir} = $gsnap_genome_dir;
  }
  return $self->{gsnap_genome_dir};
}

sub gsnapGenomeName {
  my $self = shift;
  my $gsnap_genome_name = shift;
  if(defined $gsnap_genome_name) {
    $self->{gsnap_genome_name} = $gsnap_genome_name;
  }
  return $self->{gsnap_genome_name};
}

sub gsnapNbThreads {
  my $self = shift;
  my $gsnap_nb_threads = shift;
  if(defined $gsnap_nb_threads) {
    $self->{gsnap_nb_threads} = $gsnap_nb_threads;
  }
  return $self->{gsnap_nb_threads};
}

sub gsnapSoftclipThreshold {
  my $self = shift;
  my $gsnap_softclip_threshold = shift;
  if(defined $gsnap_softclip_threshold) {
    $self->{gsnap_softclip_threshold} = $gsnap_softclip_threshold;
  }
  return $self->{gsnap_softclip_threshold};
}

sub nbReadsMapped {
  my ($self,$chimera) = @_;
  if(defined $self->{mapped}{$chimera->getKey}) {
    return $self->{mapped}{$chimera->getKey};
  } else {
    return 0;
  }
}

sub tmpDir {
  my ($self,$new_tmp_dir) = @_;
  $self->{tmp_dir} = $new_tmp_dir unless !defined $new_tmp_dir;
  return $self->{tmp_dir};
}

sub _addMappingValue {
  my ($self,$chimera_key) = @_;

  if(defined $self->{mapped}{$chimera_key}) {
    $self->{mapped}{$chimera_key}++;
  } else {
    $self->{mapped}{$chimera_key} = 1;
  }
}

sub _init {
  my $self = shift;

  # Temp fasta file for mapping read
  my $fasta_file = new File::Temp( SUFFIX => '.fasta', DIR => $self->tmpDir);#, UNLINK => 0);

  # Fill the fasta file with all read in the chimera structure
  # We store the chimeraKey and the readKey in the comment line of
  # each sequence in order to get back at them later
  foreach my $chim (values %{$self->chimeraStruct->chimeras}) {
    print $fasta_file ">",$chim->getKey,":",$_->id,"\n",$_->seq, "\n" foreach @{$chim->reads};
  }
  close $fasta_file;

  # Temp file for Gsnap output
  my $gsnap_output = new File::Temp(DIR => $self->tmpDir);#UNLINK => 0);

  if($self->keepGsnapOutput) {
    $gsnap_output->unlink_on_destroy(0);
    print STDERR "GSNAP outputs will be kept after analysis (see files with $gsnap_output prefix)\n";
  }

  # Temp file for GSNAP STDERR
  my $gsnap_STDERR = new File::Temp( SUFFIX => '.log', DIR => $self->tmpDir);

  # Run GSNAP on this file
  # TODO save error into a log file and provide a warning with a link to this file if there was an error
  # during the execution of GSNAP!
  my $command_line = $self->gsnapExe." $fasta_file -t ".$self->gsnapNbThreads." -n 1 -N 1 -A sam -D ".$self->gsnapGenomeDir." -d ".$self->gsnapGenomeName." --split-output=$gsnap_output 2> $gsnap_STDERR";

  # If everything went right we analyse GSNAP output files and verify chimeras
  if(system($command_line) == 0) {
    # Analysing GSNAP ouput
    my $mult_filename = "$gsnap_output\.unpaired_mult";
    my $uniq_filename = "$gsnap_output\.unpaired_uniq";

    $self->_verifyChimerasMapping($mult_filename) if (-e $mult_filename);
    $self->_verifyChimerasMapping($uniq_filename) if (-e $uniq_filename);
    
  # Else we print a warnings on STDERR and give the path to the log file
  } else {
    # The log file will not be automatically unlink so the user can take a look at it
    $gsnap_STDERR->unlink_on_destroy(0);
    print STDERR "Error in GSNAP execution (see GSNAP log file $gsnap_STDERR for more information)\n";
  }
}

sub _verifyChimerasMapping {
  my $self = shift;
  my $sam_file = shift;

  # we capture warnings because we are using sam_file that do not have a ".sam" extension
  # and we dont want CracTools::SAMReader to print uggly warnings
  local $SIG{__WARN__} = sub { };

  # Create a CracTools::SAMReader to easily read the file
  my $sam_reader = CracTools::SAMReader->new($sam_file);
  # Get an iterator for this file
  my $sam_it = $sam_reader->iterator();
  # loop on each alignement
  while (my $sam_line = $sam_it->()) {

    # Get the chim_id and the read_id from the sequence query name
    my ($chim_id,$read_id) = split(":",$sam_line->qname,2);

    # get the cigar string
    my $cigar = $sam_line->cigar;

    # get the chimera that correspond to this alignement, the query name (qname) is the chimera key
    my $chimera = $self->chimeraStruct->getChimera($chim_id);

    # get the read object corresponding to this sam line
    my $read;
    foreach (@{$chimera->reads}) {
      if ($_->id eq $read_id) {
        $read = $_;
        last;
      }
    }

    # compute the read length (once for all)
    my $read_length = length($read->seq);

    # get the position of the chimeric junction and reverse it (if needed) in order to have the same "orientation"
    # as the alignment give
    my $pos_junction = $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED})? $read_length - $read->posJunction - 1 : $read->posJunction;

    # get the value of the left and right softclip(s) (if any)
    my ($left_softclip) = $cigar =~ /^(\d+)[SH]/;
    my ($right_softclip) = $cigar =~ /(\d+)[SH]$/;

    # boolean that will be use to check if the chimera is mapped
    my $is_mapped = 1;

    # Check if the left part of the chimera is mapped
    # 1. if there is no softclip, we can consider this part as "mapped"
    # 2. if there is a softclip but there are enough mapped nucleotides around the left side
    #    of the junction. Enough is defined by CracTools::ChimCT::Const::GSNAP_SOFTCLIP_THRESHOLD constant
    # NB : +1 because pos_junction is 0 based
    if(defined $left_softclip && ($pos_junction - $left_softclip + 1) < $self->gsnapSoftclipThreshold) {
      $is_mapped = 0; 
    }
    
    # Check if the right part of the chimera is mapped
    # 1. if there is no softclip, we can consider this part as "mapped"
    # 2. if there is a softclip but there are enough mapped nucleotides around the right side
    #    of the junction. Enough is defined by CracTools::ChimCT::Const::GSNAP_SOFTCLIP_THRESHOLD constant
    # NB : +1 because pos_junction is 0 based
    if($is_mapped && defined $right_softclip &&  ($read_length - $pos_junction + 1 - $right_softclip) < $self->gsnapSoftclipThreshold) {
      $is_mapped = 0; 
    }

    # If read is mapped we marked it as is
    if($is_mapped) {
      $self->_addMappingValue($chim_id);
    }
  }
}

1;
