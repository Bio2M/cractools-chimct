package CracTools::ChimCT::Analyzer;
# ABSTRACT: Base class for CracTools::ChimCT Analyzers

use strict;
use warnings;

use Carp;

sub new {
  my $class = shift;
  my %args = @_;

  croak "Missing chimera structure in argument" unless defined $args{chimera_struct};

  my $self = bless {
    chimera_struct => $args{chimera_struct},
  },$class;

  return $self;
}

sub getHeaders {
  #croak "This method should be overriden and never called";
  return ();
}

sub getOutput {
  #croak "This method should be overriden and never called";
  return ();
}

# Return a score and eventually a comment about this score
sub getScore {
  #croak "This method should be overriden and never called";
  return 100;
}

# Return a binary classification and eventually a comment
sub getChimScore {
  #croak "This method should be overriden and never called";
  return 0;
}

sub chimeraStruct {
  my $self = shift;
  return $self->{chimera_struct}
}

1;
