package CracTools::ChimCT::Analyzer::Read;
# ABSTRACT: Print optional field in chimCT output for the best junction read

use strict;
use warnings;

use Carp;
use parent 'CracTools::ChimCT::Analyzer';

sub new {
  my $class = shift;

  # Creating Annotation Analyzer using the generic analyzer
  my $self  = $class->SUPER::new(@_);

  my %args = @_;
  
  if(defined $args{detailed_sam} && $args{detailed_sam} == 1) {
    $self->{detailed_sam} = 1;
  } else {
    $self->{detailed_sam} = 0;
  }

  return $self;
}

#sub getHeaders {
#  my $self = shift;
#  my @headers = ('Read id','Pos junction','Seq');
#  if($self->detailedSam) {
#    push (@headers,'P_support','P_loc');
#  }
#  return @headers;
#}

sub getOutput {
  my $self = shift;
  my $chim = shift;
  my $read = $chim->getBestRead();
  my @output = ();
  push(@output, "Read_id=".$read->id) unless !defined $read->id;
  push(@output, "Pos_junction=".$read->posJunction) unless !defined $read->posJunction;
  push(@output, "Read_seq=".$read->seq) unless !defined $read->seq;
  if($self->detailedSam) {
    my ($p_support,$p_loc) = ($read->pSupport,$read->pLoc);
    if(!defined $p_support) {
      $p_support = 'N/A';
    }
    if(!defined $p_loc) {
      $p_loc = 'N/A';
    }
    push (@output,"Profil_support=".$p_support,"Profil_location=".$p_loc);
  }
  return @output;
}

sub detailedSam {
  my $self = shift;
  return $self->{detailed_sam};
}

1;
