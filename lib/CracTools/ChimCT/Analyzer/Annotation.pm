###############################################################################
#                                                                             #
#    Copyright © 2012-2013 -- IRB/INSERM                                      #
#                            (Institut de Recherche en Biothérapie /          #
#                             Institut National de la Santé et de la          #
#                             Recherche Médicale)                             #
#                             LIRMM/UM2                                       #
#                            (Laboratoire d'Informatique, de Robotique et de  #
#                             Microélectronique de Montpellier /              #
#                             Université de Montpellier 2)                    #
#                                                                             #
#  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@inserm.fr>             #
#                   Jerome AUDOUX  <jerome.audoux@univ-montp2.fr>             #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  Ce fichier  fait partie  du Pipeline  de traitement  de données NGS de la  #
#  plateforme ATGC labélisée par le GiS IBiSA.                                #
#                                                                             #
#  Ce logiciel est régi  par la licence CeCILL  soumise au droit français et  #
#  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  #
#  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  #
#  la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur  #
#  le site "http://www.cecill.info".                                          #
#                                                                             #
#  En contrepartie de l'accessibilité au code source et des droits de copie,  #
#  de modification et de redistribution accordés par cette licence, il n'est  #
#  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  #
#  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  #
#  titulaire des droits patrimoniaux et les concédants successifs.            #
#                                                                             #
#  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  #
#  associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au  #
#  développement  et à la reproduction du  logiciel par  l'utilisateur étant  #
#  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  #
#  manipuler et qui le réserve donc à des développeurs et des professionnels  #
#  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  #
#  utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du  #
#  logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la  #
#  sécurité de leurs systêmes et ou de leurs données et,  plus généralement,  #
#  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         #
#                                                                             #
#  Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez  #
#  pris connaissance  de la licence CeCILL,  et que vous en avez accepté les  #
#  termes.                                                                    #
#                                                                             #
#  -------------------------------------------------------------------------  #
#                                                                             #
#  This File is part of the NGS data processing Pipeline of the ATGC          #
#  accredited by the IBiSA GiS.                                               #
#                                                                             #
#  This software is governed by the CeCILL license under French law and       #
#  abiding by the rules of distribution of free software. You can use,        #
#  modify and/ or redistribute the software under the terms of the CeCILL     #
#  license as circulated by CEA, CNRS and INRIA at the following URL          #
#  "http://www.cecill.info".                                                  #
#                                                                             #
#  As a counterpart to the access to the source code and rights to copy,      #
#  modify and redistribute granted by the license, users are provided only    #
#  with a limited warranty and the software's author, the holder of the       #
#  economic rights, and the successive licensors have only limited            #
#  liability.                                                                 #
#                                                                             #
#  In this respect, the user's attention is drawn to the risks associated     #
#  with loading, using, modifying and/or developing or reproducing the        #
#  software by the user in light of its specific status of free software,     #
#  that may mean that it is complicated to manipulate, and that also          #
#  therefore means that it is reserved for developers and experienced         #
#  professionals having in-depth computer knowledge. Users are therefore      #
#  encouraged to load and test the software's suitability as regards their    #
#  requirements in conditions enabling the security of their systems and/or   #
#  data to be ensured and, more generally, to use and operate it in the same  #
#  conditions as regards security.                                            #
#                                                                             #
#  The fact that you are presently reading this means that you have had       #
#  knowledge of the CeCILL license and that you accept its terms.             #
#                                                                             #
###############################################################################

=head1 DESCRIPTION

Provide annotation for chimeras analyzed by chimCT and compute an Annotation "score".

=cut

package CracTools::ChimCT::Analyzer::Annotation;
# ABSTRACT: Annotation Analyzer for ChimCT

use strict;
use warnings;

use CracTools::Annotator;
use CracTools::ChimCT::Const;
#use List::Util qw[min max];
use Carp;
use parent 'CracTools::ChimCT::Analyzer';

=head1 METHODS

=head2 new

  Arg [gff_file]    : String - GFF file to perform annotation with
  Arg [est_file]    : String (Optional) - GFF file with EST to find annotation if none
                      have been found in the 
  Arg [keep_ig]     : String (Default: 0) - Keep genes annotated IG_ or TR_
  Arg [is_stranded] : String (Default: 0) - $gff_file

  Example     : my $annotation = CracTools::GFF::Annotation->new(gff_file => $gff, is_stranded => 1);
  Description : Create a new CracTools::ChimCT::Analyzer::Annotation object
  ReturnType  : CracTools::ChimCT::Analyzer::Annotation
  Exceptions  : none

=cut

sub new {
  my $class = shift;

  # Creating Annotation Analyzer using the generic analyzer
  my $self  = $class->SUPER::new(@_);

  my %args = @_;

  croak "Missing gff_file in arguments" unless defined $args{gff_file};

  # Declare an empty hash to strore annotations
  $self->{annotations} = {};
  $self->keepIG($args{keep_ig});
  $self->isStranded($args{is_stranded});
  $self->{annotator} = CracTools::Annotator->new($args{gff_file});
  if (defined $args{est_file}){
    $self->estFile($args{est_file});
    $self->{annotator_est} = CracTools::Annotator->new($args{est_file});
  }

  $self->_init();
  return $self;
}


=head2 getExonsRank

Return exons rank of the chimeric junction separated with "---" characters

=cut

sub getExonsRank {
  my $self = shift;
  my $chimera = shift;
  croak "Missing argument chimera in getExonsRank" unless defined $chimera;
  
  my $NA = $CracTools::ChimCT::Const::NOT_AVAILABLE;
  my ($exon_left_rank,$exon_right_rank) = ($self->getAnnotationAttribute($chimera,'before','exon','exon_rank'),
                                           $self->getAnnotationAttribute($chimera,'after','exon','exon_rank'));


  $exon_left_rank = $NA unless defined $exon_left_rank;
  $exon_right_rank = $NA unless defined $exon_right_rank;
  return $exon_left_rank."---".$exon_right_rank;
}

=head2 getExonsDistance

Return the distance between the chimeric junction and the start/end associated exon separated with "---" characters

=cut

sub getExonsDistance {
  my $self = shift;
  my $chimera = shift;
  
  croak "Missing argument chimera in getExonsDistance" unless defined $chimera;

  my $exon_left = $self->getAnnotationFeature($chimera,'before','exon');
  my $exon_right = $self->getAnnotationFeature($chimera,'after','exon');
  
  my $NA = $CracTools::ChimCT::Const::NOT_AVAILABLE;
  my ($exon_left_dist,$exon_right_dist) = ($NA,$NA);
  
  ##FIRST part of the junction
  #it means that the reads is inside an exon
  #In case of strand_specific, we check if we are in the same strand
  if (defined $exon_left 
      && (!$self->isStranded() || $exon_left->strand == $chimera->strand1)){  
      my $dist;
      if ($chimera->strand1 == 1){
	$dist = abs($chimera->pos1 - $exon_left->end);
      }else{
	$dist = abs($chimera->pos1 - $exon_left->start);
      }
      if ($exon_left_dist eq $NA || $dist < $exon_left_dist){
	  $exon_left_dist = $dist;
      }
  }
  #it means that the reads is inside an intron or an intergenic region
  else{
      my $candidates;
      if ($chimera->strand1 == 1){
	  $candidates = $self->annotator->getAnnotationNearestDownCandidates($chimera->chr1,$chimera->pos1,$chimera->strand1);
	  # If there is still no annotation (even non-coding)
	  # we try with EST
	  if(defined $self->annotator_est && (scalar @$candidates == 0)){
	      $candidates = $self->annotator_est->getAnnotationNearestDownCandidates($chimera->chr1,$chimera->pos1,$chimera->strand1);
	  }
      }else{
	  $candidates = $self->annotator->getAnnotationNearestUpCandidates($chimera->chr1,$chimera->pos1,$chimera->strand1);
	  if(defined $self->annotator_est && (scalar @$candidates == 0)){
	      $candidates = $self->annotator_est->getAnnotationNearestUpCandidates($chimera->chr1,$chimera->pos1,$chimera->strand1);
	  }
      }
      foreach my $candidate (@{ $candidates}){
	  #it means that the reads is inside an exon
	  #In case of strand_specific, we check if we are in the same strand
	  if ( defined $candidate->{exon} 
	       && (!$self->isStranded() ||  $candidate->{exon}->strand == $chimera->strand1)){  
	      my $dist;
	      if ($chimera->strand1 == 1){
		  $dist= abs($chimera->pos1 - $candidate->{exon}->end);    
	      }else{
		  $dist = abs($chimera->pos1 - $candidate->{exon}->start);  
	      }
	      if ($exon_left_dist eq $NA || $dist < $exon_left_dist){
		  $exon_left_dist = $dist;
	      }
	  }
      }
  }

  ##SECOND part of the junction
  #it means that the reads is inside an exon
  #In case of strand_specific, we check if we are in the same strand
  if (defined $exon_right 
      && (!$self->isStranded() || $exon_right->strand == $chimera->strand2)){  
      my $dist;
      if ($chimera->strand2 == 1){
	  $dist = abs($chimera->pos2 - $exon_right->start);  
      }else{
	  $dist = abs($chimera->pos2 - $exon_right->end);  
      }
      if ($exon_right_dist eq $NA || $dist < $exon_right_dist){
	  $exon_right_dist = $dist;
      }
  }
  #it means that the reads is inside an intron or an intergenic region
  else{
      my $candidates;
      if ($chimera->strand2 == 1){
	  $candidates = $self->annotator->getAnnotationNearestUpCandidates($chimera->chr2,$chimera->pos2,$chimera->strand2);
	  if(defined $self->annotator_est && (scalar @$candidates == 0)){
	      $candidates = $self->annotator_est->getAnnotationNearestUpCandidates($chimera->chr2,$chimera->pos2,$chimera->strand2);
	  }
      }else{
	  $candidates = $self->annotator->getAnnotationNearestDownCandidates($chimera->chr2,$chimera->pos2,$chimera->strand2);
	  if(defined $self->annotator_est && (scalar @$candidates == 0)){
	      $candidates = $self->annotator_est->getAnnotationNearestDownCandidates($chimera->chr2,$chimera->pos2,$chimera->strand2);
	  }
      }
      foreach my $candidate (@{ $candidates}){
	  #In case of strand_specific, we check if we are in the same strand
	  if (defined $candidate->{exon} 
	      && (!$self->isStranded() || $candidate->{exon}->strand == $chimera->strand2)){  
	      my $dist;
	      if ($chimera->strand2 == 1){
		  $dist = abs($chimera->pos2 - $candidate->{exon}->start);    
	      }else{
		  $dist = abs($chimera->pos2 - $candidate->{exon}->end);  
	      }
	      if ($exon_right_dist eq $NA || $dist < $exon_right_dist){
		  $exon_right_dist = $dist;
	      }
	  }
      }
  }
  return $exon_left_dist."---".$exon_right_dist;
}


=head2 getFusionGenes

Return fusion genes names separated with "---" characters

=cut

sub getFusionGenes {
  my $self = shift;
  my $chimera = shift;
  my $upper_gene = $self->getAnnotationAttribute($chimera,'before','gene','Name');
  my $lower_gene = $self->getAnnotationAttribute($chimera,'after','gene','Name');
  # Use ID attribute if no "Name" available for this gene
  $upper_gene = $self->getAnnotationAttribute($chimera,'before','gene','ID') if $upper_gene eq $CracTools::ChimCT::Const::NOT_AVAILABLE;
  $lower_gene = $self->getAnnotationAttribute($chimera,'after','gene','ID') if $lower_gene eq $CracTools::ChimCT::Const::NOT_AVAILABLE;
  return $upper_gene."---".$lower_gene;
}

=head2 getOutput

return Output Array

=cut

sub getOutput {
  my $self = shift;
  my $chimera = shift;
  croak "Missing argument chimera in getOutput" unless defined $chimera;
  
  my @output = ("mRNA_IDs=".$self->getAnnotationAttribute($chimera,'before','mRNA','ID')."---".$self->getAnnotationAttribute($chimera,'after','mRNA','ID'),
                "Annotation_type=".$self->getAnnotationType($chimera,'before')."---".$self->getAnnotationType($chimera,'after'),
                "Description_genes=".$self->getAnnotationAttribute($chimera,'before','mRNA','type')."---". $self->getAnnotationAttribute($chimera,'after','mRNA','type'),
                "Exon_IDs=".$self->getAnnotationAttribute($chimera,'before','exon','ID')."---".$self->getAnnotationAttribute($chimera,'after','exon','ID'),
              );
  
  # Computing Exon ranks
  my $exons_rank = $self->getExonsRank($chimera);

  # Computing Exons-end distances
  my $exons_distance = $self->getExonsDistance($chimera);

  push(@output,"Exon_end_distances=$exons_distance",
               "Exon_ranks=$exons_rank");

  return @output;
}


=head2 getScore

Return the annotation score for a given chimera

=cut

sub getScore {
  my $self = shift;
  my $chimera = shift;
  my $before_type = $self->getAnnotationAttribute($chimera,'before','mRNA','type');
  my $after_type = $self->getAnnotationAttribute($chimera,'after','mRNA','type');
  my $score = 100;
  my $comment;
  if (isPseudogene($before_type) || isPseudogene($after_type)) {
    $score = 0;
    $comment = 'Annotation=pseudogene';
  }
  if (isIGgene($before_type) && isIGgene($after_type)) {
    $comment = 'Annotation=IG-gene';
  }
  return ($score,$comment);
}

=head2 getAnnotation

Return the annotation hash for a given chimera : Hash(annot => (gene => (start => 1, end => 2, ...) ), priority => 2, type => 'CDS')

=cut

sub getAnnotation {
  my $self = shift;
  my $chimera = shift;
  croak "Missing argument chimera in getAnnotation" unless defined $chimera;

  # If this chimera has not been annotated yet we init its annotation hash
  if(defined $self->chimeraStruct->getChimera($chimera->getKey) && !defined $self->{annotations}{$chimera->getKey}) {
    $self->_initChimera($chimera); 
  }

  return $self->{annotations}{$chimera->getKey};
}

=head2 getAnnotationFeature

=cut

sub getAnnotationFeature {
    my $self = shift;
    my ($chimera,$side,$feature) = @_;
    my $annot_ref = $self->getAnnotation($chimera);
    return $annot_ref->{$side}{annot}->{$feature};
}

=head2 getAnnotationAttribute

=cut

sub getAnnotationAttribute {
    my $self = shift;
    my ($chimera,$side,$feature_name,$attr) = @_;
    my $feature = $self->getAnnotationFeature($chimera,$side,$feature_name);
    if(defined $feature && defined $feature->attribute($attr)) {
	return $feature->attribute($attr);
    } else {
	return $CracTools::ChimCT::Const::NOT_AVAILABLE;
    }
}

=head2 getAnnotationType

=cut

sub getAnnotationType {
  my $self = shift;
  my ($chimera,$side) = @_;
  croak "Unknown chimera key..." unless defined $self->getAnnotation($chimera);
  my $annot_ref = $self->getAnnotation($chimera);
  my $type = $annot_ref->{$side}{type};
  defined $type? return $type : return $CracTools::ChimCT::Const::NOT_AVAILABLE;
}

=head2 getAnnotationPriority

=cut

sub getAnnotationPriority {
  my $self = shift;
  my ($chimera,$side) = @_;
  croak "Unknown chimera key..." unless defined $self->getAnnotation($chimera);
  my $annot_ref = $self->getAnnotation($chimera);
  my $priority = $annot_ref->{$side}{priority};
  defined $priority? return $priority : return $CracTools::ChimCT::Const::NOT_AVAILABLE;
}

=head1 GETTERS & SETTERS

=head2 keepIG

=cut

sub keepIG {
  my $self = shift;
  my $keep_ig = shift;
  if(defined $keep_ig) {
    $self->{keep_ig} = $keep_ig;
  } elsif(!defined $self->{keep_ig}) {
    $self->{keep_ig} = 0;
  }
  return $self->{keep_ig};
}

=head2 isStranded

=cut

sub isStranded {
  my $self = shift;
  my $is_stranded = shift;
  if(defined $is_stranded) {
    $self->{is_stranded} = $is_stranded;
  }
  return defined $self->{is_stranded}? $self->{is_stranded}: 0;
}

=head2 estFile

=cut

sub estFile {
  my $self = shift;
  my $new_est_file = shift;
  if(defined $new_est_file) {
    $self->{est_file} = $new_est_file;
  }
  return $self->{est_file};
}

=head1 PRIVATE METHODS

=head2 _init

Init all chimeras contained in the chimeraStructure associated with this CracTools::ChimCT::Analyzer

=cut

sub _init {
  my $self = shift;

  # Check if the annotators are defined
  croak("Missing annotator in Annotation init") unless (defined $self->annotator);
  croak("Missing annotator_est in Annotation init") unless (defined $self->annotator_est || !defined $self->estFile);

  # Progress bar for annotation
  #my $count = 0;
  #my $nb_chimeras_initial = $self->chimeraStruct->nbChimeras;
  
  # my $flexibility = $CracTools::ChimCT::Const::FLEXIBILITY;
  # Get the annotation for each side of each chimera
  foreach my $chim (values %{$self->chimeraStruct->chimeras}) {
    $self->_initChimera($chim);
  }
}

=head2 _initChimera

Annotate the chimera passed in arguments

=cut

sub _initChimera {
  my ($self,$chim) = @_;

  my %annot;
  my ($pos1_start,$pos1_end,$pos2_start,$pos2_end) = ($chim->pos1,$chim->pos1,$chim->pos2,$chim->pos2);
  my $compare_sub = \&CracTools::Annotator::compareTwoCandidatesDefault;
  ## TODO: First, I decide to fix a flexixity, and I retracted myself because it does not seem mandatory to me. What do you think about that bro?  (Maybe, it makes sense for readthrough cases but i'm not pretty sure.
  # Compute pos1_start, pos1_end, pos1_start, pos1_end according to the strand

  # if ($chim->strand1 == 1){
  # 	$pos1_start = max(0,$pos1_start-$flexibility);
  # 	$pos2_start = max(0,$pos2_start-$flexibility);
  # }else{
  # 	$pos1_end += $flexibility;
  # 	$pos2_end += $flexibility;
  # }

  # Get left annotation according the gff_file, first sense otherwise antisense
  # A splicing junction may be inside an intron with a few nucleotides imprecision, so we apply a FLEXIBILITY parameter
  ## TODO: should have systematically check if a mRNA is defined? What should we do if only "five" or "three" is defined? (Maybe, we should check the consistency of our gff in a first step of the cracTools with systematically (gene,mRNA and exon), thus this test will make sense.
  my ($annot_before,$priority_before,$type_before) = $self->annotator->getBestAnnotationCandidate($chim->chr1,$pos1_start,$pos1_end,$chim->strand1,\&getCandidatePriority,$compare_sub);
  if(!defined $annot_before->{mRNA}) {
    ($annot_before,$priority_before,$type_before) = $self->annotator->getBestAnnotationCandidate($chim->chr1,$pos1_start,$pos1_end,$chim->strand1*-1,\&getCandidatePriority,$compare_sub); 
    $type_before .= '_antisense' if(defined $type_before);
  }
  $annot{before} = {annot => $annot_before, priority => $priority_before, type => $type_before};

  # Get right annotation according the gff_file, first sense otherwise antisense
  my ($annot_after,$priority_after,$type_after) = $self->annotator->getBestAnnotationCandidate($chim->chr2,$pos2_start,$pos2_end,$chim->strand2,\&getCandidatePriority,$compare_sub);
  if(!defined $annot_after->{mRNA}) {
    ($annot_after,$priority_after,$type_after) = $self->annotator->getBestAnnotationCandidate($chim->chr2,$pos2_start,$pos2_end,$chim->strand2*-1,\&getCandidatePriority,$compare_sub); 
    $type_after .= '_antisense' if(defined $type_after);
  }
  $annot{after} = {annot => $annot_after, priority => $priority_after, type => $type_after};

  # If there is still no annotation (even non-coding)
  # we try with EST
  if(defined $self->annotator_est && (!defined $annot_before->{mRNA} || !defined $annot_after->{mRNA})) {
    # Get left annotation according the est_file, first sense otherwise antisense
    if(!defined $annot_before->{mRNA}) {
      ($annot_before,$priority_before,$type_before) = $self->annotator_est->getBestAnnotationCandidate($chim->chr1,$pos1_start,$pos1_end,$chim->strand1,\&getCandidatePriority,$compare_sub); 
      if(!defined $annot_before->{mRNA}) {
        ($annot_before,$priority_before,$type_before) = $self->annotator_est->getBestAnnotationCandidate($chim->chr1,$pos1_start,$pos1_end,$chim->strand1*-1,\&getCandidatePriority,$compare_sub); 
        $type_before .= '_antisense' if(defined $type_before);
      }
      $annot{before} = {annot => $annot_before, priority => $priority_before, type => $type_before};
    }
    # Get right annotation according the est_file, first sense otherwise antisense
    if(!defined $annot_after->{mRNA}) {
      ($annot_after,$priority_after,$type_after) = $self->annotator_est->getBestAnnotationCandidate($chim->chr2,$pos2_start,$pos2_end,$chim->strand2,\&getCandidatePriority,$compare_sub); 

      if(!defined $annot_after->{mRNA}) {
        ($annot_after,$priority_after,$type_after) = $self->annotator_est->getBestAnnotationCandidate($chim->chr2,$pos2_start,$pos2_end,$chim->strand2*-1,\&getCandidatePriority,$compare_sub); 
        $type_after .= '_antisense' if(defined $type_after);
      }
      $annot{after} = {annot => $annot_after, priority => $priority_after, type => $type_after};
    }
  }

  # Check if the annotation is on the same gene for both parts of the chimera
  # when the chimera is a class 2. In this case we remove the chimera because
  # it coresponds to a splice.
  ## TODO?: Check neighboorhood for close gene that would lead to a splice (case of NA --- gene). If we do that, flexibility does definitely not make sense (in my opinion)
  my $remove_cause;
  my $remove_info;
  ## TODO: it could be interesting to do or not this sub-classification of the class2 according to a parameter (I have no suggestion for the parameter name)
  if($chim->getClass() == 2) {
    # 1/ Whatever the protocol, we remove splice events (class 2, sense, same genes)
    if($self->annotator->foundSameGene($chim->chr1,$chim->pos1,$chim->pos1,$chim->pos2,$chim->pos2,$chim->strand1)
      || (defined $self->annotator_est && $self->annotator_est->foundSameGene($chim->chr1,$chim->pos1,$chim->pos1,$chim->pos2,$chim->pos2,$chim->strand1))
    ){
      $remove_info = 'validated';
    }
    # 2/ With a non-stranded protocol, we remove splice events on opposite strand with (class 2, antisense, same genes) because they are considered as a splice event
    elsif(!$self->isStranded
      && 
      ($self->annotator->foundSameGene($chim->chr1,$chim->pos1,$chim->pos1,$chim->pos2,$chim->pos2,$chim->strand1*-1)
        || (defined $self->annotator_est && $self->annotator_est->foundSameGene($chim->chr1,$chim->pos1,$chim->pos1,$chim->pos2,$chim->pos2,$chim->strand1*-1))
      )
    ) {
      $remove_info = 'validated';
    }
    # 3/ Whatever the protocol, we remove new splice events (NA---gene, gene---NA or NA---NA) if |pos1-pos2| <= SPLICE_EVENT_THRESHOLD
    ## TODO: does mRNA the best feature to check a NA candidat? 
    elsif (!defined $annot_before->{mRNA} || !defined $annot_after->{mRNA}){
      if (abs($chim->pos1-$chim->pos2) <= $CracTools::ChimCT::Const::SPLICE_EVENT_THRESHOLD){
        $remove_info = 'new';
      }
    }
    # 4/ Whatever the protocol, we remove new events with an antisense annotation (AS---gene, gene---AS) if |pos1-pos2| <= SPLICE_EVENT_THRESHOLD
    elsif (
      (
        ($annot{before}->{type} =~ /_antisense/ && $annot{after}->{type} !~ /_antisense/) 
        || ($annot{before}->{type} !~ /_antisense/ && $annot{after}->{type} =~ /_antisense/)
      )
      && (abs($chim->pos1-$chim->pos2) <= $CracTools::ChimCT::Const::SPLICE_EVENT_THRESHOLD)
    ){
      $remove_info = 'new';
    }
    # 5/ With a stranded protocol, we remove antisense splice events (AS---AS) if |pos1-pos2| <= SPLICE_EVENT_THRESHOLD
    elsif($self->isStranded
      && ($annot{before}->{type} =~ /_antisense/ && $annot{after}->{type} =~ /_antisense/)
      && (abs($chim->pos1-$chim->pos2) <= $CracTools::ChimCT::Const::SPLICE_EVENT_THRESHOLD)
    ){
      $remove_info = 'antisense';
    } 

    # In all cases, when we found an information "new, antisense, validated", we talk about a splice event
    if(defined $remove_info) {
      $remove_cause = 'splice';
    }
  }

  # IG cause
  if(!defined $remove_cause && !$self->keepIG) {
    if ((defined $annot_before->{mRNA} && isIGgene($annot_before->{mRNA}->attribute('type'))) && (defined $annot_after->{mRNA} && isIGgene($annot_after->{mRNA}->attribute('type')))) {
      $remove_cause = 'ig_gene';
    }
  }

  # If a non-chimeric cause was found, we remove the chimera and we build a patch 
  if($remove_cause) {
    $self->chimeraStruct->removeChimera($chim,$remove_cause,$remove_info);
  }else{ 
    #$count++;
    $self->{annotations}{$chim->getKey} = \%annot;
  }
}
  
=head1 STATIC METHODS

=head2 isIGgene

Return true is the type givin in parameter corrspond to an IG gene.

=cut

sub isIGgene {
  my $gene_type = shift;
  if(defined $gene_type) {
    return $gene_type =~ /IG_/ || $gene_type =~ /TR_/;
  } else {
    return 0;
  }
}

=head2 isPseudogene

Return true is the type given in parameter corrspond to a pseudogene

=cut

sub isPseudogene {
  my $gene_type = shift;
  if(defined $gene_type) {
    return $gene_type =~ /pseudogene/i;
  } else {
    return 0;
  }
}

=head2 getCandidatePriority

  Arg [1] : String - pos_start
  Arg [2] : String - pos_end
  Arg [3] : hash - candidate

  Description : Method used to give a priority to a candidate in CracTools;:Annotator
                The best priority is 0. A priority of -1 means that this candidate
                should be avoided.
  ReturnType  : ($priority,$type)

=cut

sub getCandidatePriority {
  my ($pos_start,$pos_end,$candidate) = @_;
  my ($priority,$type) = (-1,$CracTools::ChimCT::Const::NOT_AVAILABLE);
  my ($mRNA,$exon) = ($candidate->{mRNA},$candidate->{exon});
  if(defined $mRNA) {
    if(!defined $mRNA->attribute('type') || $mRNA->attribute('type') =~ /protein_coding/i) {
      if(defined $exon) {
	$priority = 1;
        if(defined $candidate->{three}) {
          $type = '3PRIM_UTR';
        } elsif(defined $candidate->{five}) {
          $type = '5PRIM_UTR';
          # } elsif(defined $candidate->{cds}) {
          #   $type = 'CDS';
        } else {
          $type = 'EXON';
        }
      } else {
        $priority = 4;
        $type = 'INTRON';
      }
    } else {
      if(defined $exon) {
        $priority = 3;
        $type = 'NON_CODING';
      }
    }
  }
  return ($priority,$type);
}

sub priorityNeighborhood5prim {
  my ($pos_start,$pos_end,$candidate) = @_;
  my ($priority,$type) = (-1,'');
  my $gene = $candidate->{gene};
  if(defined $gene) {
    $priority = $pos_start - $gene->end;
    $type = '5PRIM_NEIGHBORHOOD';
  }
  return ($priority,$type);
}

sub priorityNeighborhood3prim {
  my ($pos_start,$pos_end,$candidate) = @_;
  my ($priority,$type) = (-1,'');
  my $gene = $candidate->{gene};
  if(defined $gene) {
    $priority = $gene->start - $pos_end;
    $type = '3PRIM_NEIGHBORHOOD';
  }
  return ($priority,$type);
}


=head1 GETTERS AND SETTERS

=cut

=head2 line 
    
  Description : Getter/Setter for the annotator object.
                If there is not gff_file, this method return
                the original annotator.
  ReturnType  : CracTools::Annotator object

=cut

sub annotator {
    my $self = shift;
    my $gff_file = shift;
    if (defined $gff_file){
      $self->{annotator} = CracTools::Annotator->new($gff_file);
    }
    return $self->{annotator}; 	
}

=head2 line 
    
  Description : Getter/Setter for the annotator object.
                If there is not gff_file, this method return
                the original annotator.
  ReturnType  : CracTools::Annotator object

=cut

sub annotator_est {
    my $self = shift;
    my $est_file = shift;
    if (defined $est_file){
      $self->{annotator_est} = CracTools::Annotator->new($est_file);
    }
    return $self->{annotator_est}; 	
}


1;
