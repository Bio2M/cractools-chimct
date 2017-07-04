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
#                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               #
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

use strict;
use warnings;

package CracTools::ChimCT;
# ABSTRACT: Extract, format, classify and annote chimeras.

use CracTools 1.091;
#use Carp;
#use CracTools::ChimCT::Structure;
#use CracTools::ChimCT::OverlapStructure;
#
#=head1 METHODS
#
#=head2 new
#
#ChimCT constructor
#
#=cut
#
#sub new {
#  my $class = shift;
#
#  my %args = @_;
#
#  croak "Missing annot_analyzer argument" unless defined $args{annot_analyzer};
#  croak "Missing is_stranded argument" unless defined $args{is_stranded};
#  
#  my $self = bless {
#    is_stranded => $args{is_stranded},
#    chimera_struct => CracTools::ChimCT::Structure->new(),
#    false_chimera_struct => CracTools::ChimCT::Structure->new(),
#    paired_chimera_overlap_struct => CracTools::ChimCT::OverlapStructure->new(is_stranded => $args{is_stranded}),
#    annot_analyzer => $args{annot_analyzer},
#    nb_reads_import => 0,
#  }, $class;
#
#  return $self;
#}
#
#=head2 importFromSAM
#
#Load chimeras from the SAM file in arguments
#
#=cut
#
#sub importFromSAM {
#  my $self = shift;
#  my $sam_file = shift;
#  my $sam_reader = CracTools::SAMReader->new($sam_file);
#  my $nb_reads_import = 0;
#  # Loop over each line of the SAM file
#  my $sam_it = $sam_reader->iterator();
#  while (my $sam_line = $sam_it->()) {
#    $nb_reads_import++;
#
#    # Loop over each chimera of the current SAM line
#    foreach my $chim_event (@{$sam_line->events('chimera')}) {
#
#      # Get chimeras positions
#      my ($chr1,$pos1,$strand1) = @{$chim_event->{loc1}}{'chr','pos','strand'};
#      my ($chr2,$pos2,$strand2) = @{$chim_event->{loc2}}{'chr','pos','strand'};
#      my $chimera = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
#      my $seq = $sam_line->seq;
#      my $read = CracTools::ChimCT::Read->new(id => $sam_line->qname,
#                                                       seq => $sam_line->getOriginalSeq, 
#                                                       pos_junction => $chim_event->{pos}, 
#                                                       sam_line => $sam_line, 
#                                                       chimera_event => $chim_event);
#
#      $self->_processChimera($chimera,$read,$sam_line);
#    }
#    
#    # If "verify-splice" option is activated we also add Junction as chimeras
#    # that will be verified and discarded using the annotation process
#    if($verify_splice) {
#      foreach my $splice (@{$sam_line->events('Junction')}) {
#
#        # Get chimeras positions
#        my ($chr,$pos,$strand) = @{$splice->{loc}}{'chr','pos','strand'};
#
#        my ($pos1,$pos2);
#        if($strand == -1) {
#          $pos1 = $pos+$splice->{gap}+1;
#          $pos2 = $pos;
#        } else {
#          $pos1 = $pos;
#          $pos2 = $pos+$splice->{gap}+1;
#        }
#
#        my $chimera = CracTools::ChimCT::Chimera->new($chr,$pos1,$strand,$chr,$pos2,$strand);
#        my $seq = $sam_line->seq;
#        my $read = CracTools::ChimCT::Read->new(id => $sam_line->qname,
#                                                         seq => $sam_line->getOriginalSeq, 
#                                                         pos_junction => $splice->{pos}, 
#                                                         sam_line => $sam_line,
#                                                       );
#
#        $self->_processChimera($chimera,$read,$sam_line);
#      }
#    }
#
#    # If there is a pairedChimera we add it to the overlap structure
#    if (defined $sam_line->pairedChimera() &&
#        defined $sam_line->isClassified('duplicated') && $sam_line->isClassified('duplicated') == 0 &&
#        defined $sam_line->isClassified('multiple') && $sam_line->isClassified('multiple') == 0 &&
#        defined $sam_line->isPairedClassified('duplicated') && $sam_line->isPairedClassified('duplicated') == 0 &&
#        defined $sam_line->isPairedClassified('multiple') && $sam_line->isPairedClassified('multiple') == 0) {
#      my ($chr1,$pos1,$strand1,$chr2,$pos2,$strand2) = $sam_line->pairedChimera();
#      # Reverse strand2 (Paired-end protocol)
#      $strand2 = -$strand2;
#      my $PE_chim = CracTools::ChimCT::Chimera->new($chr1,$pos1,$strand1,$chr2,$pos2,$strand2);
#      # Because PE chimeras are only given for the first read of the pair
#      # and since stranded protocol (by illumina) is designed in order to have
#      # the good strand for the mate, we need to reverse the chimera
#      $PE_chim->reverseChimera(); # Test reversing chimera for stranded protocol
#      $paired_chimera_overlap_struct->addChimera($PE_chim);
#      
#    # If verify splice is "on" and this pair correspond to a splice
#    # We always use the second read of the pair, which is flagged 'LAST_SEGMENT'
#    # because in case of stranded protocol, this one give the rght orientation
#    } elsif($verify_splice && $sam_line->rnext eq '=' && abs($sam_line->tlen) > $CracTools::ChimCT::Const::MAX_TEMPLATE_LENGTH && 
#        $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{LAST_SEGMENT})) {
#      my $strand = 1;
#      $strand = -1 if $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{REVERSE_COMPLEMENTED});
#      my $PE_chim = CracTools::ChimCT::Chimera->new($sam_line->chr,$sam_line->pos,$strand,$sam_line->chr,$sam_line->pnext,$strand);
#      $paired_chimera_overlap_struct->addChimera($PE_chim);
#    }
#  }
#  $chimera_struct->nbReadsImported($nb_reads_import); # For normalization
#}
#
#=head1 GETTERS & SETTERS
#
#=head2 isStranded
#
#Return true is chimCT use "stranded" data.
#
#=cut
#
#sub isStranded {
#  my ($self,$new_value) = @_;
#  $self->{is_stranded} = $new_value unless !defined $new_value;
#  return $self->{is_stranded};
#}
#
#=head2 chimeraStruct
#
#Return the CracTools::ChimCT::Structure used to store chimeras.
#
#Warning: getter ONLY!
#
#=cut
#
#sub chimeraStruct {
#  my $self = shift;
#  return $self->{chimera_struct};
#}
#
#=head2 falseChimeraStruct
#
#Return the CracTools::ChimCT::Structure used to store 'false positive' chimeras
#
#Warning: getter ONLY!
#
#=cut
#
#sub falseChimeraStruct {
#  my $self = shift;
#  return $self->{false_chimera_struct};
#}
#
#=head2 annotAnalyzer
#
#Return the CracTools::ChimCT::Analyzer::Annotation used to perform annotation of chimeras
#
#Warning: getter ONLY!
#
#=cut
#
#sub annotAnalyzer {
#  my $self = shift;
#  return $self->{annot_analyzer};
#}
#
#
#=head1 PRIVATE METHODS
#
#=head2 _processChimera
#
#Add a chimera with its read and the corresponding CracTools::SAMReader::SAMline object
#
#=cut
#
## Correct the orientation of the chimera based on the "is_stranded option"
## Force Annot_analyzer to annotate the chimera
#sub _processChimera {
#
#  my ($self,$chimera,$read,$sam_line) = @_;
#  # We introduce a boolean for reversing the read in order to do it only if the we are considering this chimera
#  # (this will not be the case of a splice for example) and reversig a read is time consuming because we
#  # reverse complemente the whole sequence
#  my $reverse_read = 0;
#
#  if($self->isStranded && $sam_line->isFlagged($CracTools::SAMReader::SAMline::flags{FIRST_SEGMENT})) {
#    $reverse_read = 1;
#    $chimera->reverseChimera();
#  } elsif (!$is_stranded && defined $self->chimeraStruct->getChimera($chimera->getReverseKey)) {
#    $reverse_read = 1;
#    $chimera->reverseChimera();
#  }
#
#  my $existing_chimera = $self->chimeraStruct->getChimera($chimera->getKey());
#  if(defined $existing_chimera) {
#    $read->reverseRead() if $reverse_read;
#    $existing_chimera->addRead($read);
#  } elsif(!defined $self->falseChimeraStruct->getChimera($chimera->getKey)) {
#    $read->reverseRead() if $reverse_read;
#    $chimera->addRead($read);
#    $self->annotAnalyzer->addChimera($chimera);
#    # We force the annot analyzer to annotate the newly added chimera
#    #$annot_analyzer->_initChimera($chimera);
#    # If no annotation is available that means that this chimeras has been
#    # removed from the chimer_struct due to its annotation.
#    # We add the coordinates to the false positive chimera_struct in order to
#    # keep the information about this false positive if we found a other read that
#    # has the same chimeric coordinates
#    if(!defined $self->annotAnalyzer->getAnnotation($chimera)) {
#      $self->falseChimeraStruct->addChimera($chimera);
#    }
#  }
#}

1; 
