package CracTools::ChimCT::Const;
# ABSTRACT: Constants for package CracTools-chimCT

use strict;
use warnings;

# Symbol to use when priting a data that is not available
our $NOT_AVAILABLE = 'N/A';

# Minimum (nomalized) support to report the chimera as "low_support"
our $MIN_SUPPORT = 20;

# Maximum (nomalized) support to report the chimera as "high_support"
our $MAX_SUPPORT = 5000;

# Normalization value
our $NORMALIZATION_CONST = 1000000000;

# Number of thread to use when running GSNAP
# TODO (this should disapear when chimCT will be multi-threaded)
our $GSNAP_NB_THREADS   = 4;

# Minimum number of nuclotides that have to be mapped in each side of the junction
# in order to determine if the chimera is really mapped
our $GSNAP_SOFTCLIP_THRESHOLD    = 15;

# Maximum overlapping distance between two couples of paired-end reads to be
# considered as similar
our $PAIRED_MAX_OVERLAPPING_DISTANCE = 10000;

# Fusion distance used for computing score of Class 3 chimeras.
our $FUSION_DISTANCE_THRESHOLD = 1000;

# Maximum tolerated (genomic) distance used to considered a splice events
our $SPLICE_EVENT_THRESHOLD = 300000;

# A splicing junction may be inside an intron with a few nucleotides imprecision. This parameter is the maximum of flexibility tolerated
# our $FLEXIBILITY = 10;

# Maximum length of th DNA fragment used to sequence paired-end reads
our $MAX_TEMPLATE_LENGTH =  500;

# Maximum number of reads to keep for one chimeras
# The bigger this constant is, the more memory chimCT will consume
our $MAX_READS_BY_CHIMERA = 50;

1;
