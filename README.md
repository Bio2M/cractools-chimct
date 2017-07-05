# What is it?

The CracTools-ChimCT software (or simply ChimCT) extracts, formats, classifies and annotes chimeras identified by CRAC.


# Installation

## Requirements

* Perl5 distribution
* [cpanm](http://search.cpan.org/~miyagawa/App-cpanminus-1.7043/lib/App/cpanminus.pm) (CPAN minus)
* [CracTools-core](https://metacpan.org/release/CracTools) perl package. It will be automatically installed by cpanm along with all other CPAN dependancies.

## Install from tarball

This is the simpliest way to install chimCT.

1. Go to the [release page](https://github.com/Bio2M/cractools-chimct/releases) of the github projet.
1. Download the latest tarball (tar.gz) release : ``wget https://github.com/Bio2M/cractools-chimct/releases/download/v0.14/CracTools-chimCT-0.14.tar.gz``
2. Install the package with cpanm : ``cpanm [-l local_dir] CracTools-chimCT-0.14.tar.gz``

If you do not have admin rights, you can use the option ``-l`` to specify cpanm a local directory to install simCT.

## Install from sources

To install chimCT from the sources, you will need [Dist::Zilla](http://dzil.org/) software, which is a Perl package manager.

1. Clone chimCT repository : ``git clone https://github.com/Bio2M/cractools-chimct.git``
2. Build and install : ``dzil install --install-command 'cpanm [-l LOCAL_INSTALLATION_DIRECTORY].'`` (ommit the `-l` option if you want to install chimCT to the system).



# Documentation


### NAME  

chimCT - Extract, format, classify and annote chimeras identified by CRAC.

### VERSION  

version 0.14

### SYNOPSIS  

chimCT -g file.gff -s file.sam -n sample_name > output.csv

### DESCRIPTION
    
The chimera pipeline is a tool for both annotating and formating chimera 
produced by CRAC software as well as computing a score of beauty for each
chimera and sort them using this criteria.

### OPTIONS

**General**
      --help                  Print this help
      --man                   Open man page
      --version               Print chimCT version number

**Mandatory arguments**
  
      -g,--gff=file.gff       Specify GFF3 file to perform annotation with. (or declare ANNOTATION_GFF in conf file)
      -s,--sam=file.sam       Specify SAM file used to extract chimeras.
      -n,--sample=name        Specify Sample name

**Optionnal arguments**
  
      --stranded              Make use of stranded information to output chimeras in the right orientation. Stranded PE reads are interpreted as reverse-forward reads.
      --keep-ig               Keep chimeras that correspond to fusions IG genes.
      --est                   Specify GFF3 file with EST annotations. These annotations will be used if we don't find any matching annotation in the --gff annotation file.
      --verify-splice         Verify splices called by CRAC using the GFF3 annotation and convert them into chimera if needed.
      --tmp-dir               Specify a location were the temporary files (like GSNAP mapping files) should be kept

**Output : files & formats**
      --detailed-sam          Print the p_support and the p_loc profiles for the selected read.
      --primers               Prints the sequence of concatenated paired-end reads with junction position highlighted for easier primers design.
      --summary=summary.txt   Specify file to print summary stats about the chimeraPipeline in a separate output file.
      --spanning-reads=file   Specify basename to output chimeric spanning reads (PE and single) is fasta format (the fasta extension is automatically added to the end of the filename.

**Analyzer-specific options**

      --gsnap-softclip-threshold=Integer  Minimum number of nucleotides that have to be mapped on both sides of the breakpoint gave by CRAC in order to consider a chimera as "MAPPED" (DEFAULT 15)

**Configuration**

      --gsnap-exe                 Specify Gsnap exe file
      --gsnap-genome-directory    Specify GSNAP genome directory
      --gsnap-genome-name         GSNAP genome database
      --gsnap-nb-threads=Integer  Number of threads used to run GSNAP (DEFAULT 4)
      --keep-gsnap-output         Keep GSNAP output mapping files an print a link on STDERR to reach them

#### INPUT (SAM) FILE

Chimeras identified by CRAC are loaded from the SAM file in argument. If
option ``--stranded`` is specified, PE reads will be treat as
"reverse-forward". Chimeras identified on the first read of a pair will be
reversed.

#### ANALYZERS

The chimera pipeline is construct around a concept of "Analyzer". Each
Anlayzer process the chimeras extracted from the SAM file in order to give
additional information and to compute a score that will contribute to the
chimera "Rank". Using those different pieces of information, we can
classify and order the chimeras using this "confidence" rank.

Analyzers are described below.

####  ANNOTATION

Annotation is the first step of the chimeraPipeline. Annotation is
performed using a GFF3 file provided as command line argument. Such a file
can be generated using "buildGFF3FromEnsembl.pl" script of the
CracTools-Core distribution.

The annotation analyzer aims to find the annotations (gene -> mRNA ->
exon) of both chimera parts. We first look for Protein coding transcript,
then non-coding transcript then EST transcript when ``--est`` GFF3 file
provided.

Annotation analyzer will remove class 2 chimeras that are actually splices.

####  CRAC SCORE

CRAC SCORE analyzer takes the algorithmic score given by CRAC for each
read that contains the chimera and integrates it in the chim_value.

**Warning**: this analyzer is only available for CRAC versions greater than
(or equal to) 1.9.1.

####  CLASSIFICATION

Analyze and score each chimeras. Close to 1, the output will be very
sensitive. Close to 0, the output will be very specific but less
sensisitive.

####  GSNAP MAPPING

GSNAP Mapping analyzer will try to map chimeric reads with another mapping
tool (GSNAP) to discard false positives.

In order to use GSNAP Mapping you need to defined some values in the
configuration file, otherwise GSNAP will not be used.

      GSNAP_EXE '/usr/bin/gsnap'
      GSNAP_GENOME_DIRECTORY '/data/indexes/gsnap'
      GSNAP_GENOME_NAME 'GRCh37'
      GSNAP_NB_THREADS 2

####  PAIRED END COVERAGE

If you are using PE reads we expect to find a consistent Paired-end
coverage for the chimeras. This piece of information will be used to
compute the "chim_value".

#### STRINGENT FLAGS

#### FUSION DISTANCE (for class 3 chimeras only)
    
Fusion distance in class 3 chimeras is an important factor to discrimate
artefacts.

####  PRIMERS DESIGN
#### CLASSIFICATION & CHIM_VALUES
    We define a chim_value for each chimera combining scores given by the
    analyzers. Depending on the chimera class, weights are applied for each
    analyzers (see table below).

|                 | Class 1 | Class 2 | Class 3 | Class 4 |
| --------------- | :-----: | :-----: | :-----: | :-----: |
| Annotation      |   xxxx  |   xxxx  |   xxxx  |   xxxx  |
| GSNAP Mapping   |     xx  |         |   x     |      x  |
| PE coverage     |    xxx  |     xx  |         |     xx  |
| Stringent flags |     xx  |      x  |      x  |      x  |
| Fusion distance |         |         |     xx  |         |


#### OUTPUT FORMAT

This is a description of the format used to output analysed chimeras. Each
line correspond to a uniq chimera identified and annotated by the
chimera Pipeline. Chimeras are ordered by Score, then class, then number of
spanning reads.

This format is composed by 13 mandatory fields TAB-separeted. This is an
home-made format, because there is no standard format that is able to
store chimeras coordinates and related informations.

1. **Id** - A Uniq Id for each chimera. This id is composed by "sample_name:chimera id".  
2. **Name** - Fusion genes names separated by three dashes ('---')  
3. **Chr1** - Chromosome of the 5' part of the chimera  
4. **Pos1** - Genomic positions of the 5' part of the chimera  
5. **Strand1** - Genomic strand of the 5' part of the chimera. If sample is not
    ``--stranded`` No assumption can be made about the strand  
6. **Chr2** - Chromosome of the 3' part of the chimera. Same as Chr2, unless it
    is a class 1 chimera (translocation).  
7. **Pos2** - Genomic positions of the 3' part of the chimera  
8. **Strand2** - Genomic strand of the 3' part of the chimera. If sample is not
    ``--stranded`` No assumption can be made about the strand  
9. **Chim_value** - Confidence value about chimera's positivity (empirically
    constructed)  
10. **Spanning_junction_normalized** - Spaning junction reads coverage
    (normalized per billion of reads). A spanning junction read is the read
    that contains the chimeric junction.  
11. **Spanning_PE_normalized** - Coverage of paired-end reads (normalized per
    billion of reads) that contains the chimeric junction in the non-sequenced
    part.  
12. **Class** - Chimeric class from 1 to 4. (add more details)  
13. **Comments** - Free comments about the chimera. This field try to explain
    textualy the "Rank" value given for the chimera.  
14. **Others** - This last fields are dedicated to additional informations.
    They are based one the {Key='value'} paradigm.  

### AUTHORS

*   Jérôme AUDOUX <jaudoux@cpan.org>
*   Nicolas PHILIPPE <nicolas.philippe@inserm.fr>
*   Sacha BEAUMEUNIER <sacha.beaumeunier@live.fr>

### COPYRIGHT AND LICENSE

This software is Copyright (c) 2016 by IRB/INSERM (Institut de Recherche
en Biothérapie / Institut National de la Santé et de la Recherche
Médicale).

This is free software, licensed under:

GPL_3 SOFTWARE LICENSE AGREEMENT
