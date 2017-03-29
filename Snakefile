# Global snakefile
"""
Date: Fabruary 2017
"""

import subprocess

# Config file
configfile: "samples.yml"
configfile: "config.yml"


# variables from configfile
# fastx parameters
FASTQ_DIR = config['fastq_dir']
SAMPLES = expand(config['samples'])
# bam parameters
BAM_DIR = config['hisat2']['output_dir']
# gff - gtf parameters
ANNOTATION_FILE = config['annotation_file']
# tsv params
TSV_FILES = expand(config['chimct']['output_dir']+"/{sample}"+config['chimct']['output_suffix'], sample=SAMPLES)
# others
EMAIL = config['email']
SNAKEMAKE_VERSION = subprocess.check_output("snakemake --version", shell=True)


# includes

#include: "snakefiles/crac_pe.snakefile"
include: "snakefiles/cractools_extract.snakefile"
include: "snakefiles/chimct.snakefile"
include: "snakefiles/extract_common_chimeras.snakefile"
include: "snakefiles/hisat2_pe.snakefile"
include: "snakefiles/freebayes.snakefile"


rule all:
    input:
        #expand(config['hisat2']['output_dir']+"/{sample}.bam", sample=config['samples']['paired']),
        expand(config['freebayes']['output_dir']+"/{sample}.vcf", sample=config['samples']),
        expand(config['chimct']['output_dir']+"/{sample}_chimct.tsv", sample=SAMPLES),
        expand(config['cractools']['output_dir']['mutations']+"/{sample}.vcf", sample=SAMPLES),
        config['extract_common_chimeras']['output_dir'] + "/common_chimeras.txt",
    log:
        version = config['version_dir'] + "/snakemake.txt",
    shell:
        "echo snakemake versions $(snakemake --version) > {log.version}"
    #version:
    #    SNAKEMAKE_VERSION

onsuccess:
    print("Workflow finished, no error")
    #shell( 'mail -s "$HOSTNAME: Workflow finished, no error" ' + EMAIL + ' < {log}')

onerror:
    print("An error occured, pas glop!")
    #shell('mail -s "$HOSTNAME: an error occurred on Workflow" ' + EMAIL + ' < {log}')

