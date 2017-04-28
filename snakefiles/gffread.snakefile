# gffread

"""
On en fait ensuite un Fasta qui sera indexé par kallisto plus tard:

gffread unannotated_MSClinc_v2.gtf \
    - ­g /data/genomes/GRCh38/GRCh38.fa \
    ­-w unannotated_MSClinc_v2.fa
"""

rule gffread:
    input:
        gtf = config['gtf_exon_remove']['output_file'],
        ref = config['genome'],
    output:
        config['gffread']['output_file'],
    log:
        stderr = config['gffread']['log_file'],
        version = config['gffread']['version_file'],
    benchmark:
        config['gffread']['benchmark_file'],
    params:
        binary = config['gffread']['binary'],
        options = config['gffread']['options'],
    message: 
        "Executing gffread on {input.gtf}"
    shell:
        "{params.binary}"
        " {params.options}"
        " {input.gtf}"
        " -g {input.ref}"
        " -w {output}"
        # logs
        " 2> {log.stderr}"
        # gffread version
        # no version options, because part of cufflinks ?

