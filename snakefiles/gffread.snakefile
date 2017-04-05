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
        version = config['version_dir'] + "/gffread-version.txt",
    message: 
        "Executing gffread"
    shell:
        config['gffread']['binary'] + 
        config['gffread']['options'] +
        " {input.gtf}"
        " -g {input.ref}"
        " -w {output}"
        # logs
        " 2> {log.stderr}"
        # gffread version
        # no version options, because part of cufflinks ?

