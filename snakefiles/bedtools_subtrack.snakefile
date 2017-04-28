# Bedtools subtrack

"""
remove know genes

bedtools subtract \
­       -A \
    ­   -a unannotated_lncRNA_primarymerge.gtf \
    ­   -b /data/annotations/human/Homo_sapiens.GRCh38.86.gtf > unannotated_lincRNA_minusgenes.gtf
"""

rule bedtools_subtrack:
    input:
        a = config['stringtie']['merge_gtf_file'],
        b = config['gtf_ref'],
    output:
        config['bedtools']['subtrack_gtf_file'],
    log:
        stderr = config['bedtools']['subtrack_log_file'],
        version = config['bedtools']['version_file'],
    benchmark:
        config['bedtools']['subtrack_benchmark_file'],
    params:
        binary = config['bedtools']['binary'],
        options = config['bedtools']['subtract_options']
    message: 
        "Executing bedtools on file {input.a}"
    shell:
        "{params.binary} subtract"
        " {params.options}"
        " -a {input.a}"
        " -b {input.b}"
        " > {output}"
        # logs
        " 2> {log.stderr}"
        # bedtools version
        " && bedtools --version > {log.version}"

