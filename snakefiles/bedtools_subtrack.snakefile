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
        a = config['stringtie']['output_dir'] + "/unannotated_lncRNA_primarymerge.gtf",
        b = config['gtf_ref'],
    output:
        config['bedtools']['output_subtrack'],
    log:
        stderr = config['bedtools']['log_dir'] + "/bedtools-subtrack.err",
        version = config['version_dir'] + "/bedtools-version.txt",
    message: 
        "Executing bedtools on file unannotated_lncRNA_primarymerge"
    shell:
        config['bedtools']['binary'] + 
        config['bedtools']['subtract_options'] +
        " -a {input.a}"
        " -b {input.b}"
        " > {output}"
        # logs
        " 2> {log.stderr} |"
        # bedtools version
        " && bedtools --version > {log.version}"

