# Bedtools subtrack
rule bedtools_subtrack:
    input:
        config['stringtie']['output_dir'] + "/unannotated_lncRNA_primarymerge.gtf",
    output:
        config['bedtools']['output_dir'] + "/" + config['bedtools']['output_subtrack'],
    log:
        stderr = config['bedtools']['log_dir'] + "/bedtools-subtrack.err",
        version = config['version_dir'] + "/bedtools-version.txt",
    message: 
        "Executing bedtools on file unannotated_lncRNA_primarymerge"
    shell:
        config['bedtools']['binary'] + 
        config['bedtools']['options'] +
        " -b" + config['gtf_ref'] +
        " -a {input}"
        " > {output}"
        # logs
        " 2> {log.stderr} |"
        # stringtie version
        " && echo 'stringtie version $(bedtools --version)' > {log.version}"


"""
bedtools subtract \
­       -A \
    ­-a unannotated_lncRNA_primarymerge.gtf \
    ­-b /data/annotations/human/Homo_sapiens.GRCh38.86.gtf > unannotated_lincRNA_minusgenes.gtf
"""
