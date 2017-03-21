# Cractools extract
rule cractools_extract:
    input:
        bam = BAM_DIR+"/{sample}.bam",
    output:
        chimeras = config['cractools']['output_dir']['chimeras']+"/{sample}_chimeras.tsv",
        splices = config['cractools']['output_dir']['splices'] + "/{sample}_splices.bed",
        mutations = config['cractools']['output_dir']['mutations'] + "/{sample}.vcf",
    log:
        stderr = config['cractools']['log_dir'] + "/{.sample}_cractools.log",
        version = config['report_versions'],
    threads:
        config['nb_threads']
    message:
        "Executing cractools on {input}"
    run:
        if config['genome']:
            params += " -r " + config['genome']
        shell(
            "cractools extract "
            + config['cractools']['options'] +
            " -p {threads}"
            " -s {output.splices}"
            " -c {output.chimeras}"
            " -m {output.mutations}"
            # Cractools logs
            " 2>{log.stderr}"
            # Cractools version
            " && cractools --version >> {log.version}"
        )
    


