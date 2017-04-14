# Cractools extract
rule cractools_extract:
    input:
        bam = BAM_DIR+"/{sample}.bam",
    output:
        chimeras = config['cractools']['output_dir']['chimeras']+"/{sample}_chimeras.tsv",
        splices = config['cractools']['output_dir']['splices'] + "/{sample}_splices.bed",
        mutations = config['cractools']['output_dir']['mutations'] + "/{sample}_mutations.vcf",
    log:
        stderr = config['cractools']['log_dir'] + "/{sample}_cractools.log",
        version = config['version_dir'] + "/cractools-version.txt",
    threads:
        #config['cractools']['nb_threads'],
        THREADS,
    message:
        "Executing cractools on {input}"
    shell:
        "cractools extract"
        " {input.bam} "
        + config['cractools']['options'] +
        " -r " + GENOME +
        " -p {threads}"
        " -s {output.splices}"
        " -c {output.chimeras}"
        " -m {output.mutations}"
        # Cractools logs
        " 2>{log.stderr}"
        # Cractools version
        " && cractools --version > {log.version}"

    


