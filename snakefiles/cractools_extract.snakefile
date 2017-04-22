# Cractools extract
rule cractools_extract:
    input:
        bam = BAM_DIR+"/{sample}.bam",
    output:
        chimeras = config['cractools']['output_dir']['chimeras'] + "/{sample}" + config['cractools']['output_suffix']['chimeras'],
        splices = config['cractools']['output_dir']['splices'] + "/{sample}" + config['cractools']['output_suffix']['chimeras'],
        mutations = config['cractools']['output_dir']['mutations'] + "/{sample}" + config['cractools']['output_suffix']['chimeras'],
    log:
        stderr = config['cractools']['log_dir'] + "/{sample}-cractools.log",
        version = config['version_dir'] + "/cractools-version.txt",
    threads:
        #config['cractools']['nb_threads'],
        config['nb_threads'],
    message:
        "Executing cractools on {input}"
    shell:
        "cractools extract"
        " {input.bam} "
        + config['cractools']['options'] +
        " -r " + config['genome'] +
        " -p {threads}"
        " -s {output.splices}"
        " -c {output.chimeras}"
        " -m {output.mutations}"
        # Cractools logs
        " 2>{log.stderr}"
        # Cractools version
        " && cractools --version > {log.version}"

    


