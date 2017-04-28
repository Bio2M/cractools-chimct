# Cractools extract
rule cractools_extract:
    input:
        bam = BAM_DIR+"/{sample}.bam",
    output:
        chimeras = config['cractools']['output_dir']['chimeras'] + "/{sample}" + config['cractools']['output_suffix']['chimeras'],
        splices = config['cractools']['output_dir']['splices'] + "/{sample}" + config['cractools']['output_suffix']['splices'],
        mutations = config['cractools']['output_dir']['mutations'] + "/{sample}" + config['cractools']['output_suffix']['mutations'],
    log:
        stderr = config['cractools']['extract_log_dir'] + "/{sample}" + config['cractools']['extract_log_suffix'],
        version = config['cractools']['version_file'],
    benchmark:
        config['cractools']['extract_benchmark_dir'] + "/{sample}" + config['cractools']['extract_benchmark_suffix'],
    threads:
        #config['cractools']['nb_threads'],
        config['nb_threads'],
    params:
        binary = config['cractools']['binary'],
        options = config['cractools']['extract_options'],
        genome = config['genome'],
    message:
        "Executing cractools on {input}"
    shell:
        "{params.binary} extract"
        " {input.bam} "
        " {params.options}"
        " -r {params.genome}"
        " -p {threads}"
        " -s {output.splices}"
        " -c {output.chimeras}"
        " -m {output.mutations}"
        # Cractools logs
        " 2>{log.stderr}"
        # Cractools version
        " && cractools --version > {log.version}"

    


