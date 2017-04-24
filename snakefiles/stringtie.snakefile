# Stringtie assembling
rule stringtie:
    input:
        bam = BAM_DIR + "/{sample}.bam",
    output:
        gtf = config['stringtie']['output_dir'] + "/{sample}.gtf",
    threads: 
        config['nb_threads']
    log:
        stderr = config['stringtie']['log_dir'] + "/{sample}_stringtie.log",
        version = config['version_dir'] + "/stringtie-version.txt",
    benchmark:
        "output/benchmarks/stringtie/{sample}_stringtie.benchmark"
    message: 
        "Executing stringtie on the following files {input}"
    shell:
        config['stringtie']['binary'] +
        " {input.bam}"
        " -f " + config['stringtie']['min_iso'] +
        " -m " + config['stringtie']['min_length'] +
        " -p {threads}"
        " -G " + config['gtf_ref'] +
        " -o {output.gtf}" 
        # logs
        " 2> {log.stderr} "
        # stringtie version
        " && echo 'stringtie version $(stringtie --version)' > {log.version}"
