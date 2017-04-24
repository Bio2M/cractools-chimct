# freebayes variant caller
rule freebayes:
    input:
        bam = config['hisat2']['output_dir'] + '/{sample}.bam'
    output:
        vcf = config['freebayes']['output_dir'] + '/{sample}.vcf'
    threads: 
        config['nb_threads']
    log:
        stderr = config['freebayes']['log_dir'] + '/{sample}_freebayes.log',
        version = config['version_dir'] + '/freebayes-version.txt',
    benchmark:
        "output/benchmarks/freebayes/{sample}_freebayes.benchmark"
    message: 
        "Executing freebayes on the following files {input.bam}"
    shell:
        "freebayes"
        + config['freebayes']['options'] +
        " -t " + config['bed_file'] +
        " -f " + config['genome'] +
        " {input.bam}"
        " > {output.vcf}"
        # logs
        " 2> {log.stderr}"
        # version
        " && freebayes --version > {log.version}"
