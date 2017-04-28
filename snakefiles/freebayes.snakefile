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
        version = config['freebayes']['version_file'],
    benchmark:
        config['freebayes']['benchmark_dir'] +  '/{sample}_freebayes.benchmark',
    params:
        binary = config['freebayes']['binary'],
        options = config['freebayes']['options'],
        bed_ref = config['freebayes']['bed_ref_file'],
        genome = config['genome'],
    message: 
        "Executing freebayes on the following files {input.bam}"
    shell:
        "{params.binary}"
        " {params.options}"
        " -t {params.bed_ref}"
        " -f {params.genome}"
        " {input.bam}"
        " > {output.vcf}"
        # logs
        " 2> {log.stderr}"
        # version
        " && freebayes --version > {log.version}"
