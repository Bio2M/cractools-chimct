# Stringtie assembling
rule stringtie:
    input:
        bam = config['raw_dir'] + "/{sample}.bam",
    output:
        gtf = config['stringtie']['output_dir'] + "/{sample}_stringtie.gtf",
    threads: 
        config['nb_threads']
    log:
        stderr = config['stringtie']['log_dir'] + "/{sample}_stringtie.log",
        version = config['stringtie']['version_file'],
    benchmark:
        config['stringtie']['benchmark_dir'] + "/{sample}_stringtie.benchmark",
    params:
        binary = config['stringtie']['binary'],
        options = config['stringtie']['options'],
        min_iso = config['stringtie']['min_iso'],
        min_len = config['stringtie']['min_length'],
        gtf_ref = config['gtf_ref'],
    message: 
        'Executing stringtie on the following files {input}'
    shell:
        '{params.binary}'
        ' {params.options}'
        ' {input.bam}'
        ' -f {params.min_iso}'
        ' -m {params.min_len}'
        ' -p {threads}'
        ' -G {params.gtf_ref}'
        ' -o {output.gtf}'
        # logs
        ' 2> {log.stderr} '
        # stringtie version
        ' && echo "stringtie version $(stringtie --version)" > {log.version}'

# Stringtie merge
rule stringtie_merge:
    input:
        gtf = expand(config['stringtie']['output_dir'] + "/{sample}_stringtie.gtf", sample=SAMPLES)
    output:
        gtf = config['stringtie']['merge_gtf_file'],
    log:
        stderr = config['stringtie']['merge_log_file'],
        version = config['stringtie']['version_file'],
    benchmark:
        config['stringtie']['merge_benchmark_file'],
    params:
        binary = config['stringtie']['binary'],
        options = config['stringtie']['options'],
        min_iso = config['stringtie']['min_iso'],
        min_len = config['stringtie']['min_length'],
        gtf_ref = config['gtf_ref'],
    message: 
        "Executing stringtie on the following files {input}"
    shell:
        "{params.binary}"
        " {params.options}"
        " --merge {input.gtf}"
        " -f {params.min_iso}"
        " -m {params.min_len}"
        " -G {params.gtf_ref}"
        " -o {output.gtf}" 
        # logs
        " 2> {log.stderr}"
        # stringtie version
        " && echo 'stringtie version $(stringtie --version)' > {log.version}"
