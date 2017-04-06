# Crac alignment
rule crac_pe:
    input:
        r1 = FASTQ_DIR + "/{sample}_1.fastq.gz",
        r2 = FASTQ_DIR + "/{sample}_2.fastq.gz",
    output:
        bam = BAM_DIR + "/{sample}.bam",
        bai = BAM_DIR + "/{sample}.bam.bai",
        summary = SUMMARY_CRAC_DIR + "/{sample}.summary",
    threads: 
        THREADS,
    params:
        tmp = "/tmp/{sample}",
    log:
        stderr = config['crac']['log_dir'] + "/{sample}_crac.log",
        version = config['version_dir'] + "/crac-version.txt",
    benchmark:
        "benchmarks/{sample}.bwa.benchmark.txt"
    message: 
        "Executing crac on the following files {input}"
    shell:
        config['crac']['binary'] + " "
        + config['crac']['options'] +
        " -i " + config['crac']['index'] +
        " -k " + config['crac']['kmer_len'] +
        " --summary {output.summary}"
        " --nb-threads {threads}"
        " -r {input.r1} {input.r2}"
        " -o - "
        # crac logs
        " 2> {log.stderr}"
        # Samtools sort
        " | samtools sort "
        " -T {params.tmp} "
        " -@" + config['samtools']['threads'] +
        " -m " + config['samtools']['mem'] +
        " - -o {output.bam}"
        # samtools index
        " && samtools index {output.bam} {output.bai} "
        # crac version
        " && crac -v | head -1 > {log.version}"
    
