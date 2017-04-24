# Crac alignment
rule crac_pe:
    input:
        r1 = FASTQ_DIR + "/{sample}_1.fastq.gz",
        r2 = FASTQ_DIR + "/{sample}_2.fastq.gz",
    output:
        bam = BAM_DIR + "/{sample}.bam",
        bai = BAM_DIR + "/{sample}.bam.bai",
        summary = config['crac']['summary_dir'] + "/{sample}.summary",
    threads: 
        config['nb_threads'],
    params:
        tmp = "/tmp/{sample}",
    log:
        stderr = config['crac']['log_dir'] + "/{sample}_crac.log",
        version = config['version_dir'] + "/crac-version.txt",
    benchmark:
        "output/benchmarks/crac/{sample}_crac.benchmark"
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
        # Samtools sam to bam
        " | samtools view -bS -"
        # Samtools sort
        " | samtools sort "
        " -T {params.tmp} "
        " -@ " + config['samtools']['threads'] +
        " -m " + config['samtools']['mem'] +
        " -" 
        " -o {output.bam}"
        # samtools index
        " && samtools index {output.bam} {output.bai} "
        # crac version
        " && crac -v | head -1 > {log.version}"
    
