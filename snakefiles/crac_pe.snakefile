# Crac alignment
rule crac_pe:
    input:
        r1 = config['reads_dir']+"/{sample}_1.fastq.gz",
        r2 = config['reads_dir']+"/{sample}_2.fastq.gz",
    output:
        bam = config['crac']['output_dir']+"/{sample}.bam",
        bai = config['crac']['output_dir']+"/{sample}.bam.bai",
        summary = config['crac']['summary_dir'] + "/{sample}_summary.txt",
    threads: 
        config['nb_threads']
    log:
        stderr = config['crac']['log_dir'] + "/{sample}_crac.log",
        version = config['report_versions'],
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
        # Samtools sort & index
        " | samtools sort -@4 - -o {output.bam}"
        " && samtools index {output.bam} {output.bai} "
        # crac version
        " && crac -v | head -1 >> {log.version}"
    
