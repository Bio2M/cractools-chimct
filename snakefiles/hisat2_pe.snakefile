# Hisat alignment
rule hisat_pe:
    input:
        r1 = config['fastq_dir']+"/{sample}_1.fastq.gz",
        r2 = config['fastq_dir']+"/{sample}_2.fastq.gz",
    output:
        bam = config['hisat2']['output_dir']+"/{sample}.bam",
        bai = config['hisat2']['output_dir']+"/{sample}.bam.bai",
    threads: 
        config['nb_threads']
    log:
        stderr = config['hisat2']['log_dir'] + "/{sample}_hisat.log",
        version = config['hisat2']['version_file'],
    benchmark:
        "output/benchmarks/hisat/{sample}_hisat.benchmark"
    params:
        binary = config['hisat2']['binary'],
        options = config['hisat2']['options'],
        index = config['hisat2']['index'],
        samtools_threads = config['samtools']['threads'],
        samtools_mem = config ['samtools']['mem'],
        tmp = "/tmp/{sample}",
    message: 
        "Executing hisat on the following files {input.r1} {input.r2}"
    shell:
        "{params.binary}"
        " {params.options}"
        " -p {threads}"
        " -x {params.index}"
        " -1 {input.r1} -2 {input.r2}"
        # logs
        " 2> {log.stderr} |"
        # Samtools sort
        " samtools view -bS - |"
        " samtools sort "
        " -T {params.tmp} "
        " -@{params.samtools_threads}"
        " -m {params.samtools_mem}"
        " - -o {output.bam}"
        # index
        " && samtools index {output.bam} {output.bai} "
        # hisat version
        " && hisat2 --help | head -1 | cut -d' ' -f1-3 > {log.version}"

