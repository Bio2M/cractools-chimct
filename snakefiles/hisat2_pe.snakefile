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
        version = config['version_dir'] + "/hisat-version.txt",
    message: 
        "Executing hisat on the following files {input}"
    shell:
        "hisat2 "
        + config['hisat2']['options'] +
        " -p {threads}"
        " -x " + config['hisat2']['index'] +
        " -1 {input.r1} -2 {input.r2}"
        # logs
        " 2> {log.stderr} |"
        # Samtools sort
        " samtools view -bS - |"
        " samtools sort "
        " -@" + config['samtools']['threads'] +
        " -m " + config ['samtools']['mem'] +
        " - -o {output.bam}"
        # index
        " && samtools index {output.bam} {output.bai} "
        # hisat version
        " && hisat2 --help | head -1 | cut -d' ' -f1-3 > {log.version}"

