# Hisat alignment
rule hisat_pe:
    input:
        r1 = config['reads_dir']+"/{sample}_1.fastq.gz",
        r2 = config['reads_dir']+"/{sample}_2.fastq.gz",
    output:
        bam = config['hisat2']['output_dir']+"/{sample}.bam",
        bai = config['hisat2']['output_dir']+"/{sample}.bam.bai",
    threads: 
        config['nb_threads']
    log:
        stderr = config['hisat2']['log_dir'] + "/{sample}_hisat.log",
        version = config['report_versions'],
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
        # Samtools sort & index
        " samtools view -bS - |"
        " samtools sort -@4 - -o {output.bam}"
        " && samtools index {output.bam} {output.bai} "
        # crac version
        " && hisat2 --help | head -1 | cut -d' ' -f1-3 >> {log.version}"


"""
miseq () {
    # Aligner les reads Illumina
    shortname="$(basename ${mates1%R1.*}MiSeq)"
    hisat2 -p 10 -x $genome -1 $mates1 -2 $mates2 2> ${shortname}.log | \
    samtools view -bS - | \
    samtools sort - ${shortname} &&
    samtools index ${shortname}.bam ${shortname}.bai
}

"""
