# http_download

# To get files throught HTTP
from snakemake.remote.HTTP import RemoteProvider as HTTPRemoteProvider
HTTP = HTTPRemoteProvider()


rule http_download_fastq:
    input:
        R1 = HTTP.remote(config['http']['input_dir'] + "/{sample}_01.fastq.gz", keep_local=config['http']['keep_local']),
        R2 = HTTP.remote(config['http']['input_dir'] + "/{sample}_02.fastq.gz", keep_local=config['http']['keep_local']),
    output:
        R1 = temp(config['raw_dir'] + "/{sample}_01.fastq.gz"),
        R2 = temp(config['raw_dir'] + "/{sample}_02.fastq.gz"),
    shell:
        "mv {input.R1} {output.R1} ;"
        "mv {input.R2} {output.R2}"


rule http_download_bam:
    input:
        bam = HTTP.remote(config['http']['input_dir'] + "/{sample}.bam", keep_local=config['http']['keep_local']),
        bai = HTTP.remote(config['http']['input_dir'] + "/{sample}.bam.bai", keep_local=config['http']['keep_local']),
    output:
        bam = temp(config['raw_dir'] + "/{sample}.bam"),
        bai = temp(config['raw_dir'] + "/{sample}.bam.bai"),
    shell:
        "mv {input.bam} {output.bam} ;"
        "mv {input.bai} {output.bai}"
