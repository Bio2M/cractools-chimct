

# ChimCT
rule chimct:
    input:
        bam = BAM_DIR+"/{sample}.bam",
        gff = config['chimct']['gff'],
    output:
        tsv = config['chimct']['output_dir']+"/{sample}_chimct.tsv",
    log: 
        stderr = config['chimct']['log_dir'] + "/{sample}_chimct.log",
        version = config['version_dir'] + "/chimCT-version.txt",
    threads: 
        config['nb_threads']
    message:
        "Executing chimCT on {input.bam}",
    shell:
        config['chimct']['binary']
        + config['chimct']['options'] +
        " --gsnap-nb-threads {threads}"
        " --gsnap-genome-name " + config['chimct']['gsnap_genome_name'] +
        " --gsnap-genome-directory " + config['chimct']['gsnap_genome_dir'] +
        " --gsnap-exe " + config['chimct']['gsnap_binary'] + 
        " -n {wildcards.sample}"
        " -g " + ANNOTATION_FILE +
        " -s " + BAM_DIR + "/{wildcards.sample}.bam"
        " > " + config['chimct']['output_dir'] + "/{wildcards.sample}" + config['chimct']['output_suffix'] +
        # ChimCT logs
        " 2> {log.stderr}"
        # ChimCT version
        " && chimCT --version > {log.version}"

"""
chimCT  
    -n om100011 
    -g /data/annotations/human/ensembl.gff 
    -s mapping/crac/om100011.bam 
    --verify-splice 
    --gsnap-exe /data/bin/gsnap 
    --gsnap-genome-directory /data/indexes/gsnap/GRCh38 
    --gsnap-genome-name GRCh38 
    --detailed-sam 
    --gsnap-nb-threads 20 
    > chimct/om100011_chimct.tsv 
    2> chimct/om100011_chimct.log
"""
