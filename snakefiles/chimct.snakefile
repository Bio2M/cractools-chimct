

# ChimCT
rule chimct:
    input:
        bam = BAM_DIR + "/{sample}.bam",
        gff = config['chimct']['gff'],
    output:
        tsv = config['chimct']['output_dir']+"/{sample}" + config['chimct']['output_suffix'],
        summary = config['chimct']['summary_dir'] + "/{sample}-chimct.summary",
    threads: 
        config['nb_threads'],
    log: 
        stderr = config['chimct']['log_dir'] + "/{sample}-chimct.log",
        version = config['version_dir'] + "/chimCT-version.txt",
    benchmark:
        "benchmarks/chimct/{sample}_chimct.benchmark"
    message:
        "Executing chimCT on {input.bam}",
    shell:
        config['chimct']['binary'] + " "
        + config['chimct']['options'] +
        " --tmp-dir " + config['chimct']['tmp_dir'] +
        " --summary={output.summary}"
        " --gsnap-nb-threads {threads}"
        " --gsnap-genome-name " + config['chimct']['gsnap_genome_name'] +
        " --gsnap-genome-directory " + config['chimct']['gsnap_genome_dir'] +
        " --gsnap-exe " + config['chimct']['gsnap_binary'] + 
        " -n {wildcards.sample}"
        " -g " + config['chimct']['gff'] +
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
