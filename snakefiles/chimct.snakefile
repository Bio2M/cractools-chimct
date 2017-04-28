# ChimCT
rule chimct:
    input:
        bam = BAM_DIR + "/{sample}.bam",
        gff_ref = config['chimct']['gff_ref_file'],
    output:
        tsv = config['chimct']['output_dir'] + "/{sample}" + config['chimct']['output_suffix'],
        summary = config['chimct']['summary_dir'] + "/{sample}" + config['chimct']['summary_suffix'],
    threads: 
        config['nb_threads'],
    log: 
        stderr = config['chimct']['log_dir'] + "/{sample}_chimct.log",
        version = config['chimct']['version_file'],
    benchmark:
        config['chimct']['benchmark_dir'] + "/{sample}" + config['chimct']['benchmark_suffix']
    params:
        binary = config['chimct']['binary'],
        options = config['chimct']['options'],
        tmp_dir = config['chimct']['tmp_dir'],
        gsnap_genome_name = config['chimct']['gsnap_genome_name'],
        gsnap_genome_dir = config['chimct']['gsnap_genome_dir'],
        gsnap_binary = config['chimct']['gsnap_binary'],
    message:
        "Executing chimCT on {input.bam}",
    shell:
        "{params.binary}"
        " {params.options}"
        " --tmp-dir {params.tmp_dir}"
        " --summary={output.summary}"
        " --gsnap-nb-threads {threads}"
        " --gsnap-genome-name {params.gsnap_genome_name}"
        " --gsnap-genome-directory {params.gsnap_genome_dir}"
        " --gsnap-exe {params.gsnap_binary}"
        " -n {wildcards.sample}"
        " -g {input.gff_ref}"
        " -s {input.bam}"
        " > {output.tsv}"
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
