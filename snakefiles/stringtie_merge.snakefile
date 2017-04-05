# Stringtie merge
rule stringtie_merge:
    input:
        gtf = expand(config['stringtie']['output_dir'] + "/{sample}.gtf", sample=SAMPLES)
    output:
        # gtf = config['stringtie']['output_dir'] + "/unannotated_lncRNA_primarymerge.gtf",
        gtf = config['stringtie']['output_dir'] + "/" + config['stringtie']['output_merge'],
    log:
        stderr = config['stringtie']['log_dir'] + "/stringtie-merge.log",
        version = config['version_dir'] + "/stringtie-version.txt",
    message: 
        "Executing stringtie on the following files {input}"
    shell:
        config['stringtie']['binary'] + 
        " --merge {input.gtf}"
        " -f " + config['stringtie']['min_iso'] +
        " -m " + config['stringtie']['min_length'] +
        " -G " + config['gtf_ref'] +
        " -o {output.gtf}" 
        # logs
        " 2> {log.stderr}"
        # stringtie version
        " && echo 'stringtie version $(stringtie --version)' > {log.version}"
