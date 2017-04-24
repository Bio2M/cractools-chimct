# Extract common Chimeras
# COMMON_CHIMERAS_FILE = config['extract_common_chimeras']['output_dir'] + "/common_chimeras.txt"
# COMMON_CHIMERAS_SUMMARY_FILE = config['extract_common_chimeras']['output_dir'] + "/common_chimeras_summary.txt"
# BINARY = config['extract_common_chimeras']['binary']
# LOG = config['extract_common_chimeras']['log_file']


rule extract_common_chimeras:
    input:
        tsv = TSV_FILES,
    output:
        ext_common_chim = config['extract_common_chimeras']['output_file'],
    params:
        summary = config['extract_common_chimeras']['summary_file'],
    log:
        stderr = config['extract_common_chimeras']['log_file'],
        version = config['version_dir'] + "/extractCommonChimeras-version.txt",
    benchmark:
        "output/benchmarks/extract_common_chimeras.benchmark"
    message:
        "Executing extractCommonChimeras on {input.tsv}"
    shell:
        config['extract_common_chimeras']['binary'] +
        " --summary {params.summary}"
        " --output {output.ext_common_chim}"
        " {input.tsv}"
        # extractCommonChimeras.pl logs
        " 2>{log.stderr}"
        # version
        " && echo 'extractCommonChimeras.pl: no version available' > {log.version}"

    
"""
    extractCommonChimeras.pl    
    --summary Common_Chimeras_Summary.txt   
    chimct/om100011_chimct.tsv chimct/om110223_chimct.tsv > Common_Chimeras.txt
"""
