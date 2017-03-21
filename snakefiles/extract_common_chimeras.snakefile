# Extract common Chimeras
# COMMON_CHIMERAS_FILE = config['extract_common_chimeras']['output_dir'] + "/common_chimeras.txt"
# COMMON_CHIMERAS_SUMMARY_FILE = config['extract_common_chimeras']['output_dir'] + "/common_chimeras_summary.txt"
# BINARY = config['extract_common_chimeras']['binary']
# LOG = config['extract_common_chimeras']['log_file']


rule extract_common_chimeras:
    input:
        tsv = TSV_FILES,
    output:
        common_chim = config['extract_common_chimeras']['output_dir'] + "/common_chimeras.txt",
        summary = config['extract_common_chimeras']['output_dir'] + "/common_chimeras_summary.txt",
    log:
        stderr = config['extract_common_chimeras']['log_file'],
        version = config['report_versions'],
    message:
        "Executing extractCommonChimeras on {input.tsv}"
    shell:
        config['extract_common_chimeras']['binary'] +
        " --summary {output.summary}"
        " {input.tsv}"
        " > " + config['extract_common_chimeras']['output_dir'] + "/common_chimeras.txt"
        # extractCommonChimeras.pl logs
        " 2>{log.stderr}"
        # version
        " && echo 'extractCommonChimeras.pl: no version available' >>{log.version}"

    
"""
    extractCommonChimeras.pl    
    --summary Common_Chimeras_Summary.txt   
    chimct/om100011_chimct.tsv chimct/om110223_chimct.tsv > Common_Chimeras.txt
"""
