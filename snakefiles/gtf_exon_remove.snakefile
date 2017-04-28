# gtf_exon_remove

"""
Après le subtract, il reste des lignes correspondant aux exons des gènes enlevés, 
le script perl de sébastient affiche les lignes de transcrit et seulement leurs 
exons correspondants, les exons seuls sont donc enlevés.

perl ~/gtf_exon_remove.pl unannotated_MSClinc_minusgenes.gtf > unannotated_MSClinc_v2.gtf

On a au final un gtf qui ne contient que ce qui ne croise pas quelque chose de connus !
"""

rule gtf_exon_remove:
    input:
        config['bedtools']['subtrack_gtf_file'],
    output:
        config['gtf_exon_remove']['output_file'],
    log:
        stderr = config['gtf_exon_remove']['log_file'],
        version = config['gtf_exon_remove']['version_file'],
    benchmark:
        config['gtf_exon_remove']['benchmark_file'],
    params:
        binary = config['gtf_exon_remove']['binary'],
        options = config['gtf_exon_remove']['options'],
    message: 
        "Executing gtf_exon_remove"
    shell:
        "{params.binary}"
        " {params.options}"
        " {input}"
        " > {output}"
        # logs
        " 2> {log.stderr}"
        # gtf_exon_remove version


