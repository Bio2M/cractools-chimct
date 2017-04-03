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
        config['bedtools']['output_subtrack'],
    output:
        config['gtf_exon_remove']['output_file'],
    log:
        stderr = config['gtf_exon_remove']['log_file'],
        version = config['version_dir'] + "/gtf_exon_remove-version.txt",
    message: 
        "Executing gtf_exon_remove"
    shell:
        config['gtf_exon_remove']['binary'] + 
        config['gtf_exon_remove']['options'] +
        " {input}"
        " > {output}"
        # logs
        " 2> {log.stderr} |"
        # gtf_exon_remove version


