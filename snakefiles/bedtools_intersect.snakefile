#

"""
J'ai essayé de créer la branche pour les lncRNA chevauchants antisens.

bedtools intersect \
    -­wa \
    ­-u \
    ­-S \ 
    -a unannotated_MSClnc_primarymerge.gtf \
    ­-b /data/annotations/human/Homo_sapiens.GRCh38.86.gtf > test
    
Devrait donner les transcripts reconstruits pour lequels il y a overlap antisens 
(meme sur une très petite partie) entre mes recontructions et les gènes connus.
Comme ça on devrait avoir des liste pour chaque catégories.
"""

rule bedtools_intersect:
    input:
        a = config['stringtie']['merge_gtf_file'],
        b = config['gtf_ref'],
    output:
        config['bedtools']['intersect_gtf_file'],
    log:
        stderr = config['bedtools']['intersect_log_file'],
        version = config['bedtools']['version_file'],
    benchmark:
        config['bedtools']['intersect_benchmark_file'],
    params:
        binary = config['bedtools']['binary'],
        options = config['bedtools']['intersect_options'],
    message: 
        "Executing bedtools on file {input.a}"
    shell:
        "{params.binary} intersect"
        " {params.options}"
        " -a {input.a}"
        " -b {input.b}"
        " > {output}"
        # logs
        " 2> {log.stderr}"
        # bedtools version
        " && bedtools --version > {log.version}"

