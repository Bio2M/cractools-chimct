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
        a = config['stringtie']['output_dir'] + "/unannotated_lncRNA_primarymerge.gtf",
        b = config['gtf_ref'],
    output:
        config['bedtools']['output_intersect'],
    log:
        stderr = config['bedtools']['log_dir'] + "/bedtools-intersect.err",
        version = config['version_dir'] + "/bedtools-version.txt",
    message: 
        "Executing bedtools on file unannotated_lncRNA_primarymerge"
    shell:
        config['bedtools']['binary'] + 
        " intersect "
        + config['bedtools']['intersect_options'] +
        " -a {input.a}"
        " -b {input.b}"
        " > {output}"
        # logs
        " 2> {log.stderr}"
        # bedtools version
        " && bedtools --version > {log.version}"

