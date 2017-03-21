# freebayes variant caller
rule freebayes:
    input:
        bam = config['hisat2']['output_dir'] + "/{sample}.bam"
    output:
        vcf = config['freebayes']['output_dir'] + "/{sample}.vcf"
    threads: 
        config['nb_threads']
    log:
        stderr = config['freebayes']['log_dir'] + "/{sample}_freebayes.log",
        version = config['report_versions'],
    message: 
        "Executing freebayes on the following files {input.bam}"
    shell:
        "freebayes"
        + config['freebayes']['options'] +
	" -t " + config['freebayes']['bed'] +
	" -f " + config['genome'] +
	" {input.bam}"
	" > {output.vcf}"
        # logs
        " 2> {log.stderr}"
        # version
        " && freebayes --version >> {log.version}"


"""
genome="/data/genomes/GRCh37/GRCh37_fused.fa"
# genome="/data/genomes/GRCh37/GRCh38.fa"
basedir="/data/nas/external/TP53"
bed="../manifest.bed"

options="-C 10 -F 0.008"

for (( i=0; $i<$nbfiles; i+=1 ))
do
    destfile=$(basename ${1%.*}.vcf)
    shortname=$(basename ${1%.*})
    [ -f $destfile ] && read -p "$destfile already exists, overwrite ? [Y,n] " over
    [ ${over:=y} == "y" ] && time ( freebayes $options -f $genome $1 > $destfile ) 2>${shortname}.log
    shift
done
"""
