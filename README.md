# Utiliser Snakemake

Lien vers la [documentation officielle](https://snakemake.readthedocs.io/en/stable/)


# Préparer snakemake

git clone git@gitlab.montp.inserm.fr:benoit/snakemake-bio2m.git  monSnake

On se retrouve avec l’arborescence suivante : 

.  
├── lib/  
├── snakefiles/  
├── templates/  
├── cluster.yml  
├── config.yml  
├── samples.yml  
├── Snakefile -> templates/Snakefile.crac-chimct-ctextract-extcomchims  
└── snake-slurm.sh  

__lib :__

Pour y mettre du code à nous
snakefiles : 
contient des fichiers Snakefile, chacun avec une "règle". Ces fichiers seront invoqués par des instructions includes du fichier Snakefile "maître". En usage courant, on ne modifie pas ces fichiers

__Templates :__

sont stockées des fichier Snakefiles "maîtres", que l’on peut réutiliser tel quels ou en les modifiant. Ils chargent et invoquent des règles contenues dans le répertoire snakefiles.

__Cluster.yml :__
 
Liés à slurm, ce fichier permet de récupérer des expressions snakemake pour les réutiliser dans slurm. Il permet aussi d’avoir des paramétrages slurm en fonction des règles snakemake (par exemple, la règle http_download n’utilisera que 4 nœuds maximum).
Config.yml :
fichier au format yaml, il contient les paramètres à utiliser pour son pipeline, règle par règle. Par exemple, on y défini le chemin des fichiers fastq, le nombre de threads à utiliser, les options pour les différentes commandes, etc. Ce fichier sera systématiquement a adapter en fonction de ses besoins.

__Samples.yml :__

Un autre fichier au format yaml, contenant le nom des échantillons, sans le répertoire ni l’extension, et pour les fichiers fastq pairés sans le numéro de paire ni le signe de liaison (undersore souvent)

__Snakefile :__

lien symbolique vers un des fichiers du répertoire template. On peut aussi créer sont propre fichier "Snakefile". On utilise le template le plus proche de ses besoins, en le modifiant au besoin. 

__Snake-slurm.sh :_

script shell invoquant snakemake dans un environnement de cluster de calcul slurm. En général, il n’est pas nécessaire de la modifier.


La commande snakemake

Tester, sans rien faire : 

```
Snakemake -np
```
    -n : tester sans rien faire
    -p : afficher la commande


Forcer l’exécution de tout le pipeline (en réel ou simulation)

```
snakemake -p --forceall
snakemake -np --forceall
```

Après interruption forcée

```
snakemake --unlock
```

faire un joli graphique du pipeline (pas d’éxecution du pipeline)

```
snakemake --dag | dot -Tsvg > deg.svg
```

Forcer l’execution à un certain stade

```
snakemake --forcerun <nom_de_la_regle>
```


# Le fichier Snakefile

```
# Global snakefile

# Config file
configfile: "samples.yml"
configfile: "config.yml"


# variables from configfile
SAMPLES = expand(config['samples'])

# bam parameters
BAM_DIR = config['bam_dir']

# others
EMAIL = config['email']

# includes
include: 'snakefiles/stringtie.snakefile'
include: 'snakefiles/stringtie_merge.snakefile'
include: 'snakefiles/bedtools_subtract.snakefile'
include: 'snakefiles/gtf_exon_remove.snakefile'
include: 'snakefiles/gffread.snakefile'
include: 'snakefiles/bedtools_intersect.snakefile'


rule all:
    input:
        # expand(config['stringtie']['output_dir'] + "/{sample}_stringtie.gtf", sample=config['samples']),
        # config['stringtie']['merge_gtf_file'],
        # config['bedtools']['subtract_gtf_file'],
        # config['gtf_exon_remove']['output_file'],
        config['gffread']['output_file'],
        config['bedtools']['intersect_gtf_file'],
    log:
        version = "output/snakemake.version",
    shell:
        "echo snakemake versions $(snakemake --version) > {log.version}"


onsuccess:
    print("Workflow finished, no error")
    #shell( 'mail -s "$HOSTNAME: Workflow finished, no error" ' + EMAIL + ' < {log}')

onerror:
    print("An error occured, pas glop!")
    #shell('mail -s "$HOSTNAME: an error occurred on Workflow" ' + EMAIL + ' < {log}')
```

# Le fichier config.yml

A modifier en fonction de ses besoins.
