# tuto-snakemake-Bio2M

Lien vers la [documentation officielle de snakemake](https://snakemake.readthedocs.io/en/stable/)

#### Les étapes de création d'un pipeline


1. Initialiser un workflow
2. Créer un fichier Snakefile depuis un template
3. Créer un fichier samples.yml
4. Modifier le fichier config.yml
5. La commande Snakemake
6. slurm pour paralléliser les taches
7. Les fichiers créés




# 1. Initialiser un workflow

Un master permettant d'initialiser un workflow est présent sur le serveur gitlab

    git clone git@gitlab.montp.inserm.fr:benoit/snakemake-bio2m.git MyProject
    cd MyProject

Un certain nombre de fichiers sont présents

    $ tree -L 1
    ├── cluster.yml
    ├── config.yml
    ├── lib/
    ├── README.md
    ├── samples.yml
    ├── Snakefile
    ├── snakefiles/
    ├── snake-slurm.sh
    └── templates/

# 2. Créer un fichier Snakefile depuis un template

Le répertoire ``templates/`` contient la configuration pour un certain nombre de pipelines préconfigurés.

    $ ls templates/
    Snakefile.crac
    Snakefile.crac_chimct_ctextract_extcomchims
    Snakefile.hisat_freebayes
    Snakefile.hisat_freeb_extcomchim_chimct_cractools
    Snakefile.stringtie_bedtools_gtfExonRemove_gffread


### Utiliser un template existant
Il suffit de créer un lien symbolique

    ln -s templates/Snakefile.un_template Snakefile

Sinon, il faut créer un nouveau Snakefile, de préférence à partir d'un template existant


### Dans un template
    $ cat templates/Snakefile.stringtie_bedtools_gtfExonRemove_gffread

Quelques informations générales commentées

    """
    Author: B. Guibert - S. Riquier
    Affiliation: Bio2M
    Aim: projet LncRNA
    Date: May 2017
    Run: ./snake-slurm.sh

Suivies de la description "visuelle" du pipeline
    
    files.bam
        |--stringtie
            |--stringtie_merge
                |--bedtools_subtrack
                    |   |--gtf_exon_remove
                        |         |--gffread
                        |-- bedtools_intersect                                       
    """

Les instructions ``configfile`` chargent des fichiers de paramètes utiles pour notre pipeline

    # Config file
    configfile: "samples.yml"
    configfile: "config.yml


On peut passer des variables, en respectant les règles propres à python 3. Celles présentes ci-dessous exploitent des informations chargées grâce aux instructions ``configfile``.

    # variables from configfile
    SAMPLES = expand(config['samples'])
    BAM_DIR = config['raw_dir']
    EMAIL = config['email']

On charge des fichiers contenant des règles. Il y a un fichier par règle. Ceci permet la réutilisation des règles déjà créées et la maintenabilité.

    # includes
    include: 'snakefiles/stringtie.snakefile'
    include: 'snakefiles/bedtools_subtract.snakefile'
    include: 'snakefiles/gffread.snakefile'
    include: 'snakefiles/bedtools_intersect.snakefile'


La première règle du fichier, ici ``all`` permet à Snakemake de déduire les autres règles à appliquer.

    rule all:
        input:
            config['gffread']['output_file'],
            config['bedtools']['intersect_gtf_file'],
        log:
            version = "output/snakemake.version",
        shell:
             "echo snakemake versions $(snakemake --version) > {log.version}"


* snakemake va chercher la présence des fichiers dans le bloc ``input`` de la règle;
* s'ils ne sont pas présents, il les cherche dans les blocs ``output`` des autres règles ;
* si les fichiers sont définis dans le bloc ``output`` d'une règle, il vérifie la présence des fichiers définis en ``input`` de cette même règle ;
* et ainsi de suite jusqu'à trouver les fichiers en ``input`` d'une règle.




# 3. Créer un fichier samples.yml

Le fichier ``samples.yml`` contient :


* le répertoire contenant les données initiales (en général des fastq) ;
* la liste des échantillons.


    raw_dir:
        # For rule http_download, url: http://getdata.montp.inserm.fr/
        # 'getdata.montp.inserm.fr/Bio2M/external/sra/GSE62190/mapping/crac'
        input/raw
    samples:
        GSM1521606: ''
        GSM1521607: ''
        GSM1521608: ''

Le fichier ``samples.yml`` peut être créé manuellement ou généré grâce la commande ``lib/snaketools.py ``:

    $ lib/snaketools.py --help
    usage: snaketools.py [-h] [-o output_file] [--version] fastq_dir
                               
    positional arguments:
      fastq_dir             fastq directory
                                
    optional arguments:
      -h, --help            show this help message and exit
      -o output_file, --output_file output_file
                        yaml formatted output file, default : sample.yml
      --version, -v         show version


# 4. Modifier le fichier config.yml

Le fichier ``config.yml`` contient les paramètres généraux et les paramètres propre à chaque outil/commande.
Il faudra vérifier et modifier le cas échéant ces paramètres.

### Les paramètres généraux
    ################################################################################
    #   GLOBAL
    ################################################################################
    
    email:      'benoit.guibert@free.fr'
    nb_threads: 12
    genome:     '/data/genomes/GRCh38/GRCh38.fa'
    gtf_ref:    '/data/annotations/human/Homo_sapiens.GRCh38.88.gtf'


### Détail des paramètres pour la commande crac
    crac:
        binary:         'crac'
        index:          '/data/indexes/crac/GRCh38/GRCh38'
        kmer_len:       '22'
        output_dir:     'output/crac/bam'
        benchmark_dir:  'output/crac/benchmark'
        summary_dir:    'output/crac/summary'
        log_dir:        'output/crac/logs'
        version_file:   'output/crac/version/crac.version'
        options:        '--nb-tags-info-stored 10000 --stranded --detailed-sam'

**Nota** : bien noter que les paramètres sont liées aux commandes, plutôt qu'aux règles. par exemple, pour crac, il y aura deux règles, ``crac_pe`` et ``crac_se`` pour le paired-end et le single-end, mais il n'y a qu'une seule section ``crac`` comprenant les paramètres des deux règles.


# 5. La commande Snakemake

Toutes les options ne sont pas référencées ici, pour plus de détail, consulter la [doc en ligne](https://snakemake.readthedocs.io/en/stable/executable.html) de snakemake.

#### Lancer le pipeline
    $ snakemake

snakemake tente de trouver le fichier ``Snakefile`` depuis le répertoire courant, s'il ne le trouve pas, il renvoie une erreur

#### Tester, sans rien faire :
Très pratique pour tester son pipeline

    $ snakemake -np


* ``-n`` : tester sans rien faire
* ``-p`` : afficher la commande


#### Indiquer le nombre de cœurs à utiliser
    $ snakemake --cores 12

Cet argument est en corrélation avec l'instruction ``threads:`` définies dans les règles snakemake.
imaginons qu'une règle comporte l'instruction ``threads: 8``.
avec la commande ``snakemake`` sans argument, seul 1 cœur est utilisé par la règle.
avec la commande ``snakemake --cores 6``, 6 cœurs sont utilisés par la règle.
avec la commande ``snakemake --cores 12``, 8 cœurs sont utilisés par la règle.

#### Forcer l’exécution de tout le pipeline
    $ snakemake -p --forceall
    $ snakemake -np --forceall    


* ``-p`` : forcer l'exécution
* ``-np`` : simuler l'exécution forcée du pipeline


#### Après interruption forcée
Lorsque le précédent lancement de snakemake a été interrompu brutalement

    $ snakemake --unlock


#### Faire un graphique du pipeline
Le pipeline n'est pas exécuté dans ce cas (``-n`` implicite)

    $ snakemake --dag | dot -Tsvg > deg.svg

__Nota__ : créez un graphique avec très peu de samples, pour la lisibilité.

#### Forcer l’execution à un certain stade
Comme pour un --forceall, mais le pipeline sera exécuté à partir de la règle invoquée

    $ snakemake --forcerun une_regle


#### Utiliser un fichier Snakefile alternatif
    $ snakemake -s autreFichierSnakefile


# 6. slurm pour paralléliser les taches

### Lancer le pipeline sur le cluster de calcul
Le serveur **Marygold** est un serveur de calcul, comprenant 5 nœuds. Snakemake est capable d'en utiliser les fonctionalités.

Sur Marigold, on ne lance pas la commande snakemake, mais le script fournit ``snake-slurm.sh`` :

    $ snake-slurm.sh


### Les fichiers pour utiliser le serveur de calcul

* ``snake-slurm.sh`` : exécute snakemake en utilisant slurm
* ``cluster.yml`` :
    * applique des paramétrages par règle.
    * peut utiliser les mots génériques de snakemake (wildcards)


### Le fichier snake-slurm.sh
Il suffit d'exécuter le fichier ``snake-slurm.sh`` pour lancer le workflow.

    #!/bin/bash
    # loading environment
    [ -f /data/bin/init.sh ] && . /data/bin/init.sh 
                                
    { time snakemake -ap -j 99 \
        --cluster-config cluster.yml \
        --cluster "sbatch \
        -A {cluster.account} \
        -p {cluster.partition} \
        -J {cluster.jobname} \
        -e $slurm_log_dir/{cluster.output} \
        -o /dev/null \
        -n {cluster.n}" ;
    } 2>output/snakemake.log  

Sauf pour un usage particulier, il n'est pas nécessaire de le modifier

### Le fichier cluster.yml
Le fichier ``cluster.yml`` permet de réutiliser des mots génériques de snakemake, par exemple ``wildcards.sample`` pour les utiliser avec slurm. Il permet également d'utiliser des paramètres spécifiques pour certaines règles.

    __default__ :
        account:    "benoit"
        partition:  "std"
        jobname:    "{wildcards.sample}"
        output:     "{wildcards.sample}_%j.log"
        n :         10
    
    histat2_pe :
        n : 16

la section ``__default__`` s'applique de manière globale, mais les paramètres peuvent être spécifique pour certaines règles.

Dans l'exemple ci-dessus, les paramètres par défaut sont :

* Le nom login l'utilisateur est ``benoit`` ;
* la partition de cluster utilisée est ``std`` ;
* le nom du job affiché sera celui de l'échantillon (``wildcards.sample``) ;
* le nom du fichier de log sera le nom de l'échantillon suivi du numéro de job (``%j``) ;
* le nombre de cœurs maximum par nœud sera de 10.


Mais pour la règle ``hisat2_pe``, le nombre maxi de cœurs sera de 16.


# 7. Les fichiers créés

Par défaut, tous les fichiers générés par le workflow sont stockés dans un répertoire ``output/``.

Par exemple : 

    $ tree -L 1 output/output/
    ├── chimct
    ├── crac
    ├── cractools
    ├── logs
    └── snakemake.log


### Le fichier snakemake.log
contient le rapport renvoyé par snakemake. On y trouve : 

* des informations générales telles que le nombre de tâches exécutées ;
* le détail de chaque tâche (avec la commande complète) ;
* l'état d'avancement de tâche par tâche ;
* le temps d'exécution total (par la commande time) si le pipeline s'est achevé correctement ;
* éventuellement, on y trouve les erreurs produites par snakemake.


### Le répertoire logs
Contient les logs liés aux outils utilisés pour généré le pipeline, à l'exeption du fichier ``snakemake.log``.
Actuellement, il contient le sous-répertoire ``slurm/``.

#### Le sous-répertoire logs/slurm/
Contient un fichier par tâche exécutée, générée par slurm.

#### Les sous-répertoires faisant référence à des commandes
Dans notre exemple, les sous répertoires ``chimct``, ``crac`` et ``cractools`` représentent chacun les commandes utilisées par le pipeline. 
En plus du ou des répertoires contenant les fichiers de résultats de la commande, on y trouvera d'autre répertoires.
Par exemple, voici le contenu du répertoire ``chimct`` : 

* ``tsv`` : répertoire propre à la commande, contient les fichiers de sortie ;
* ``summary`` : un autre répertoire propre à la commande, qui contient ici des résumés des résultats obtenus ;
* ``benchmark`` : contient la durée prise pour l'exécution de la tâche, en secondes et en heures/minutes/secondes ;
* ``logs`` : contient les informations normalement renvoyées à l'écran par la commande utilisée. Une taille différente de quelques d'un ou quelques uns des fichiers peut vouloir dire que des erreurs ont été consignées dans le fichier ;
* ``version`` : contient la version de la commande utilisée.



