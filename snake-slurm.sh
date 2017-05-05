#!/bin/bash

# loading environment
[ -f /data/bin/init.sh ] && . /data/bin/init.sh 

# slurm log directory
slurm_log_dir=output/logs/slurm
[ -d "$slurm_log_dir" ] || mkdir -p $slurm_dir


{ time snakemake \
    -ap \
    -j 99 \
    --cluster-config cluster.yml \
    --cluster 'sbatch \
        -A {cluster.account} \
        -p {cluster.partition} \
        -J {cluster.jobname} \
        -e $slurm_log_dir/{cluster.output} \
        -o /dev/null \
        -n {cluster.n}' ;
} 2>output/snakemake.log

