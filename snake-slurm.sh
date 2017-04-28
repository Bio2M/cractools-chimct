#!/bin/bash

{ time snakemake \
    -ap \
    -j 99 \
    --cluster-config cluster.yml \
    --cluster "sbatch \
        -A {cluster.account} \
        -p {cluster.partition} \
        -n {cluster.n}" ;
} 2>output/snakemake.log


# Move slurm output in logs directory

slurm_dir=output/logs/slurm

[ -d "$slurm_dir" ] || mkdir -p $slurm_dir
mv slurm-*.out $slurm_dir
