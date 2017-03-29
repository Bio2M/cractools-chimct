#!/bin/bash

{ time snakemake \
    -j 99 \
    --cluster-config cluster.yml \
    --cluster "sbatch \
        -A {cluster.account} \
        -p {cluster.partition} \
        -e {log_dir}/slurm-%A.err \
        -o {log_dir}/slurm-%A.out \
        -n {cluster.n}" ;
} 2> log.err >log.out
