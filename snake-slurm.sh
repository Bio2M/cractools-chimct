#!/bin/bash

snakemake \
	-np \
	-j 99 \
	--cluster-config cluster.yml \
	--cluster "sbatch \
        -A {cluster.account} \
        -p {cluster.partition} \
    	-n {cluster.n}"
