#!/bin/bash

#load proxy to enable comunication of the node with the internet
module load eth_proxy

#load R and needed functionalities
module load stack/2024-06 r/4.4.0 cmake/3.27.7 udunits/2.2.28 openssl/3.1.3-zhfub4o gdal/3.4.3 geos/3.9.1 proj/9.2.1 sqlite/3.43.2

# srun --ntasks=1 --cpus-per-task=20 --time=4:00:00 --mem-per-cpu=1024 --x11 --pty bash
# sbatch --ntasks=1 --cpus-per-task=20 --time=4:00:00 --mem-per-cpu=1024 Rscript development_scripts/cluster_comp/cluster_setup.sh                                                                                                                         │** testing if installed package can be loaded from final location
Rscript Rscript /cluster/home/merler/USE.MCMC/development_scripts/cluster_comp/distributed_precomputation.R --seed 42
