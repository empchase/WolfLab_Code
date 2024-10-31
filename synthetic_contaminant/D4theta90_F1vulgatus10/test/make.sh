#!/usr/bin/env bash
#Author: empchase@berkeley.edu

#SBATCH --job-name=pear
#
# Account:
#SBATCH --account=ac_wolflab
# Partition:
#SBATCH --partition=savio2
#
# Wall Clock Limit:
#SBATCH --time=01:00:00
## Commands to run:

make --jobs $(nproc)