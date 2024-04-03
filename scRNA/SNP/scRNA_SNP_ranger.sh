#!/bin/bash
#SBATCH -p fat4way
#SBATCH -N 1 
#SBATCH -c 20
#SBATCH -J LMY_10x_ln
#SBATCH -o ln.out
#SBATCH -e ln.err
#SBATCH --no-requeue
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchouf4w

export PATH=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/wangrui/Software/NewInstall_Software/Ranger/ranger-1.0.1:/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/wangrui/Software/NewInstall_Software/Ranger/ranger-1.0.1/freebayes/v1.0.2:$PATH
export MROPATH=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/wangrui/Software/NewInstall_Software/single-cell-3prime-snp-clustering/mro:$MROPATH
export PYTHONPATH=/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/wangrui/Software/NewInstall_Software/single-cell-3prime-snp-clustering/lib/python/:/gpfs1/tangfuchou_pkuhpc/tangfuchou_test/wangrui/Software/NewInstall_Software/Ranger/ranger-1.0.1/ranger-cs/1.0.1/tenkit/lib/python:$PYTHONPATH

ranger mrp ln_snp.mro ln_snp
