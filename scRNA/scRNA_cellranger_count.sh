#!/bin/bash
#SBATCH -p fat4way
#SBATCH -N 1 
#SBATCH -c 24
#SBATCH -J LMY_10x_ln
#SBATCH -o ln.out
#SBATCH -e ln.err
#SBATCH --no-requeue
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchouf4w

/lustre1/tangfuchou_pkuhpc/Software/cellranger/cellranger-3.0.2/cellranger-cs/3.0.2/bin/cellranger count \
    --id=LN_outs \
    --transcriptome=/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/liumengya/Pipeline/database/refdata-cellranger-mm10-3.0.0 \
    --fastqs=/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/liumengya/Immune_Rawdata/CHI_LN/LN \
    --sample=LN \
    --expect-cells=7000 \
    --localcores=23 \
    --localmem=800

