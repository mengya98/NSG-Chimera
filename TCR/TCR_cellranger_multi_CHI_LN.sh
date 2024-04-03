#!/bin/bash
#SBATCH -p cn-long
#SBATCH -N 1
#SBATCH -c 20
#SBATCH -J LMY_10x
#SBATCH -o multi_cnl.out
#SBATCH -e multi_cnl.err
#SBATCH --no-requeue
#SBATCH -A tangfuchou_g1
#SBATCH --qos=tangfuchoucnl

cellranger=/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/liumengya/software/cellranger-7.1.0/cellranger

$cellranger multi \
    --id=CHI_LN \
    --csv=/gpfs1/tangfuchou_pkuhpc/tangfuchou_cls/liumengya/Immune_Rawdata/20230208_vdj/CHI_LN/File_CHI_LN.csv \
    --localcores=23 \
    --localmem=800
