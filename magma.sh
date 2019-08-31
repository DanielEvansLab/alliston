#!/bin/bash -l

#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=10G
#SBATCH --time=9-00:00:00   
#SBATCH --output=magma_annot_bmd.out
#SBATCH --mail-user=danev1@gmail.com
#SBATCH --mail-type=ALL
#SBATCH --job-name="magma"
#SBATCH -p batch

cd $SLURM_SUBMIT_DIR

## Fracture
#annotation
#~/bigdata/programs/magma_1.07b/magma --annotate window=50 --snp-loc ../data/snp_loc_frac.txt --gene-loc ~/bigdata/programs/magma_1.07b/aux/gene_loc_b37/NCBI37.3.gene.loc --out ../results/frac

#gene analysis
#~/bigdata/programs/magma_1.07b/magma --bfile ~/bigdata/programs/magma_1.07b/ref/EUR/g1000_eur  --gene-annot ../results/frac.genes.annot --pval ../data/pvalfile_frac.txt ncol=N --gene-model snp-wise=top --out ../results/frac

#gene-set analysis
#~/bigdata/programs/magma_1.07b/magma --gene-results ../results/frac.genes.raw --set-annot ../data/serra_sets.txt --out ../results/frac

## BMD
#annotation
~/bigdata/programs/magma_1.07b/magma --annotate window=50 --snp-loc ../data/snp_loc_bmd.txt --gene-loc ~/bigdata/programs/magma_1.07b/aux/gene_loc_b37/NCBI37.3.gene.loc --out ../results/bmd

#gene analysis
~/bigdata/programs/magma_1.07b/magma --bfile ~/bigdata/programs/magma_1.07b/ref/EUR/g1000_eur  --gene-annot ../results/bmd.genes.annot --pval ../data/pvalfile_bmd.txt ncol=N --gene-model snp-wise=top --out ../results/bmd

#gene-set analysis
~/bigdata/programs/magma_1.07b/magma --gene-results ../results/bmd.genes.raw --set-annot ../data/serra_sets.txt --out ../results/bmd

