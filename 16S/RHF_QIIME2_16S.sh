#!/bin/bash

# *** Script for QIIME2 Analysis of RHF Bacterial Samples (16S) ***

# Unzip all .fasta.qz files
# sudo gunzip -k *.fastq.gz
# Start QIIME
source activate qiime2-2020.11

# %%% Import fastq data %%%
 qiime tools import \
 --input-path /home/akoontz11/bear-share/Pearse_seq/16S/RHF_16S_manifest.csv \
 --type SampleData[PairedEndSequencesWithQuality] \
 --output-path /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Demux/RHF_16S.qza \
 --input-format PairedEndFastqManifestPhred33

# %%% Denoising with DADA2 %%%
qiime dada2 denoise-paired \
--i-demultiplexed-seqs /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Demux/RHF_16S.qza \
--p-trunc-len-f 250 --p-trunc-len-r 250 \
--p-trim-left-f 23 --p-trim-left-r 23 \
--p-n-threads 12 \
--o-table /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseTable \
--o-representative-sequences /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseRepSeqs \
--o-denoising-stats /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseStats

# Generating FeatureTable and FeatureData visualizations
qiime feature-table summarize \
--i-table /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseTable.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseTable
qiime feature-table tabulate-seqs \
--i-data /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseRepSeqs.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseRepSeqs

# %%% Clustering at 99% %%%
# Clustering sequences de novo (no reference) at 99%
qiime vsearch cluster-features-de-novo \
--i-sequences /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseRepSeqs.qza \
--i-table /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Denoising/RHF_16S_DenoiseTable.qza \
--p-perc-identity 0.99 \
--p-threads 24 \
--o-clustered-table /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringTable \
--o-clustered-sequences /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringSeqs

# Generating FeatureTable and FeatureData visualizations
qiime feature-table summarize \
--i-table /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringTable.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringTable
qiime feature-table tabulate-seqs \
--i-data /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringSeqs.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringSeqs

# %%% Phylogenetic Tree for Diversity Analyses %%%
qiime phylogeny align-to-tree-mafft-raxml \
--i-sequences /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringSeqs.qza \
--p-n-threads 24 \
--o-alignment /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Phylogeny/RHF_16S_AligedRepSeqs.qza \
--o-masked-alignment /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Phylogeny/RHF_16S_MaskedAlignedRepSeqs.qza \
--o-tree /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Phylogeny/RHF_16S_UnrootedTree.qza \
--o-rooted-tree /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Phylogeny/RHF_16S_RootedTree.qza

# %%% Diversity Analyses %%%
qiime diversity core-metrics-phylogenetic \
--i-phylogeny /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Phylogeny/RHF_16S_RootedTree.qza \
--i-table /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/Clustering/RHF_16S_99ClusteringTable.qza \
--p-sampling-depth 43363 \
--p-n-jobs-or-threads 24 \
--m-metadata-file /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/16S_RHF_metadata.csv \
--output-dir /home/akoontz11/bear-share/Pearse_seq/16S/QIIME2/DiversityAnalyses 

# Close out
conda deactivate
