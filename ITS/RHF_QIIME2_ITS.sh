#!/bin/bash

# *** Script for QIIME2 Analysis of RHF Fungal Samples (ITS) ***

# Unzip all .fasta.qz files
# sudo gunzip -k *.fastq.gz
# Start QIIME
source activate qiime2-2020.11

# %%% Import fastq data %%%
qiime tools import \
 --input-path /home/akoontz11/bear-share/Pearse_seq/ITS/RHF_ITS_manifest.csv \
 --type SampleData[PairedEndSequencesWithQuality] \
 --output-path /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Demux/RHF_ITS.qza \
 --input-format PairedEndFastqManifestPhred33

# %%% Denoising with DADA2 %%%
# Last base (251) somewhat lower quality, so truncated by one; first sections of reads rather noisy, for first ~50 bases trimmed
qiime dada2 denoise-paired \
--i-demultiplexed-seqs /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Demux/RHF_ITS.qza \
--p-trunc-len-f 250 --p-trunc-len-r 250 \
--p-trim-left-f 52 --p-trim-left-r 61 \
--p-n-threads 24 \
--o-table /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseTable \
--o-representative-sequences /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseRepSeqs \
--o-denoising-stats /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseStats

# Generating FeatureTable and FeatureData visualizations
qiime feature-table summarize \
--i-table /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseTable.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseTable
qiime feature-table tabulate-seqs \
--i-data /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseRepSeqs.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseRepSeqs

# Clustering sequences de novo (no reference) at 99%
qiime vsearch cluster-features-de-novo \
--i-sequences /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseRepSeqs.qza \
--i-table /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Denoising/RHF_ITS_DenoiseTable.qza \
--p-perc-identity 0.99 \
--p-threads 16 \
--o-clustered-table /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringTable \
--o-clustered-sequences /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringSeqs

# Generating FeatureTable and FeatureData visualizations
qiime feature-table summarize \
--i-table /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringTable.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringTable
qiime feature-table tabulate-seqs \
--i-data /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringSeqs.qza \
--o-visualization /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringSeqs

# %%% Phylogenetic Tree for Diversity Analyses %%%
qiime phylogeny align-to-tree-mafft-raxml \
--i-sequences /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringSeqs.qza \
--p-n-threads 24 \
--o-alignment /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Phylogeny/RHF_ITS_AligedRepSeqs.qza \
--o-masked-alignment /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Phylogeny/RHF_ITS_MaskedAlignedRepSeqs.qza \
--o-tree /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Phylogeny/RHF_ITS_UnrootedTree.qza \
--o-rooted-tree /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Phylogeny/RHF_ITS_RootedTree.qza

# %%% Diversity Analyses %%%
qiime diversity core-metrics-phylogenetic \
--i-phylogeny /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Phylogeny/RHF_ITS_RootedTree.qza \
--i-table /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/Clustering/RHF_ITS_99ClusteringTable.qza \
--p-sampling-depth 34170 \
--p-n-jobs-or-threads 24 \
--m-metadata-file /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/ITS_RHF_metadata.csv \
--output-dir /home/akoontz11/bear-share/Pearse_seq/ITS/QIIME2/DiversityAnalyses

# Close out
conda deactivate
