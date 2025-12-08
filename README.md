# ArcGlacierPRJ
Scripts and pipelines for the Arctic Glacier Microbiome Project (ArcGlacierPRJ)

## Metagenomic Analysis Pipeline
End-to-end bioinformatics pipeline for processing metagenomic data from Arctic glacier samples, encompassing the complete workflow:
- Raw data quality control and adapter/trimming (TrimGalore)
- Metagenomic assembly (MEGAHIT)
- Open Reading Frame (ORF) prediction (Prodigal)
- Metagenome-Assembled Genome (MAG) binning using multiple tools:
  - MetaBAT2
  - MaxBin2
  - CONCOCT
  - VAMB
- DAS Tool integration for consolidated MAG binning results
- MAG quality assessment and dereplication:
  - CheckM2 (quality evaluation)
  - dRep (dereplication to remove redundant MAGs)
- MAG taxonomic classification (GTDB-Tk)
- Coverage calculation for MAGs (CoverM)
- Functional annotation of MAGs:
  - KEGG pathway annotation
  - Metabolic pathway reconstruction

## Network Analysis Pipeline
Pipeline for constructing and analyzing co-occurrence networks from MAG abundance data using sparkCC:
- FastSPRCC calculation (spearman correlation for microbial co-occurrence)
- Network construction and visualization in R (e.g., igraph, ggplot2)
- Network topology analysis (degree, modularity, centrality metrics)

