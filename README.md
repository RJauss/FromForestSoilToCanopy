# From forest soil to the canopy: Increased habitat diversity does not increase species richness of Cercozoa and Oomycota in tree canopies

Welcome to the **From Forest Soil to the Canopy** repository!

This repository is a collection of several scripts and mini-tutorials guiding you through the methods of metabarcoding analyses which were performed in the paper by [Jauss & Walden et al., 2020]().

Most scripts here deal with the data for oomycetes, but they were also applied for Cercozoa if not stated otherwise. The raw data can be downloaded [here](), plots and figures were generated with the final OTU tables and annotation files accessible in the folder [00_Data](00_Data/). 

## Table of Content
### 00 Data
This folder contains the final OTU table, taxonomic annotation, sample metadata as well as the oligosheet for demultiplexing the raw data. Intermediate files are not provided here, but you can generate them yourself by following the next steps.

### 01 Metabarcoding Pipeline
**[In this pipeline](01_Metabarcoding-Pipeline/Metabarcoding-Pipeline.md)**, you find the neccessary scripts to generate the (unfiltered) OTU table from raw .fastq files.

### 02 Taxonomic Annotation and Visualisation
This is a collection of several scripts neccessary for the taxonomic annotation of our OTUs. What we do first is to **[BLAST against the NCBI nt database to remove contaminants](02_Taxonomic_Annotation_and_Visualisation/BLAST-against-NCBI-nt-Database.md)**. After that, we **[download and process public oomycete/cercozoan sequences](02_Taxonomic_Annotation_and_Visualisation/Downloading-&-Processing-ITS-Sequences.md)**, which we use as a reference database for our **[taxonomic annotation with `vsearch`](02_Taxonomic_Annotation_and_Visualisation/Annotate-with-vsearch-and-the-ITS1-reference-database.md)**.

The visualisation of the taxonomy then includes **[plotting the sequence similarity to reference sequences]()** and a diagram showing the **[total taxonomic composition]()**.

### 03 Postprocessing the OTU Table
This section provides scripts on how to **[import, explore and filter the OTU table with `Qiime`](03_Postprocessing_OTU-Table/Importing-and-Filtering-OTU-Table.md)**, how to **[extract sequences from the filtered table](03_Postprocessing_OTU-Table/Postprocessing-the-OTU-Table.md#Extract-Sequences-from-Filtered-Table)** and last but not least how to **[paste the filtered metadata into the filtered table](03_Postprocessing_OTU-Table/Postprocessing-the-OTU-Table.md#Paste-Filtered-OTU-Table-and-Filtered-Metadata)**.

### 04 Determining Alpha Diversity
Here we deal with the methods of **[plotting rarefaction curves](04_Alpha_Diversity/)** and of course how to **[plot alpha diversity indices in a boxplot](04_Alpha_Diversity/)**.

### 05 Exploring Beta Diversity
One of the most straightforward methods of visualising beta diversity is an NMDS plot, the script is provided **[here](04_Beta_Diversity/)**. But we also perform and plot a **[redundancy analysis](04_Beta_Diversity/)**. 

### 06 Incidence based diversity
This section deals with the **[visualisation of shared OTUs](06_Incidence-based_Diversity/)** and how to plot a **[species accumulation curve](06_Incidence-based_Diversity/)**. Furthermore, we explore the OTU table to **[check for nestedness](06_Incidence-based_Diversity/)** in our OTU distribution.

