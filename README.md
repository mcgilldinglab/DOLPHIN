# DOLPHIN
**D**eep Ex**o**n-**l**evel Gra**ph** Neural Network for S**i**ngle-cell Representatio**n** Learning and Alternative Splicing

## Overview
<img title="DeepExonas Overview" alt="Alt text" src="DeepExonas_overview.png">
The advent of single-cell sequencing has revolutionized the study of cellular dynamics, providing unprecedented resolution into the molecular states and heterogeneity of individual cells. However, the rich potential of exon-level information and junction reads within single cells remains underutilized. Conventional gene-count methods overlook critical exon and junction data, limiting the quality of cell representation and downstream analyses such as subpopulation identification and alternative splicing detection. To address this, we introduce DeepExonas, a deep learning method that integrates exon-level and junction read data, representing genes as graph structures. These graphs are processed by a variational autoencoder to improve cell embeddings. Compared to conventional gene-based methods, DeepExonas shows superior performance in cell clustering, biomarker discovery, and alternative splicing detection, providing deeper insights into cellular processes. By examining cellular dynamics with enhanced resolution, DeepExonas detects subtle differences often missed at the gene level, offering new insights into disease mechanisms and potential therapeutic targets.

## Key Capabilities of DeepExonas:

- Exon-Level Quantification: It represents genes as graphs, where nodes are exons and edges are junction reads, capturing detailed transcriptomic information at the exon level.
- Better Cell Embedding: DeepExonas leverages exon and junction read data to significantly improve the accuracy of cell embeddings, providing better resolution and resulting in more precise, biologically meaningful cell clusters compared to conventional gene-count based approaches.
- Enhanced Alternative Splicing Detection: By aggregating exon and junction reads from neighboring cells, DeepExonas significantly enhances the detection of alternative splicing events, providing deeper insights into cell-specific splicing patterns.
- Superior Performance in Downstream Analysis: DeepExonas consistently outperforms conventional gene-count methods in multiple downstream tasks, including the identification of differential exon markers and alternative splicing events. This high-resolution approach allows DeepExonas to uncover biologically significant exon markers that are often missed by traditional methods.

## Installation

Create a new conda environment
```
conda create -n deepexonas_Env python=3.9
conda activate deepexonas_Env
```

### Install from Github

Installing DeepExonas directly from GitHub ensures you have the latest version. 
**(Please install directly from GitHub to use the provided Jupyter notebooks for tutorials and walkthrough examples.)**

```
git clone https://github.com/mcgilldinglab/DeepExonas.git
cd DeepExonas
pip install .
```

## Tutorials:

### Dataset preparation
Please refer to the step0_preprocess.ipynb file in the [tutorial](https://github.com/mcgilldinglab/DeepExonas/tree/main/tutorial).
For the processed dataset, please download it from this link and save it to the test_dataset folder.
[Prepare datasets to run DeepExonas.](https://mcgill-my.sharepoint.com/:f:/g/personal/kailu_song_mail_mcgill_ca/EvZtHeW7qjJJs_RHc2-327ABeLXafa-ruvfk9Vs134crig?e=jNygC6)

### Training on an example dataset
```
cd tutorial
bash step1_run_model.sh
```
