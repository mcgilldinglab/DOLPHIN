<img title="DOLPHIN Logo" alt="Alt text" src="DOLPHIN_logo.png">

## Overview
<img title="DOLPHIN Overview" alt="Alt text" src="Overview_DOLPHIN.png">
The advent of single-cell sequencing has revolutionized the study of cellular dynamics, providing unprecedented resolution into the molecular states and heterogeneity of individual cells. However, the rich potential of exon-level information and junction reads within single cells remains underutilized. Conventional gene-count methods overlook critical exon and junction data, limiting the quality of cell representation and downstream analyses such as subpopulation identification and alternative splicing detection. To address this, we introduce DOLPHIN, a deep learning method that integrates exon-level and junction read data, representing genes as graph structures. These graphs are processed by a variational autoencoder to improve cell embeddings. Compared to conventional gene-based methods, DOLPHIN shows superior performance in cell clustering, biomarker discovery, and alternative splicing detection, providing deeper insights into cellular processes. By examining cellular dynamics with enhanced resolution, DOLPHIN detects subtle differences often missed at the gene level, offering new insights into disease mechanisms and potential therapeutic targets.

## Key Capabilities of DOLPHIN:

- Exon-Level Quantification: It represents genes as graphs, where nodes are exons and edges are junction reads, capturing detailed transcriptomic information at the exon level.
- Better Cell Embedding: DOLPHIN leverages exon and junction read data to significantly improve the accuracy of cell embeddings, providing better resolution and resulting in more precise, biologically meaningful cell clusters compared to conventional gene-count based approaches.
- Enhanced Alternative Splicing Detection: By aggregating exon and junction reads from neighboring cells, DOLPHIN significantly enhances the detection of alternative splicing events, providing deeper insights into cell-specific splicing patterns.
- Superior Performance in Downstream Analysis: DOLPHIN consistently outperforms conventional gene-count methods in multiple downstream tasks, including the identification of differential exon markers and alternative splicing events. This high-resolution approach allows DOLPHIN to uncover biologically significant exon markers that are often missed by traditional methods.

## Installation

Installing DOLPHIN directly from GitHub ensures you have the latest version. 
**(Please install directly from GitHub to use the provided Jupyter notebooks for tutorials)**

```
git clone https://github.com/mcgilldinglab/DOLPHIN.git
cd DOLPHIN
```

Creating and Activating the Conda Environment

```
conda env create -f environment.yaml
conda activate DOLPHIN
```

Installing the DOLPHIN Package
1. Standard Installation
```
pip install .
```

2. Developer Mode Installation
```
pip install -e .
```

Validate That DOLPHIN Is Successfully Installed
```
import DOLPHIN
```

## Tutorials:

### Dataset Preparation

1. First, generate the exon-level reference GTF file by following the instructions in the [exon_gtf_generation](https://github.com/mcgilldinglab/DOLPHIN/blob/main/tutorial/step0_generate_exon_gtf.ipynb) tutorial.

2. Then, use the following tutorials to align the raw RNA-seq data and generate exon read counts and junction read counts:

   - For **Full-length scRNA-seq**, refer to the [Full-length scRNA-seq tutorial](https://github.com/mcgilldinglab/DOLPHIN/blob/main/tutorial/step1_1_preprocess_full_length.md).

   - For **10X RNA-seq**, refer to the [10X tutorial](https://github.com/mcgilldinglab/DOLPHIN/blob/main/tutorial/step1_2_preprocess_10X.md).


### Model Training and Cell Embedding Visualization
[DOLPHIN Training and Cell Embedding](https://github.com/mcgilldinglab/DOLPHIN/blob/main/tutorial/run_DOLPHIN.ipynb)

#### Run on example dataset:
You can download the processed dataset from [here](https://mcgill-my.sharepoint.com/my?id=%2Fpersonal%2Fkailu%5Fsong%5Fmail%5Fmcgill%5Fca%2FDocuments%2FDeepExonas%5Fgithub%5Fexample%2Fprocessed%5Fdataset)
and follow the [example](https://github.com/mcgilldinglab/DOLPHIN/blob/main/example/run_DOLPHIN.ipynb) to run the model.

