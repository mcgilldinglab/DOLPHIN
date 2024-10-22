# DeepExonas
Advancing Cell Representation Beyond Gene-Level by Integrating Exon-Level Quantification and Junction Reads with Deep Neural Networks.

## Overview
<img title="DeepExonas Overview" alt="Alt text" src="DeepExonas_overview.png">
DeepExonas is an advanced deep learning framework designed for exon-level quantification and alternative splicing detection in single-cell RNA sequencing (scRNA-seq) data. Traditional gene-count methods often overlook critical exon and junction read information, missing intricate splicing mechanisms crucial for understanding complex cellular processes. DeepExonas addresses these limitations by integrating exon read counts and junction reads into a graph-based neural network model, allowing for more accurate cell representation and alternative splicing analysis.

## Key Capabilities of DeepExonas:

- Exon-Level Quantification: It represents genes as graphs, where nodes are exons and edges are junction reads, capturing detailed transcriptomic information at the exon level.
- Enhanced Cell Clustering: The integration of exon and junction read data leads to improved cell clustering accuracy compared to gene-count based methods.
- Alternative Splicing Detection: By aggregating exon and junction reads from neighboring cells, DeepExonas significantly enhances the detection of alternative splicing events, providing deeper insights into cell-specific splicing patterns.
- Superior Performance: DeepExonas consistently outperforms conventional gene-count methods in identifying biologically significant features, including differential gene expression and splicing variations, particularly in cancer-related studies.

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
Please refer to the smart_seq_pre_process.sh file under the tutorial.
For processed dataset, please download from here:
[Prepare datasets to run DeepExonas.] (https://mcgill-my.sharepoint.com/:f:/g/personal/kailu_song_mail_mcgill_ca/EvZtHeW7qjJJs_RHc2-327ABeLXafa-ruvfk9Vs134crig?e=jNygC6)

### Training on an example dataset
```
python ./DeepExonas/model/smart_seq/run_hyper_seed1.py
```