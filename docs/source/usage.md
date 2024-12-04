# Usage

(installation)=
## Installation

Install from Github

```bash
  $ git clone https://github.com/mcgilldinglab/DOLPHIN.git
  $ cd DOLPHIN
```

Create a conda environment using provided yaml file
```bash
  $ conda env create -f environment.yaml
  $ conda activate DOLPHIN
```

### Prerequisites
Align raw RNA-seq dataset and count the exon reads and junction reads
-   STAR >= 2.7.3a
-   cellranger >= 7.0.1
-   featurecounts

## Example data

- For quick test on the python packages succeffully installed, please download test dataset: [One Drive.](https://outlook.office.com/host/377c982d-9686-450e-9a7c-22aeaf1bc162/7211f19f-262a-42eb-a02e-289956491741).
- For run the tutorial and reproduce the results in the manuscripts, please download the processed dataset here: [One Drive.](https://outlook.office.com/host/377c982d-9686-450e-9a7c-22aeaf1bc162/7211f19f-262a-42eb-a02e-289956491741)

**Required files**

For human alignment to the exon level, please download the pre-processed gtf file here:[One Drive.](https://outlook.office.com/host/377c982d-9686-450e-9a7c-22aeaf1bc162/7211f19f-262a-42eb-a02e-289956491741).


