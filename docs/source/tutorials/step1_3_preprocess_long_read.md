# Preprocessing Long-read Single-Cell RNA-Seq for Exon and Junction Read Counts

Here is the brief pipeline for long-read single-cell RNA-seq shown:

![preprocess pipeline](../_static/preprocess_pipeline_long.png)

## Step 1: Download Required Tools

Before starting the alignment process, make sure to download and install the following tools:

[STAR](https://github.com/alexdobin/STAR) >=2.7.3a

[featurecounts](https://sourceforge.net/projects/subread/files/subread-2.0.8/) >=2.0.3

[cellranger](https://www.10xgenomics.com/support/software/cell-ranger/latest/tutorials/cr-tutorial-in) >= 7.0.1

[subset-bam](https://github.com/10XGenomics/subset-bam)

[bamtools](https://github.com/pezmaster31/bamtools)

[scNanoGPS](https://github.com/gaolabtools/scNanoGPS)

## Step 2: Follow the scNanoGPS pipeline to obtain the curated BAM file.

## Step 3: Count Exon Reads and Junction Reads

Get exon gene count using the modified exon GTF file. This will generate the gene count (`${ID_SAMPLE}.exongene.count.txt`), which will be used later for HVG identification.

```bash
mkdir ./04_exon_gene_cnt
featureCounts -t exon -O -M \
    -a ./dolphin_exon_gtf/dolphin_exon_gtf.gtf \
    -o ./04_exon_gene_cnt/${ID_SAMPLE}.exongene.count.txt \
    ./00_scNanogps/MEL2/curated.minimap2.bam/$ID_SAMPLE.curated.minimap2.bam
```

Run the following command to get the exon and junction counts. This step will generate the following files:
- `${ID_SAMPLE}.exon.count.txt`: Exon read counts.
- `${ID_SAMPLE}.exon.count.txt.jcounts`: Junction read counts.

```bash
mkdir ./05_exon_junct_cnt
featureCounts -t exon -f -O -J -M \
    -a ./dolphin_exon_gtf/dolphin_exon_gtf.gtf \
    -o ./05_exon_junct_cnt/${ID_SAMPLE}.exon.count.txt \
    ./00_scNanogps/MEL2/curated.minimap2.bam/$ID_SAMPLE.curated.minimap2.bam
```