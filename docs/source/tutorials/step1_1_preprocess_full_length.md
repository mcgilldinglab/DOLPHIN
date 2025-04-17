# Preprocessing Full-Length Single-Cell RNA-Seq for Exon and Junction Read Counts

Here is the brief pipeline for full-length and 10x single-cell RNA-seq shown:

<img title="preprocess pipeline" alt="Alt text" src="preprocess_pipeline.png" />
<img src="../_static/preprocess_pipeline.png" alt="preprocess" width="600"/>

### Step 1: Download Required Tools

Before starting the alignment process, make sure to download and install the following tools:

[STAR](https://github.com/alexdobin/STAR) >=2.7.3a

[featurecounts](https://sourceforge.net/projects/subread/files/subread-2.0.8/) >=2.0.3

optional:
[Trimmomatic](http://www.usadellab.org/cms/index.php?page=trimmomatic) >=0.39

### Step 2: Create a Reference Genome

Run the following command to generate a reference genome for alignment using STAR. 
- `ensembl_mod_indx` is the directory where the reference genome index will be stored.
- `Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa` can be downloaded [here](https://ftp.ensembl.org/pub/release-113/fasta/homo_sapiens/dna/).
- `Homo_sapiens.GRCh38.107.exon.gtf` is generated using the [file](./step0_generate_exon_gtf.ipynb).


```bash
STAR --runMode genomeGenerate \
    --genomeDir /mnt/data/kailu/STAR_example/ensembl_mod_indx/ \
    --genomeFastaFiles Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa \
    --sjdbGTFfile Homo_sapiens.GRCh38.107.exon.gtf \
    --runThreadN 16
```

### Step 3: Download the Raw RNA-Seq Files

Download the raw RNA-seq files from the provided sources. For the links to the raw data, please refer to the original [study](https://www.nature.com/articles/s41587-022-01312-3#data-availability).

For full-length single-cell RNA-seq, each cell is stored in a separate FASTQ file. In the following steps, we will process one cell at a time. For example, the codes below processe cell ${ID_SAMPLE}

### Step4: Trim 
```bash
# location of the timmomatic tools
trim = "/mnt/data/kailu/Apps/Trimmomatic-0.39/trimmomatic-0.39.jar"

java -jar $trim SE ${ID_SAMPLE}.fastq.gz ${ID_SAMPLE}.trim.fastq.gz ILLUMINACLIP:/mnt/data/kailu/Apps/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 	
```

### Step5: STAR Alignment - align to modifed exon gtf file and standard reference genome
```bash
## `ID_SAMPLE` is the Cell Barcode Name
mkdir ./03_exon_star/${ID_SAMPLE}
STAR --runThreadN 4 \
    --genomeDir /mnt/data/kailu/STAR_example/ensembl_mod_indx/ \
    --readFilesIn ${ID_SAMPLE}.trim.fastq.gz  \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ./03_exon_star/${ID_SAMPLE}/${ID_SAMPLE}.

mkdir ./02_exon_std/${ID_SAMPLE}
STAR --runThreadN 4 \
    --genomeDir /mnt/data/kailu/STAR_example/ensembl_indx/ \
    --readFilesIn ${ID_SAMPLE}.trim.fastq.gz  \
    --readFilesCommand gunzip -c \
    --outSAMtype BAM SortedByCoordinate \
    --outFileNamePrefix ./02_exon_std/${ID_SAMPLE}/${ID_SAMPLE}.std.
```

### Step 6: Count Exon Reads and Junction Reads

Get exon gene count using the modified exon GTF file. This will generate the gene count (`${ID_SAMPLE}.exongene.count.txt`), which will be used later for HVG identification.

```bash
mkdir ./04_exon_gene_cnt
featureCounts -t exon -O -M \
    -a Homo_sapiens.GRCh38.107.exon.gtf \
    -o ./04_exon_gene_cnt/${ID_SAMPLE}.exongene.count.txt \
    ./03_exon_star/${ID_SAMPLE}/${ID_SAMPLE}.Aligned.sortedByCoord.out.bam
```

Run the following command to get the exon and junction counts. This step will generate the following files:
- `${ID_SAMPLE}.exon.count.txt`: Exon read counts.
- `${ID_SAMPLE}.exon.count.txt.jcounts`: Junction read counts.

```bash
mkdir ./05_exon_junct_cnt
featureCounts -t exon -f -O -J -M \
    -a Homo_sapiens.GRCh38.107.exon.gtf \
    -o ./05_exon_junct_cnt/${ID_SAMPLE}.exon.count.txt \
    ./03_exon_star/${ID_SAMPLE}/${ID_SAMPLE}.Aligned.sortedByCoord.out.bam
```