# Alternative Splicing Events Detection

In our manuscript, the default method we used is [Expedition](https://pubmed.ncbi.nlm.nih.gov/28673540/). Since we are aggregating 
single-cell RNA raw reads (BAM) file, you can try with other alternative splicing detection tools. Please install [Outrigger](https://github.com/YeoLab/outrigger) and [Anchor](https://github.com/YeoLab/anchor) modules from Expedition follow the tutorial on their website, DOLPHIN used both these functions to do downstream analysis. Outrigger is used for Alternative splicing detection and anchor is used for splicing modality quantification. In addtion, we also recommended use [MARVEL](https://pmc.ncbi.nlm.nih.gov/articles/PMC10018366/) for Alternative splicing detection and it also support splicing modality analysis.  

## Step 1: STAR Alignment

From [step4](./step4_cell_aggregation.ipynb) we get aggregated cells, which will be used for alternative splicing analysis. We will frist align the aggregated bam file to reference genome to get the SJ.out.tab files, you can use the script below

```bash

STAR --runThreadN 16 \
	--genomeDir ./ensembl_indx/ \
	--readFilesIn ./cell_aggregation/$ID_SAMPLE.aggr.final.bam \
	--readFilesCommand samtools view -F 0x100 \
	--outSAMtype BAM Unsorted \
	--readFilesType SAM SE \
	--outFileNamePrefix ./cell_aggregation/STAR/$ID_SAMPLE/${ID_SAMPLE}.aggr.
	
scp -r ./cell_aggregation/STAR/$ID_SAMPLE/${ID_SAMPLE}.aggr.SJ.out.tab ./cell_aggregation/STAR/SJ/
```

## Step 2: Run Outrigger
You can follow the Outrigger [tutorial](https://yeolab.github.io/outrigger/subcommands/outrigger_index.html) to detect the alternative splicing events. Below is an running example, you can modify the parameters based on your data.

```bash
cd ./cell_aggregation/STAR/SJ/
outrigger index --sj-out-tab *SJ.out.tab --gtf Homo_sapiens.GRCh38.107.gtf
outrigger validate --genome hg38 --fasta ./Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa 
outrigger psi
```

After running this step, you will get the outrigger output and we will use psi folder for downstream analysis. 