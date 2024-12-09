# Detect Alternative Splicing Events with Outrigger

```bash

#!/bin/bash

for ID_SAMPLE in $(tail -n +2 ../your_metaData.csv | cut -d $'\t' -f1)
do  
	echo "$ID_SAMPLE"

    mkdir ./combine_KNN10_STAR/$ID_SAMPLE
    STAR --runThreadN 16 \
		--genomeDir /mnt/data/kailu/STAR_example/ensembl_indx/ \
		--readFilesIn ./DOLPHIN_aggregation/final_bam/$ID_SAMPLE.final.bam \
		--readFilesCommand samtools view -F 0x100 \
		--outSAMtype BAM Unsorted \
		--readFilesType SAM SE \
		--outFileNamePrefix ./combine_KNN10_STAR/$ID_SAMPLE/${ID_SAMPLE}.exon.

	scp -r ./combine_KNN10_STAR/$ID_SAMPLE/${ID_SAMPLE}.exon.SJ.out.tab ./combine_KNN10_SJ/

	rm -r ./combine_KNN10_STAR/$ID_SAMPLE
done

cd ./combine_KNN10_SJ/
outrigger index --sj-out-tab *SJ.out.tab --gtf Homo_sapiens.GRCh38.107.gtf
### modify the bed file to add "chr" in the beginning of each bed file
for file in ./outrigger_output/index/mxe/*.bed
do
    echo $file
    sed -i -e 's/^/chr/' $file
done
for file in ./outrigger_output/index/se/*.bed
do
    echo $file
    sed -i -e 's/^/chr/' $file
done
## using modified fasta file to validate the event
outrigger validate --genome hg38 --fasta ./Homo_sapiens.GRCh38.dna_sm.primary_assembly.addchr.fa
outrigger psi
```