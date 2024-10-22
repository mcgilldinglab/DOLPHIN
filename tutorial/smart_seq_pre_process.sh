#!/bin/bash
#echo "####################Processing Start#################################"

#binaries
prefetch="/mnt/data/kailu/Apps/sratoolkit/bin/prefetch"
fastdump="/mnt/data/kailu/Apps/sratoolkit/bin/fastq-dump"
trim="/mnt/data/kailu/Apps/Trimmomatic-0.39/trimmomatic-0.39.jar"
featurecount="/mnt/data/kailu/Apps/subread-2.0.3-source/bin/featureCounts"

#run the pipeline for every sample
# main_folder="/mnt/data/kailu/STAR_example/flash_seq_raw" #everything stored under main folder
input_sraList="$main_folder/SRR_Acc_List_PBMC_flash_seq_SE.txt"
output_mappingrate="$main_folder/mappingRateResults.txt"
#anno_gtf="/mnt/data/kailu/STAR_example/gencode_download/gencode.v41.primary_assembly.annotation.gtf"
anno_gtf="/mnt/data/kailu/STAR_example/ensembl_mod/Homo_sapiens.GRCh38.107.exon.gtf"

helpFunction()
{
   echo ""
   echo "Usage: $0 -M .. -Q .. -O .. --exon_star_idx .. --std_star_idx"
   echo -e "\t-M sample metaData and first column is sample names"
   echo -e "\t-Q fastq path"
   echo -e "\t-O output path"
   echo -e "\t--exon_star_idx modified gtf file star genome index"
   echo -e "\t--std_star_idx standard gtf file star genome index"
   exit 1 # Exit script after printing help
}

while getopts "M:Q:O:exon_star_idx:std_star_idx" opt
do
   case "$opt" in
      M ) metaData="$OPTARG" ;;
      Q ) fastq="$OPTARG" ;;
      O ) output="$OPTARG" ;;
	  exon_star_idx ) exon_idx="$OPTARG" ;;
	  std_star_idx ) std_idx="$OPTARG" ;;
      ? ) helpFunction ;; # Print helpFunction in case parameter is non-existent
   esac
done

# Print helpFunction in case parameters are empty
if [ -z "$metaData" ] || [ -z "$fastq" ] || [ -z "$output" ] || [ -z "$exon_idx" ] || [ -z "$std_idx" ]
then
   echo "Some or all of the parameters are empty";
   helpFunction
fi

# Begin script in case all parameters are correct
# echo "$metaData"
# echo "$fastq"
# echo "$output"

### create required folders
mkdir $output/01_std_star
mkdir $output/02_std_gene_cnt
mkdir $output/03_exon_star
mkdir $output/04_exon_cnt

# echo "$metaData"
# echo "$fastq"
# echo "$output"
echo "$exon_idx"
echo "$std_idx"

for ID_SAMPLE in $(cut -f1 $metaData)
do  
	# cd $main_folder
	echo "$ID_SAMPLE"
	
	#download sra file
	mkdir $line
	$prefetch $line
	echo "$line downloaded"

	convert to fastq files
	cd $line
	$fastdump --gzip --split-3 $line
	echo "$line converted to fastq"

	## trim
	java -jar $trim SE $main_folder/$line/$line.fastq.gz $main_folder/$line/$line.trim.fastq.gz ILLUMINACLIP:/mnt/data/kailu/Apps/Trimmomatic-0.39/adapters/TruSeq3-SE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 	
	echo "$line trimmed Success"

	### alignment
	STAR --runThreadN 4 \
		 --genomeDir /mnt/data/kailu/STAR_example/ensembl_mod_indx/ \
		 --readFilesIn $fastq_path/$ID_SAMPLE.trim.fastq.gz \
		 --readFilesCommand gunzip -c \
		 --outSAMtype BAM SortedByCoordinate \
		 --outFileNamePrefix $main_folder/$line/${line}.

	# echo "$line alignment Success"

	#store the mapping rate
	echo "$line" >> $output_mappingrate
	awk '{if(NR == 10){print $6}}' $main_folder/$line/${line}.Log.final.out >> $output_mappingrate
	awk '{if(NR == 25){print $9}}' $main_folder/$line/${line}.Log.final.out >> $output_mappingrate

	#count the exon and junction
	$featurecount -t exon -f -O -J -M \
	              -a $anno_gtf \
	              -o $main_folder/$line/${line}.count.txt \
	              $main_folder/$line/${line}.Aligned.sortedByCoord.out.bam

	script $line.feature_count_log.txt #save the log of featurecounts to txt file
   $featurecount -t exon -f -O -J -M -R SAM \
	             -a $anno_gtf \
	             -o $main_folder/$line/${line}.count.txt \
	             $main_folder/$line/${line}.Aligned.sortedByCoord.out.bam

	#feature count from paper
	$featurecount -T 1 -t exon -g gene_name --fracOverlap 0.25 \
				  -a $anno_gtf \
				  -o $main_folder/$line/$line_count.txt \
				  $main_folder/$line/$line.Aligned.sortedByCoord.out.bam

 	#featurecounts - exon gene
    /mnt/data/kailu/Apps/subread-2.0.3-source/bin/featureCounts -t exon -O -J -M \
        -a /mnt/data/kailu/STAR_example/ensembl_mod/Homo_sapiens.GRCh38.107.exon.gtf \
        -o ./12_exon_gene_count/$ID_SAMPLE.gene.count.txt \
        ./11_single_exon_bam/filtered.TAG_CB:Z_$ID_SAMPLE-1.bam

    #featurecounts - junction count
    /mnt/data/kailu/Apps/subread-2.0.3-source/bin/featureCounts -t exon -f -O -J -M \
        -a /mnt/data/kailu/STAR_example/ensembl_mod/Homo_sapiens.GRCh38.107.exon.gtf \
        -o ./13_exon_junction_count/$ID_SAMPLE.count.txt \
        ./11_single_exon_bam/filtered.TAG_CB:Z_$ID_SAMPLE-1.bam

	# move files
    mv $main_folder/$line/${line}.count.txt $main_folder/exon_count_table
    mv $main_folder/$line/${line}.count.txt.jcounts $main_folder/jun_count_table
    mv $main_folder/$line/${line}.count.txt.summary $main_folder/count_summary_file
    mv $main_folder/$line/${line}.Log.final.out $main_folder/alignment_log_file

    #delete folders
    rm -r $main_folder/$line

done
