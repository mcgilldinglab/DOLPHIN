### step1: ### This function processes exon count and junction raw count data for each cell and converts them into
### flattened feature and adjacency vectors.
from DOLPHIN.graph_generation.preprocess_raw_reads import run_parallel_gene_processing

run_parallel_gene_processing(
    metadata_path="/fsla_meta.csv",
    gtf_path="./dolphin_exon_gtf/dolphin.exon.pkl",
    adj_index_path="./dolphin_exon_gtf/dolphin_adj_index.csv",
    main_folder="./",
    n_processes=8
)

### step2: This function combines feature vectors and constructs the final feature matrix.
from DOLPHIN.graph_generation.process_feature_matrix import run_feature_combination

run_feature_combination(
    metadata_path="./fsla_meta.csv",
    graph_directory="./06_graph_mtx",
    gene_annotation="./dolphin_exon_gtf/dolphin_gene_meta.csv",
    gtf_pkl_path="./dolphin_exon_gtf/dolphin.exon.pkl",
    out_directory=".",
    out_name= "fsla",
    clean_temp=False
)

### step3: This function combines adjacency vectors and constructs the Adjacency matrix.
from DOLPHIN.graph_generation.process_adjacency_matrix import run_adjacency_combination

run_adjacency_combination(
    metadata_path="./fsla_meta.csv",
    graph_directory="./06_graph_mtx",
    adj_meta_file="./dolphin_exon_gtf/dolphin_adj_metadata_table.csv",
    out_directory=".",
    out_name= "fsla",
    clean_temp=False,
    adj_run_num=50,
    parallel=True 
)

## step4: This function convert adjacency matrix to compressed Adjacency matrix.
from DOLPHIN.graph_generation.process_adjacency_matrix_compress import run_adjacency_compression

run_adjacency_compression(
    metadata_path="./fsla_meta.csv",
    out_name= "fsla",
    out_directory=".",
    num_processes=50
)

## step5: This function combine comprssed adjacency matrix.
from DOLPHIN.graph_generation.process_adjacency_matrix_compress_combine import run_adjacency_compress_combination

run_adjacency_compress_combination(
    metadata_path= "./fsla_meta.csv",
    out_name= "fsla",
    out_directory="./",
    adj_run_num=50,
    clean_temp=False,
    parallel= True,
)

### step6: This function clean the final adjacency matrix.
from DOLPHIN.graph_generation.process_adjacency_matrix_final import run_adjacency_matrix_final

run_adjacency_matrix_final(
    out_name="fsla",
    out_directory="./"
)

# step7: This step is primarily used to select highly variable genes (HVGs) based on the gene count table.
# The final feature matrix and adjacency matrix will retain only these HVGs as input for the DOLPHIN model.
# This step is optional—it provides one way to select HVGs, but you are free to use your own selection method.

from DOLPHIN.graph_generation.process_raw_gene import run_raw_gene

run_raw_gene(
    metadata_path= "./fsla_meta.csv",
    featurecount_path= "./04_exon_gene_cnt",
    gtf_path="./dolphin_exon_gtf/dolphin.exon.gtf",
    out_name="fsla",
    out_directory="./")

# step8: This step is generate hvg selected feature matrix

from DOLPHIN.graph_generation.process_feature_hvg import run_feature_hvg

run_feature_hvg(
    out_name="fsla",
    out_directory="./")

# step9: This step generates the HVG-selected adjacency matrix.

from DOLPHIN.graph_generation.process_adjacency_hvg import run_adjacency_hvg

run_adjacency_hvg(
    out_name="fsla",
    out_directory="./")

## step10: This final step generates the model input
from DOLPHIN.graph_generation.process_graph_final import run_model_input

run_model_input(
    metadata_path= "./fsla_meta.csv",
    out_name="fsla",
    out_directory="./",
    celltypename="celltype1")