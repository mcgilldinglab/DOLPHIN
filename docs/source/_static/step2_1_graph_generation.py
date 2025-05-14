### step1: ### This function processes exon count and junction raw count data for each cell and converts them into
### flattened feature and adjacency vectors.
# from DOLPHIN.graph_generation.preprocess_raw_reads import run_parallel_gene_processing

# run_parallel_gene_processing(
#     metadata_path="/fsla_meta.csv",
#     gtf_path="./dolphin_exon_gtf/dolphin.exon.pkl",
#     adj_index_path="./dolphin_exon_gtf/dolphin_adj_index.csv",
#     main_folder="./",
#     n_processes=8
# )

# ### step2: This function combines feature vectors and constructs the final feature matrix.
# from DOLPHIN.graph_generation.process_feature_matrix import run_feature_combination

# run_feature_combination(
#     metadata_path="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/fsla_meta.csv",
#     graph_directory="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/06_graph_mtx",
#     gene_annotation="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/dolphin_exon_gtf/dolphin_gene_meta.csv",
#     gtf_pkl_path="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/dolphin_exon_gtf/dolphin.exon.pkl",
#     out_directory="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial",
#     out_name= "fsla",
#     clean_temp=False
# )

# ### step3: This function combines adjacency vectors and constructs the Adjacency matrix.
# from DOLPHIN.graph_generation.process_adjacency_matrix import run_adjacency_combination

# run_adjacency_combination(
#     metadata_path="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/fsla_meta.csv",
#     graph_directory="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/06_graph_mtx",
#     adj_meta_file="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/dolphin_exon_gtf/dolphin_adj_metadata_table.csv",
#     out_directory="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial",
#     out_name= "fsla",
#     clean_temp=False,
#     adj_run_num=50,
#     parallel=True 
# )

### step4: This function convert adjacency matrix to compressed Adjacency matrix.
from DOLPHIN.graph_generation.process_adjacency_matrix_compress import run_adjacency_compression

run_adjacency_compression(
    metadata_path="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial/fsla_meta.csv",
    out_name= "fsla",
    out_directory="/mnt/md1/kailu/DOLPHIN_run_input_output/DOLPHIN_tutorial",
    num_processes=50
)