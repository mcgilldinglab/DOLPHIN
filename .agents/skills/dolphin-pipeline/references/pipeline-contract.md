# DOLPHIN Pipeline Contract

Use this reference when changing or running the finalized DOLPHIN pipeline.

## Package Roots

Final code package:

- `DOLPHIN/preprocess`
- `DOLPHIN/graph_generation`
- `DOLPHIN/model`
- `DOLPHIN/alternative_splicing`
- `DOLPHIN/EDEG`

Primary docs:

- `README.md`
- `docs/PIPELINE_OUTPUTS.md`
- `docs/ALTERNATIVE_SPLICING_OUTPUTS.md`
- `docs/ORIGINAL_FUNCTION_COVERAGE.md`
- `docs/preprocessing/README.md`
- `docs/alternative_splicing/README.md`
- `docs/edeg/README.md`
- `docs/model/README.md`

## End-To-End Flow

1. Reference preparation
2. Raw preprocessing for full-length or 10x
3. Graph generation
4. Lightweight graph-store model input
5. DOLPHIN model training
6. Embedding export to `DOLPHIN_Z.h5ad`
7. Downstream alternative splicing and EDEG/JDEG analyses

## Preprocessing And Graph Outputs

Reference preparation retains:

- `reference_manifest.json`
- `dolphin.exon.gtf`
- `dolphin.exon.pkl`
- `dolphin_adj_index.csv`
- `dolphin_adj_metadata_table.csv`
- `dolphin_gene_meta.csv`
- optional `star_index/`

Raw preprocessing retains grouped outputs:

- `cell_manifest.tsv`
- grouped gene count table
- grouped exon count table
- grouped junction count table

Graph generation retains:

- `Feature_<out_name>.h5ad`
- `FeatureComp_<out_name>.h5ad`
- `Adjacency_<out_name>.h5ad`
- `AdjacencyComp_<out_name>.h5ad`
- `AdjacencyCompRe_<out_name>.h5ad`
- `ExonGene_<out_name>.h5ad`
- `ExonGene_hvg_<out_name>.h5ad`
- `FeatureCompHvg_<out_name>.h5ad`
- `AdjacencyCompReHvg_<out_name>.h5ad`
- `AdjacencyCompReHvgEdge_<out_name>.h5ad`

Model input retains:

- `model_<out_name>.graph.json`
- `model_<out_name>.edges.h5`

Model training retains:

- `DOLPHIN_Z.h5ad`

Do not treat per-cell `*_fea.csv`, per-cell `*_adj.csv`, or per-cell graph H5AD fragments as final public outputs.

## Alternative Splicing Contract

The finalized AS package is `DOLPHIN/alternative_splicing`.

Preserve four routes:

- full-length BAM aggregation
- full-length direct junction aggregation
- 10x BAM aggregation
- 10x direct junction aggregation

Shared AS logic:

- find nearest neighbors from DOLPHIN embeddings
- count reads
- aggregate evidence by BAM or junction route
- run Outrigger `index`, `validate`, and `psi`
- convert to `PSI.h5ad`, `PSI_random.h5ad`, and `PSI_DAS.h5ad`
- impute missing PSI per event within cell type, then run two-sided Wilcoxon rank-sum tests between two biological groups with BH correction

BAM route:

- Aggregates neighbor-aware BAM evidence.
- Runs STAR realignment on aggregated BAM.
- Retains aggregated BAM outputs under `cell_aggregation/`.

Direct junction route:

- Aggregates neighbor-aware junction counts directly.
- Does not create aggregated BAMs.
- Does not perform STAR realignment.
- Retains STAR-style aggregated junction summaries under `aggregated_star/`.

Full-length input assumptions:

- per-cell BAMs already exist
- per-cell STAR `SJ.out.tab` files already exist
- no barcode normalization is needed

10x input assumptions:

- BAM route consumes pooled RG-tagged BAM and prepares a per-cell BAM cache with `samtools split -d RG`
- direct junction route consumes per-cell `.jcounts`
- direct junction route backfills `annotated`, `motif`, and `strand` from the reference/GTF because `.jcounts` do not carry STAR splice semantics
- `max_overhang` remains `0` for direct-junction 10x; do not infer it from coordinates

Outrigger boundary:

- Do not change Outrigger event/index/PSI semantics.
- Compatibility changes can patch pandas, SQLite timeout, gffutils DB handling, and package import stability.
- If `outrigger index` fails, capture direct stderr/stdout and inspect partial `outrigger_output` before rerunning.

AS retained outputs:

- `metadata_prepared.tsv`
- `metadata_preparation_summary.json`
- `N_<out_name>_<k>.csv`
- `<out_name>_read_counts.csv`
- `aggregated_star/`
- `outrigger_output/`
- `alternative_splicing/<out_name>_PSI.h5ad`
- `alternative_splicing/<out_name>_PSI_random.h5ad`
- `alternative_splicing/<out_name>_PSI_DAS.h5ad`
- `alternative_splicing/<out_name>_DAS.csv`
- `<logs_root>/<out_name>.status.json`
- `<logs_root>/<out_name>.summary.json`

When comparing BAM and junction AS routes, retain:

- `bam_vs_junction_comparison.json`

## EDEG/JDEG Contract

The finalized EDEG package is `DOLPHIN/EDEG`.

Public workflow:

1. Convert `.h5ad` to Seurat `.rds`.
2. Run Seurat `FindMarkers(..., test.use = "MAST")` in one-cluster-vs-rest mode.
3. Write one `<cluster>_MAST_DE.csv` per cluster.
4. Aggregate exon-level markers to gene-level EDEG.
5. Aggregate junction-level markers to gene-level JDEG.

Do not expose pairwise MAST as the default public route.

EDEG historical input/output:

- input `FeatureComp_<out_name>.h5ad`
- converted `FeatureComp_<out_name>.rds`
- MAST output `<cluster>_MAST_DE.csv`
- final output `<cluster>_EDEG.csv`

JDEG historical input/output:

- input is the junction-analysis matrix used for Seurat/MAST
- MAST output `<cluster>_MAST_DE.csv`
- final output `<cluster>_JDEG.csv`

Interpretation:

- MAST identifies differential exon-level or junction-level features first.
- `generate_EDEG.py` combines exon evidence across exons under the same gene.
- `generate_JDEG.py` combines junction evidence across junction features under the same gene.

## GitHub-Ready Packaging

The clean package should not contain runtime result files:

- no `.h5ad`
- no `.bam` or `.bai`
- no `.pt`
- no `.rds`
- no `.h5`
- no `.so`
- no `.pyc`
- no `__pycache__/`

`setup.py` exposes:

- `dolphin-prepare-reference`
- `dolphin-preprocess`
- `dolphin-as`
- `dolphin-edeg`

`setup.py` includes package data needed for R scripts and AS runtime support.
