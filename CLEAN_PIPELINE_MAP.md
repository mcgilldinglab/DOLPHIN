# Clean Pipeline Map

This file is the current single-entry map for the cleaned DOLPHIN preprocessing workflow.

It records:

- where the cleaned code now lives
- where the cleaned documentation lives
- what each step is called
- what input each step expects from the user
- what each step writes out
- what has already been validated in this workspace
- what is still blocked by missing external dependencies

## 1. Clean Code Location

The cleaned preprocessing code is now centralized under:

- [DOLPHIN/preprocess](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess)

The main files are:

- [generate_exon_gtf.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/generate_exon_gtf.py)
- [generate_adj_index.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/generate_adj_index.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)
- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [__init__.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/__init__.py)

The cleaned execution entrypoints are now:

- `dolphin-preprocess prepare-reference`
- `dolphin-preprocess plan-preprocess`
- `dolphin-preprocess run-preprocess`
- `dolphin-preprocess recommend-parallelism`
- `dolphin-preprocess cb-to-rg`
- `dolphin-preprocess split-featurecounts`

## 2. Clean Documentation Location

The cleaned preprocessing documentation is currently split between the API docs and the tutorial docs:

- [docs/source/API/preprocess.rst](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/docs/source/API/preprocess.rst)
- [docs/source/tutorials/step1_preprocess_strategy.md](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/docs/source/tutorials/step1_preprocess_strategy.md)
- [docs/source/tutorials/step1_preprocess_cli.md](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/docs/source/tutorials/step1_preprocess_cli.md)
- [docs/source/tutorials/step1_1_preprocess_full_length.md](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/docs/source/tutorials/step1_1_preprocess_full_length.md)
- [docs/source/tutorials/step1_2_preprocess_10X.md](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/docs/source/tutorials/step1_2_preprocess_10X.md)

The tutorial index is:

- [docs/source/tutorials/index.rst](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/docs/source/tutorials/index.rst)

## 3. Current Reference Output Location

The validated step0 reference output that already exists in this workspace is:

- [reference_inputs/grch38_release107/dolphin_exon_gtf](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf)

Important files there:

- [reference_manifest.json](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/reference_manifest.json)
- [step_timing.tsv](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/step_timing.tsv)
- [dolphin_exon_gtf.gtf](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/dolphin_exon_gtf.gtf)
- [dolphin.exon.pkl](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/dolphin.exon.pkl)
- [dolphin_adj_index.csv](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/dolphin_adj_index.csv)
- [dolphin_adj_metadata_table.csv](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/dolphin_adj_metadata_table.csv)

## 4. Current Preserved Example Artifact

The preserved example full-length preprocess plan generated from the current example raw FASTQ directory is:

- [example_runs/full_length_rawfastq_only_plan](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/example_runs/full_length_rawfastq_only_plan)

Important files there:

- [preprocess_plan.json](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/example_runs/full_length_rawfastq_only_plan/preprocess_plan.json)
- [run_preprocess.sh](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/example_runs/full_length_rawfastq_only_plan/run_preprocess.sh)

## 5. Current Clean Execution Runs

The active clean execution roots in this workspace are:

- [execution_runs/step0_reference_grch38_release107_star](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/step0_reference_grch38_release107_star)
- [execution_runs/full_length_flashseq_clean_v1](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1)
- [execution_runs/tenx_srr8513796_clean_v1](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v1)

Important files already created there:

- [execution_runs/step0_reference_grch38_release107_star/run.log](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/step0_reference_grch38_release107_star/run.log)
- [execution_runs/full_length_flashseq_clean_v1/preprocess_plan.json](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/preprocess_plan.json)
- [execution_runs/full_length_flashseq_clean_v1/run_preprocess.sh](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1/run_preprocess.sh)
- [execution_runs/tenx_srr8513796_clean_v1/preprocess_plan.json](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v1/preprocess_plan.json)
- [execution_runs/tenx_srr8513796_clean_v1/run_preprocess.sh](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v1/run_preprocess.sh)

## 6. Step Inventory

### Step 0: Prepare Reference Bundle

Command:

```bash
dolphin-preprocess prepare-reference
```

Code:

- [generate_exon_gtf.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/generate_exon_gtf.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Required user input:

- `--input-gtf`

Optional user input:

- `--genome-fasta`
- `--star-binary`
- `--star-threads`
- `--star-index-dir`
- `--genome-sa-sparse-d`
- `--sjdb-overhang`

Writes:

- merged DOLPHIN exon GTF
- exon pickle
- adjacency tables
- reference manifest
- timing log
- optional STAR index if FASTA and STAR are available

Current status:

- Implemented
- Real GTF run completed
- Output preserved in [reference_inputs/grch38_release107/dolphin_exon_gtf](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf)

Current blocker:

- none for the cleaned workspace inputs
- the current in-progress execution is building the shared STAR index in [execution_runs/step0_reference_grch38_release107_star](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/step0_reference_grch38_release107_star)

### Step 1A: Plan Full-Length Preprocessing

Command:

```bash
dolphin-preprocess plan-preprocess --mode full_length
```

Code:

- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Minimum required user input:

- `--reference-manifest`
- `--raw-fastq-dir`
- `--output-dir`

Optional user input:

- `--metadata-path`
- `--read-layout`
- `--fastq-suffix`
- `--trim`
- `--java-binary`
- `--trimmomatic-jar`
- `--adapter-fasta`
- `--star-index-dir`
- `--star-binary`
- `--featurecounts-binary`
- `--star-threads`
- `--featurecounts-threads`

Additional planning artifact:

- full-length plans now also store system-aware parallel cell recommendations for `conservative`, `balanced`, and `aggressive` execution profiles

If no metadata is provided:

- sample IDs are inferred directly from FASTQ filenames

Writes:

- `preprocess_plan.json`
- `run_preprocess.sh`

Expected execution outputs when the plan is actually run:

- `01_trim/`
- `03_exon_star/`
- `04_exon_gene_cnt/`
- `04_exon_gene_cnt_grouped/`
- `05_exon_junct_cnt/`
- `05_exon_junct_cnt_grouped/`
- `logs/`

Current status:

- Implemented
- Plan generation validated with the current `flashseq/raw_fastq` directory
- Preserved example plan stored in [example_runs/full_length_rawfastq_only_plan](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/example_runs/full_length_rawfastq_only_plan)

Current blocker:

- no code blocker
- actual execution still depends on the shared STAR index under [reference_inputs/grch38_release107/dolphin_exon_gtf/star_index](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/star_index)

### Step 1A-Run: Execute Full-Length Preprocessing

Command:

```bash
dolphin-preprocess run-preprocess
```

Code:

- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Required user input:

- `--plan-manifest`

Optional user input:

- `--log-dir`
- `--timing-path`
- `--summary-path`
- `--resume`
- `--continue-on-error`
- `--max-parallel-cells`
- `--parallel-profile`

Writes:

- per-stage log files
- `execution_timing.tsv`
- `execution_summary.json`

Current workspace target:

- [execution_runs/full_length_flashseq_clean_v1](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/full_length_flashseq_clean_v1)

Current status:

- Implemented
- Backward-compatible with older full-length plans by inferring per-cell STAR and featureCounts thread settings from the saved commands when needed

### Step 1B: Plan 10x Preprocessing

Command:

```bash
dolphin-preprocess plan-preprocess --mode tenx
```

Code:

- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Minimum required user input:

- `--reference-manifest`
- `--raw-fastq-dir`
- `--output-dir`
- `--solo-whitelist`
- `--chemistry`

Optional user input:

- `--metadata-path`
- `--metadata-cb-suffix`
- `--sample-name`
- `--solo-cell-filter`
- `--samtools-binary`
- `--barcode-tag`
- `--rg-tag`
- `--append-suffix`
- `--star-index-dir`
- `--star-binary`
- `--featurecounts-binary`
- `--star-threads`
- `--featurecounts-threads`

If no metadata is provided:

- the plan assumes that `STARsolo` filtered barcodes will define the kept cells

Writes:

- `preprocess_plan.json`
- `run_preprocess.sh`

Expected execution outputs when the plan is actually run:

- `03_exon_star/<sample>/`
- `04_rg_tagged_bam/`
- `05_exon_gene_cnt_grouped/`
- `05_exon_gene_cnt/`
- `06_exon_junct_cnt_grouped/`
- `06_exon_junct_cnt/`
- `logs/`

Preferred data layout:

- grouped exon and junction matrices are the canonical count outputs
- full-length grouped exon-gene outputs are also preserved for downstream gene-level processing
- split per-cell files are kept only as backward-compatible exports for older downstream readers

Current status:

- Implemented
- Plan generation works with the local whitelist path
- A real plan is now preserved in [execution_runs/tenx_srr8513796_clean_v1](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v1)

Current blocker:

- no code blocker
- actual execution still depends on the shared STAR index under [reference_inputs/grch38_release107/dolphin_exon_gtf/star_index](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107/dolphin_exon_gtf/star_index)

### Step 1B-Run: Execute 10x Preprocessing

Command:

```bash
dolphin-preprocess run-preprocess
```

Code:

- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Required user input:

- `--plan-manifest`

Optional user input:

- `--log-dir`
- `--timing-path`
- `--summary-path`
- `--resume`
- `--continue-on-error`

Writes:

- per-stage log files
- `execution_timing.tsv`
- `execution_summary.json`

Current workspace target:

- [execution_runs/tenx_srr8513796_clean_v1](/mnt/md0/kailu/DOLPHIN_codex/execution_runs/tenx_srr8513796_clean_v1)

### Utility: CB to RG Conversion

Command:

```bash
dolphin-preprocess cb-to-rg
```

Code:

- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Required user input:

- SAM input on `stdin`

Optional user input:

- `--barcode-tag`
- `--rg-tag`
- `--append-suffix`
- `--stats-path`

Current status:

- Implemented
- Small synthetic test validated in this workspace

### Utility: Split Grouped featureCounts Output

Command:

```bash
dolphin-preprocess split-featurecounts
```

Code:

- [pipeline.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/pipeline.py)
- [cli.py](/mnt/md0/kailu/DOLPHIN_codex/DOLPHIN_optimized/DOLPHIN/preprocess/cli.py)

Required user input:

- `--mode`
- `--input-path`
- `--output-dir`

Optional user input:

- `--metadata-path`
- `--cell-id-list-path`

Current status:

- Implemented
- Small synthetic tests validated for both metadata-based and barcode-list-based splitting

## 7. Example Input Locations

Current example raw inputs inside the allowed workspace:

- full-length raw FASTQ directory: [data_example/flashseq/raw_fastq](/mnt/md0/kailu/DOLPHIN_codex/data_example/flashseq/raw_fastq)
- 10x raw FASTQ directory: [data_example/10x/raw_fastq](/mnt/md0/kailu/DOLPHIN_codex/data_example/10x/raw_fastq)
- example 10x metadata: [data_example/10x/10x_metadata.csv](/mnt/md0/kailu/DOLPHIN_codex/data_example/10x/10x_metadata.csv)
- example full-length metadata: [data_example/flashseq/metadata.csv](/mnt/md0/kailu/DOLPHIN_codex/data_example/flashseq/metadata.csv)

## 8. Current Missing External Inputs or Tools

As of the current execution state, the cleaned workspace already contains the required local tools and inputs:

- local runtime env with `STAR`, `featureCounts`, `samtools`, `pysam`, `pybedtools`
- human genome FASTA under [reference_inputs/grch38_release107](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/grch38_release107)
- 10x whitelist under [reference_inputs/10x_whitelists](/mnt/md0/kailu/DOLPHIN_codex/reference_inputs/10x_whitelists)

The remaining runtime dependency is simply the completion of the shared STAR index build.

## 9. Immediate Next Execution Order

To move from planning to actual execution, the next order should be:

1. Finish the in-progress shared STAR index build.
2. Run the preserved full-length plan with `dolphin-preprocess run-preprocess`.
3. Run the preserved 10x plan with `dolphin-preprocess run-preprocess`.
4. Compare outputs and timings against the original workflow.

## 10. Skill-Oriented Note

For later skill creation, the intended minimum preprocess contract is now:

- raw FASTQ input
- source GTF input
- optional genome FASTA input
- optional metadata

Metadata is no longer treated as a required preprocess input.
