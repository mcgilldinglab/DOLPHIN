import os
from typing import Iterable, Optional

import pandas as pd

JUNCTION_GENEID_COLUMN = "__dolphin_geneid"
JUNCTION_CHR_COLUMN = "__dolphin_chr"
JUNCTION_START_COLUMN = "__dolphin_start"
JUNCTION_END_COLUMN = "__dolphin_end"
JUNCTION_VALID_COLUMN = "__dolphin_valid"


GROUPED_FEATURECOUNTS_CONFIG = {
    "gene": {
        "comment": "#",
        "count_column_start": 6,
        "annotation_columns": (
            "Geneid",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
        ),
    },
    "exon": {
        "comment": "#",
        "count_column_start": 6,
        "annotation_columns": (
            "Geneid",
            "Chr",
            "Start",
            "End",
            "Strand",
            "Length",
        ),
    },
    "junction": {
        "comment": None,
        "count_column_start": 14,
        "annotation_columns": (
            "Gene_SP1",
            "Gene_SP2",
            "Transcript",
            "Status",
            "Donor",
            "Acceptor",
            "Chr_SP1",
            "Location_SP1",
            "Strand_SP1",
            "Chr_SP2",
            "Location_SP2",
            "Strand_SP2",
            "NearestExonBoundary_SP1",
            "NearestExonBoundary_SP2",
        ),
    },
}


def derive_grouped_cell_id(column_name: str) -> str:
    basename = os.path.basename(str(column_name))
    if ":" in basename:
        return basename.rsplit(":", 1)[-1]
    if basename.endswith(".Aligned.sortedByCoord.out.bam"):
        return basename.split(".Aligned.sortedByCoord.out.bam")[0]
    if basename.endswith(".bam"):
        return basename[:-4]
    return basename.split(".")[0]


def _normalize_count_columns(
    dataframe: pd.DataFrame,
    count_column_start: int,
) -> pd.DataFrame:
    annotation_columns = list(dataframe.columns[:count_column_start])
    original_count_columns = list(dataframe.columns[count_column_start:])
    renamed_count_columns = [
        derive_grouped_cell_id(column_name)
        for column_name in original_count_columns
    ]
    if len(set(renamed_count_columns)) != len(renamed_count_columns):
        duplicate_columns = sorted(
            {
                column_name
                for column_name in renamed_count_columns
                if renamed_count_columns.count(column_name) > 1
            }
        )
        raise ValueError(
            "Duplicate grouped featureCounts cell IDs after normalizing column "
            f"names: {duplicate_columns}"
        )

    dataframe.columns = annotation_columns + renamed_count_columns
    for column_name in renamed_count_columns:
        dataframe[column_name] = pd.to_numeric(
            dataframe[column_name],
            errors="coerce",
        ).fillna(0)
    return dataframe


def _build_grouped_usecols(mode: str, requested_cell_ids) -> Optional[callable]:
    if requested_cell_ids is None:
        return None

    requested_cell_ids = set(requested_cell_ids)
    annotation_columns = set(GROUPED_FEATURECOUNTS_CONFIG[mode]["annotation_columns"])

    def _usecol(column_name: str) -> bool:
        if column_name in annotation_columns:
            return True
        return derive_grouped_cell_id(column_name) in requested_cell_ids

    return _usecol


def load_grouped_featurecounts_table(
    input_path: str,
    mode: str,
    requested_cell_ids=None,
) -> pd.DataFrame:
    if mode not in GROUPED_FEATURECOUNTS_CONFIG:
        raise ValueError(
            f"Unsupported grouped featureCounts mode: {mode}. "
            f"Supported modes: {sorted(GROUPED_FEATURECOUNTS_CONFIG)}"
        )

    config = GROUPED_FEATURECOUNTS_CONFIG[mode]
    input_path = os.path.abspath(input_path)
    dataframe = pd.read_csv(
        input_path,
        sep="\t",
        comment=config["comment"],
        low_memory=False,
        usecols=_build_grouped_usecols(mode, requested_cell_ids),
    )
    dataframe = _normalize_count_columns(
        dataframe=dataframe,
        count_column_start=config["count_column_start"],
    )
    if requested_cell_ids is not None:
        count_columns = list(dataframe.columns[config["count_column_start"] :])
        if count_columns:
            active_row_mask = (
                dataframe[count_columns]
                .to_numpy(copy=False)
                .any(axis=1)
            )
            dataframe = dataframe.loc[active_row_mask].copy()
    if mode == "junction":
        dataframe = _prepare_grouped_junction_table(dataframe)
    return dataframe


def _split_gene_tokens(value: object) -> Iterable[str]:
    if value is None or pd.isna(value):
        return []
    tokens = []
    for token in str(value).split(","):
        token = token.strip()
        if not token or token == "NA":
            continue
        tokens.append(token)
    return tokens


def infer_primary_gene(gene_sp1: object, gene_sp2: object) -> Optional[str]:
    tokens1 = set(_split_gene_tokens(gene_sp1))
    tokens2 = set(_split_gene_tokens(gene_sp2))
    if not tokens1 and not tokens2:
        return None
    if not tokens1:
        return next(iter(tokens2)) if len(tokens2) == 1 else None
    if not tokens2:
        return next(iter(tokens1)) if len(tokens1) == 1 else None

    overlap = tokens1 & tokens2
    if len(overlap) == 1:
        return next(iter(overlap))

    if len(tokens1) == 1 and len(tokens2) == 1 and tokens1 == tokens2:
        return next(iter(tokens1))
    return None


def _prepare_grouped_junction_table(dataframe: pd.DataFrame) -> pd.DataFrame:
    dataframe = dataframe.copy()
    if JUNCTION_VALID_COLUMN in dataframe.columns:
        return dataframe

    if "PrimaryGene" in dataframe.columns:
        geneid_series = dataframe["PrimaryGene"].copy()
        chr_series = dataframe["Site1_chr"].copy()
        start_series = pd.to_numeric(
            dataframe["Site1_location"],
            errors="coerce",
        )
        end_series = pd.to_numeric(
            dataframe["Site2_location"],
            errors="coerce",
        )
    else:
        geneid_series = pd.Series(
            [
                infer_primary_gene(gene_sp1, gene_sp2)
                for gene_sp1, gene_sp2 in zip(
                    dataframe["Gene_SP1"].tolist(),
                    dataframe["Gene_SP2"].tolist(),
                )
            ],
            index=dataframe.index,
            dtype="object",
        )
        chr_series = dataframe["Chr_SP1"].copy()
        start_series = pd.to_numeric(
            dataframe["Location_SP1"],
            errors="coerce",
        )
        end_series = pd.to_numeric(
            dataframe["Location_SP2"],
            errors="coerce",
        )

    valid_mask = (
        geneid_series.notna()
        & start_series.notna()
        & end_series.notna()
    )

    dataframe[JUNCTION_GENEID_COLUMN] = geneid_series
    dataframe[JUNCTION_CHR_COLUMN] = chr_series
    dataframe[JUNCTION_START_COLUMN] = start_series.fillna(0).astype(int)
    dataframe[JUNCTION_END_COLUMN] = end_series.fillna(0).astype(int)
    dataframe[JUNCTION_VALID_COLUMN] = valid_mask.to_numpy(dtype=bool, copy=False)
    return dataframe


def normalize_junction_count_dataframe(dataframe: pd.DataFrame) -> pd.DataFrame:
    dataframe = dataframe.copy()

    if "Count" not in dataframe.columns:
        dataframe = dataframe.rename(columns={dataframe.columns[-1]: "Count"})
    dataframe["Count"] = pd.to_numeric(dataframe["Count"], errors="coerce").fillna(0)

    if JUNCTION_VALID_COLUMN in dataframe.columns:
        dataframe = dataframe[dataframe[JUNCTION_VALID_COLUMN]]
        dataframe = dataframe[
            [
                JUNCTION_GENEID_COLUMN,
                JUNCTION_CHR_COLUMN,
                JUNCTION_START_COLUMN,
                JUNCTION_END_COLUMN,
                "Count",
            ]
        ]
        dataframe = dataframe.rename(
            columns={
                JUNCTION_GENEID_COLUMN: "Geneid",
                JUNCTION_CHR_COLUMN: "Chr",
                JUNCTION_START_COLUMN: "Start",
                JUNCTION_END_COLUMN: "End",
            }
        )
    elif "PrimaryGene" in dataframe.columns:
        dataframe = dataframe.dropna(subset=["PrimaryGene"])
        dataframe = dataframe[["PrimaryGene", "Site1_chr", "Site1_location", "Site2_location", "Count"]]
        dataframe = dataframe.rename(
            columns={
                "PrimaryGene": "Geneid",
                "Site1_chr": "Chr",
                "Site1_location": "Start",
                "Site2_location": "End",
            }
        )
    else:
        dataframe["Geneid"] = dataframe.apply(
            lambda row: infer_primary_gene(row.get("Gene_SP1"), row.get("Gene_SP2")),
            axis=1,
        )
        dataframe = dataframe.dropna(subset=["Geneid"])
        dataframe = dataframe.rename(
            columns={
                "Chr_SP1": "Chr",
                "Location_SP1": "Start",
                "Location_SP2": "End",
            }
        )
        dataframe = dataframe[["Geneid", "Chr", "Start", "End", "Count"]]

    dataframe["Start"] = pd.to_numeric(dataframe["Start"], errors="coerce")
    dataframe["End"] = pd.to_numeric(dataframe["End"], errors="coerce")
    dataframe = dataframe.dropna(subset=["Start", "End"])
    dataframe["Start"] = dataframe["Start"].astype(int)
    dataframe["End"] = dataframe["End"].astype(int)
    dataframe = dataframe[dataframe["Count"] > 0]
    dataframe["Type"] = "Junction"
    return dataframe


def extract_grouped_cell_exon_counts(
    grouped_exon_counts: pd.DataFrame,
    cell_id: str,
) -> pd.DataFrame:
    if cell_id not in grouped_exon_counts.columns:
        raise KeyError(f"Cell '{cell_id}' not found in grouped exon counts")

    count_values = grouped_exon_counts[cell_id].to_numpy(copy=False)
    positive_mask = count_values > 0
    dataframe = grouped_exon_counts.loc[
        positive_mask,
        ["Geneid", "Chr", "Start", "End"],
    ].copy()
    dataframe["Count"] = count_values[positive_mask]
    dataframe["Type"] = "Exon"
    return dataframe


def extract_grouped_cell_junction_counts(
    grouped_junction_counts: pd.DataFrame,
    cell_id: str,
) -> pd.DataFrame:
    if cell_id not in grouped_junction_counts.columns:
        raise KeyError(f"Cell '{cell_id}' not found in grouped junction counts")

    if JUNCTION_VALID_COLUMN in grouped_junction_counts.columns:
        count_values = grouped_junction_counts[cell_id].to_numpy(copy=False)
        valid_mask = grouped_junction_counts[JUNCTION_VALID_COLUMN].to_numpy(
            dtype=bool,
            copy=False,
        )
        active_mask = valid_mask & (count_values > 0)
        dataframe = grouped_junction_counts.loc[
            active_mask,
            [
                JUNCTION_GENEID_COLUMN,
                JUNCTION_CHR_COLUMN,
                JUNCTION_START_COLUMN,
                JUNCTION_END_COLUMN,
            ],
        ].copy()
        dataframe = dataframe.rename(
            columns={
                JUNCTION_GENEID_COLUMN: "Geneid",
                JUNCTION_CHR_COLUMN: "Chr",
                JUNCTION_START_COLUMN: "Start",
                JUNCTION_END_COLUMN: "End",
            }
        )
        dataframe["Count"] = count_values[active_mask]
        dataframe["Type"] = "Junction"
        return dataframe

    selected_columns = list(grouped_junction_counts.columns[:14]) + [cell_id]
    dataframe = grouped_junction_counts[selected_columns].copy()
    dataframe = dataframe.rename(columns={cell_id: "Count"})
    return normalize_junction_count_dataframe(dataframe)
