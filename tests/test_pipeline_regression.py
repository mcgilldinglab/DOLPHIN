import ast
import tempfile
import unittest
from pathlib import Path

import anndata
import numpy as np
import pandas as pd

from DOLPHIN.alternative_splicing.generate_differential_as import (
    calculate_differential_as,
    impute_psi_by_cell_type,
    run_differential_as,
)
from DOLPHIN.graph_generation.func_preprocess_raw_reads import (
    build_raw_adjacency_flat,
)
from DOLPHIN.graph_generation._anndata_compat import enable_nullable_string_writes
from DOLPHIN.graph_generation.process_adjacency_hvg import _normalize_within_gene


class GraphPipelineRegressionTests(unittest.TestCase):
    def test_adjacency_normalization_preserves_flattened_var_order(self):
        adata = anndata.AnnData(
            X=np.array([[1.0, 2.0, 7.0, 3.0, 1.0]], dtype=np.float32),
            var=pd.DataFrame(
                {"gene_id": ["geneA", "geneA", "geneA", "geneB", "geneB"]},
                index=["geneA-1", "geneA-2", "geneA-10", "geneB-1", "geneB-2"],
            ),
        )

        normalized = _normalize_within_gene(adata).toarray()

        np.testing.assert_allclose(
            normalized,
            np.array([[0.1, 0.2, 0.7, 0.75, 0.25]], dtype=np.float32),
        )

    def test_raw_adjacency_contains_only_observed_supported_junctions(self):
        adjacency = build_raw_adjacency_flat(
            exon_start=[100, 300, 500],
            exon_end=[199, 399, 599],
            exon_count=[5, 0, 4],
            junction_start=[199, 199, 199],
            junction_end=[300, 500, 500],
            junction_weight=[11, 7, 2],
        ).reshape(3, 3)

        expected = np.zeros((3, 3), dtype=float)
        expected[0, 2] = 9
        np.testing.assert_array_equal(adjacency, expected)

    def test_no_observed_junctions_does_not_create_dummy_edges(self):
        adjacency = build_raw_adjacency_flat(
            exon_start=[100, 300],
            exon_end=[199, 399],
            exon_count=[5, 4],
            junction_start=[],
            junction_end=[],
            junction_weight=[],
        )

        np.testing.assert_array_equal(adjacency, np.zeros(4))

    def test_gat_layers_are_configured_for_scalar_edge_features(self):
        model_path = Path(__file__).parents[1] / "DOLPHIN" / "model" / "model.py"
        tree = ast.parse(model_path.read_text())
        gat_calls = [
            node
            for node in ast.walk(tree)
            if isinstance(node, ast.Call)
            and isinstance(node.func, ast.Name)
            and node.func.id == "GATConv"
        ]

        self.assertGreater(len(gat_calls), 0)
        for call in gat_calls:
            edge_dim = next(
                (keyword.value for keyword in call.keywords if keyword.arg == "edge_dim"),
                None,
            )
            fill_value = next(
                (keyword.value for keyword in call.keywords if keyword.arg == "fill_value"),
                None,
            )
            self.assertIsInstance(edge_dim, ast.Constant)
            self.assertEqual(edge_dim.value, 1)
            self.assertIsInstance(fill_value, ast.Constant)
            self.assertEqual(fill_value.value, 0.0)


class DifferentialASPipelineRegressionTests(unittest.TestCase):
    def test_imputation_is_event_specific_within_cell_type(self):
        psi = pd.DataFrame(
            {
                "event1": [0.2, np.nan, 0.8, np.nan],
                "event2": [0.4, 0.6, np.nan, 1.0],
            },
            index=["c1", "c2", "c3", "c4"],
        )
        cell_types = pd.Series(["A", "A", "B", "B"], index=psi.index)

        imputed = impute_psi_by_cell_type(psi, cell_types)

        self.assertEqual(imputed.loc["c2", "event1"], 0.2)
        self.assertEqual(imputed.loc["c3", "event2"], 1.0)
        self.assertEqual(imputed.loc["c4", "event1"], 0.8)

    def test_differential_as_reports_wilcoxon_and_bh_values(self):
        psi = pd.DataFrame(
            {
                "event1": [0.1, 0.2, 0.8, 0.9],
                "event2": [0.5, 0.5, 0.5, 0.5],
            },
            index=["c1", "c2", "c3", "c4"],
        )
        groups = pd.Series(["control", "control", "case", "case"], index=psi.index)

        results = calculate_differential_as(
            psi,
            groups,
            group1="case",
            group2="control",
        )

        self.assertEqual(set(results["event_id"]), {"event1", "event2"})
        self.assertIn("p_value_adj_bh", results)
        event1 = results.set_index("event_id").loc["event1"]
        self.assertAlmostEqual(event1["delta_psi"], 0.7)
        self.assertTrue(np.isfinite(event1["p_value_adj_bh"]))

    def test_end_to_end_writes_imputed_matrix_and_statistics(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            source = Path(tmpdir) / "PSI.h5ad"
            adata = anndata.AnnData(
                X=np.array(
                    [
                        [0.1, 0.2],
                        [np.nan, 0.3],
                        [0.8, np.nan],
                        [0.9, 0.7],
                    ],
                    dtype=np.float32,
                ),
                obs=pd.DataFrame(
                    {
                        "celltype1": ["T", "T", "T", "T"],
                        "Condition": ["control", "control", "case", "case"],
                    },
                    index=["c1", "c2", "c3", "c4"],
                ),
                var=pd.DataFrame(index=["event1", "event2"]),
            )
            enable_nullable_string_writes()
            adata.write_h5ad(source)

            output = run_differential_as(
                outrigger_psi_data=str(source),
                out_name="test",
                cluster_name="celltype1",
                out_directory=tmpdir,
                n_cell=2,
                group_column="Condition",
                group1="case",
                group2="control",
            )

            self.assertFalse(np.isnan(output.X).any())
            self.assertTrue(
                (Path(tmpdir) / "alternative_splicing" / "test_PSI_DAS.h5ad").exists()
            )
            self.assertTrue(
                (Path(tmpdir) / "alternative_splicing" / "test_DAS.csv").exists()
            )


if __name__ == "__main__":
    unittest.main()
