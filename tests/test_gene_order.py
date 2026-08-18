import unittest

from DOLPHIN._gene_order import build_gene_order_map, map_gene_orders


class GeneOrderTests(unittest.TestCase):
    def test_legacy_gene_ids_keep_numeric_order(self):
        gene_ids = ["gene10", "gene2", "gene1", "gene2"]

        self.assertEqual(
            build_gene_order_map(gene_ids),
            {"gene1": 0, "gene2": 1, "gene10": 2},
        )

    def test_legacy_subset_mapping_matches_numeric_sort(self):
        reference_ids = ["gene10", "gene2", "gene1"]

        self.assertEqual(
            map_gene_orders(["gene10", "gene1"], reference_ids),
            [2, 0],
        )

    def test_zebrafish_gene_ids_keep_reference_order(self):
        gene_ids = [
            "gene0304170_df_nrg",
            "ENSDARG00000012345",
            "LOC12345",
        ]

        self.assertEqual(
            build_gene_order_map(gene_ids),
            {
                "gene0304170_df_nrg": 0,
                "ENSDARG00000012345": 1,
                "LOC12345": 2,
            },
        )

    def test_mixed_gene_ids_do_not_partially_parse_numeric_ids(self):
        gene_ids = ["gene10", "gene2", "gene0304170_df_nrg"]

        self.assertEqual(
            build_gene_order_map(gene_ids),
            {"gene10": 0, "gene2": 1, "gene0304170_df_nrg": 2},
        )

    def test_mapping_uses_reference_order(self):
        reference_ids = ["gene0304170_df_nrg", "ENSDARG00000012345"]

        self.assertEqual(
            map_gene_orders(
                ["ENSDARG00000012345", "gene0304170_df_nrg"],
                reference_ids,
            ),
            [1, 0],
        )

    def test_missing_gene_id_has_clear_error(self):
        with self.assertRaisesRegex(ValueError, "not found in the reference"):
            map_gene_orders(["missing_gene"], ["gene1"])


if __name__ == "__main__":
    unittest.main()
