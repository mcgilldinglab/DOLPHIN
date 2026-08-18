import os
import stat
import sys
import tempfile
import unittest
from pathlib import Path

from DOLPHIN.alternative_splicing import pipeline
from DOLPHIN.alternative_splicing.presets import FULL_LENGTH_BAM_DEFAULTS
from DOLPHIN.alternative_splicing.presets import FULL_LENGTH_JUNCTION_DEFAULTS
from DOLPHIN.alternative_splicing.presets import TENX_BAM_DEFAULTS
from DOLPHIN.alternative_splicing.presets import TENX_JUNCTION_DEFAULTS


class AlternativeSplicingConfigurationTests(unittest.TestCase):
    def _base_overrides(self, route_defaults):
        return {
            **route_defaults,
            "embedding_h5ad": None,
            "metadata_path": None,
            "bam_root": None,
            "junction_root": None,
            "gtf_path": None,
            "gffutils_db": None,
            "genome_sizes_path": None,
            "fasta_path": None,
            "star_index_dir": None,
            "rg_bam_path": None,
        }

    def _configured_args(self, root, route_defaults, include_star=False):
        root = Path(root)
        bam_root = root / "bam"
        junction_root = root / "junction"
        bedtools_bin = root / "bin"
        for directory in (bam_root, junction_root, bedtools_bin):
            directory.mkdir()

        bedtools = bedtools_bin / "bedtools"
        bedtools.write_text("#!/bin/sh\nexit 0\n")
        bedtools.chmod(bedtools.stat().st_mode | stat.S_IXUSR)

        files = {}
        for name in ("embedding.h5ad", "metadata.tsv", "genes.gtf", "genome.sizes", "genome.fa"):
            path = root / name
            path.touch()
            files[name] = path

        argv = [
            "--embedding-h5ad",
            str(files["embedding.h5ad"]),
            "--metadata-path",
            str(files["metadata.tsv"]),
            "--bam-root",
            str(bam_root),
            "--junction-root",
            str(junction_root),
            "--gtf-path",
            str(files["genes.gtf"]),
            "--genome-sizes-path",
            str(files["genome.sizes"]),
            "--fasta-path",
            str(files["genome.fa"]),
            "--samtools-binary",
            sys.executable,
            "--outrigger-python",
            sys.executable,
            "--bedtools-bin-dir",
            str(bedtools_bin),
            "--results-root",
            str(root / "results"),
            "--logs-root",
            str(root / "logs"),
            "--prepared-inputs-root",
            str(root / "prepared"),
            "--outrigger-work-root",
            str(root / "outrigger_work"),
        ]
        if include_star:
            star_index = root / "star_index"
            star_index.mkdir()
            argv.extend(
                [
                    "--star-index-dir",
                    str(star_index),
                    "--star-binary",
                    sys.executable,
                ]
            )

        parser = pipeline.build_parser(
            default_overrides=self._base_overrides(route_defaults)
        )
        return parser.parse_args(argv)

    def test_missing_inputs_report_cli_and_environment_options(self):
        parser = pipeline.build_parser(
            default_overrides=self._base_overrides(FULL_LENGTH_JUNCTION_DEFAULTS)
        )
        args = parser.parse_args([])

        with self.assertRaisesRegex(
            pipeline.ASConfigurationError,
            r"--embedding-h5ad \(or DOLPHIN_AS_EMBEDDING_H5AD\)",
        ):
            pipeline.validate_runtime_config(args)

    def test_junction_route_does_not_require_star_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = self._configured_args(
                tmpdir,
                FULL_LENGTH_JUNCTION_DEFAULTS,
            )
            validated = pipeline.validate_runtime_config(args)

        self.assertEqual(validated.aggregation_mode, "junction")
        self.assertIsNone(validated.star_index_dir)

    def test_bam_route_requires_star_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = self._configured_args(tmpdir, FULL_LENGTH_BAM_DEFAULTS)
            with self.assertRaisesRegex(
                pipeline.ASConfigurationError,
                r"--star-index-dir \(or DOLPHIN_AS_STAR_INDEX_DIR\)",
            ):
                pipeline.validate_runtime_config(args)

    def test_bam_route_accepts_user_configured_star_index(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            args = self._configured_args(
                tmpdir,
                FULL_LENGTH_BAM_DEFAULTS,
                include_star=True,
            )
            validated = pipeline.validate_runtime_config(args)

        self.assertEqual(validated.aggregation_mode, "bam")
        self.assertTrue(os.path.isabs(validated.star_index_dir))

    def test_executable_sources_have_no_machine_specific_paths(self):
        source_paths = [
            Path(pipeline.__file__),
            Path(pipeline.__file__).with_name("presets.py"),
            Path(pipeline.__file__).with_name("build_gffutils_db.py"),
            Path(pipeline.__file__).parent
            / "runtime_support"
            / "pandas_compat"
            / "sitecustomize.py",
        ]
        for source_path in source_paths:
            source = source_path.read_text()
            self.assertNotIn("/mnt/md0", source)
            self.assertNotIn("/home/", source)

    def test_route_presets_only_describe_modality_and_method(self):
        user_path_keys = {
            "embedding_h5ad",
            "metadata_path",
            "bam_root",
            "junction_root",
            "gtf_path",
            "gffutils_db",
            "genome_sizes_path",
            "fasta_path",
            "star_index_dir",
            "rg_bam_path",
        }
        route_presets = (
            FULL_LENGTH_BAM_DEFAULTS,
            FULL_LENGTH_JUNCTION_DEFAULTS,
            TENX_BAM_DEFAULTS,
            TENX_JUNCTION_DEFAULTS,
        )

        for preset in route_presets:
            self.assertTrue(user_path_keys.isdisjoint(preset))


if __name__ == "__main__":
    unittest.main()
