import anndata
import numpy as np
import pandas as pd
import tempfile
import unittest
from pathlib import Path

from DOLPHIN.graph_generation._anndata_compat import write_h5ad_preserve_strings


class AnnDataCompatibilityTests(unittest.TestCase):
    def test_write_h5ad_enables_nullable_string_support(self):
        previous = anndata.settings.allow_write_nullable_strings
        anndata.settings.allow_write_nullable_strings = False
        try:
            obs = pd.DataFrame(
                {"label": pd.array(["alpha"], dtype="string")},
                index=pd.Index(pd.array(["cell-1"], dtype="string")),
            )
            var = pd.DataFrame(
                index=pd.Index(pd.array(["feature-1"], dtype="string"))
            )
            adata = anndata.AnnData(X=np.ones((1, 1)), obs=obs, var=var)
            with tempfile.TemporaryDirectory() as tmpdir:
                output = Path(tmpdir) / "nullable_strings.h5ad"
                write_h5ad_preserve_strings(adata, output)
                restored = anndata.read_h5ad(output)

            self.assertEqual(restored.obs_names.tolist(), ["cell-1"])
            self.assertEqual(restored.obs["label"].tolist(), ["alpha"])
        finally:
            anndata.settings.allow_write_nullable_strings = previous


if __name__ == "__main__":
    unittest.main()
