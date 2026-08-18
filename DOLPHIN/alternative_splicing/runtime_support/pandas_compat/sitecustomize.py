import pandas as pd
import os
from pathlib import Path
import sqlite3


if not hasattr(pd.Series, "iteritems"):
    pd.Series.iteritems = pd.Series.items


if hasattr(pd, "DataFrame") and hasattr(pd.DataFrame, "groupby"):
    _df_groupby = pd.DataFrame.groupby

    def _patched_df_groupby(self, *args, **kwargs):
        if kwargs.get("axis", None) == 0:
            kwargs = dict(kwargs)
            kwargs.pop("axis", None)
        return _df_groupby(self, *args, **kwargs)

    pd.DataFrame.groupby = _patched_df_groupby


if hasattr(pd, "Series") and hasattr(pd.Series, "groupby"):
    _ser_groupby = pd.Series.groupby

    def _patched_ser_groupby(self, *args, **kwargs):
        if kwargs.get("axis", None) == 0:
            kwargs = dict(kwargs)
            kwargs.pop("axis", None)
        return _ser_groupby(self, *args, **kwargs)

    pd.Series.groupby = _patched_ser_groupby


BEDTOOLS_BIN = os.environ.get("DOLPHIN_AS_BEDTOOLS_BIN")
if BEDTOOLS_BIN:
    BEDTOOLS_BIN = str(Path(BEDTOOLS_BIN).expanduser())

if BEDTOOLS_BIN and Path(BEDTOOLS_BIN).is_dir():
    os.environ["PATH"] = BEDTOOLS_BIN + os.pathsep + os.environ.get("PATH", "")
    try:
        from pybedtools import helpers as _pybedtools_helpers

        _pybedtools_helpers.set_bedtools_path(BEDTOOLS_BIN)
    except Exception:
        pass


_sqlite_connect = sqlite3.connect


def _patched_sqlite_connect(*args, **kwargs):
    kwargs.setdefault("timeout", 300)
    return _sqlite_connect(*args, **kwargs)


sqlite3.connect = _patched_sqlite_connect


try:
    from gffutils import interface as _gffutils_interface

    _featuredb_update = _gffutils_interface.FeatureDB.update

    def _patched_featuredb_update(self, *args, **kwargs):
        kwargs.setdefault("pragmas", self.pragmas)
        return _featuredb_update(self, *args, **kwargs)

    _gffutils_interface.FeatureDB.update = _patched_featuredb_update
except Exception:
    pass
