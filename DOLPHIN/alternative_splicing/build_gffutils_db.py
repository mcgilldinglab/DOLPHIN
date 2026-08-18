import argparse
import os
import time
from pathlib import Path


DEFAULT_GTF = os.environ.get("DOLPHIN_AS_GTF_PATH")
DEFAULT_DB = os.environ.get("DOLPHIN_AS_GFFUTILS_DB")


def build_db(gtf_path, db_path):
    from outrigger.io.gtf import create_db

    gtf_path = Path(gtf_path).expanduser().resolve()
    db_path = Path(db_path).expanduser().resolve()
    db_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    db = create_db(str(gtf_path), str(db_path))
    elapsed = time.time() - start
    print(f"Built gffutils DB at {db_path}")
    print(f"Elapsed seconds: {elapsed:.2f}")
    return db


def main():
    parser = argparse.ArgumentParser(description="Build a reusable gffutils DB for the AS pipeline.")
    parser.add_argument(
        "--gtf-path",
        default=DEFAULT_GTF,
        required=DEFAULT_GTF is None,
        help="Input GTF path (or set DOLPHIN_AS_GTF_PATH).",
    )
    parser.add_argument(
        "--db-path",
        default=DEFAULT_DB,
        required=DEFAULT_DB is None,
        help="Output gffutils DB path (or set DOLPHIN_AS_GFFUTILS_DB).",
    )
    args = parser.parse_args()

    gtf_path = Path(args.gtf_path).expanduser()
    db_path = Path(args.db_path).expanduser()

    if not gtf_path.exists():
        raise FileNotFoundError(f"Missing GTF: {gtf_path}")

    build_db(gtf_path, db_path)


if __name__ == "__main__":
    main()
