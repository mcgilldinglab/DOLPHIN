import argparse
import os
import site
import sys
import time
from pathlib import Path


DEFAULT_GTF = "/mnt/md0/kailu/microbiota_fungi_codex/04_human_gtf_star_index/inputs/Homo_sapiens.GRCh38.107.gtf"
DEFAULT_DB = "/mnt/md0/kailu/DOLPHIN_codex/AS_run/reference_cache/Homo_sapiens.GRCh38.107.gtf.db"
AS_SITE = "/mnt/md0/kailu/DOLPHIN_codex/runtime_envs/as_site_manual"


def build_db(gtf_path, db_path):
    site.addsitedir(AS_SITE)
    from outrigger.io.gtf import create_db

    db_path.parent.mkdir(parents=True, exist_ok=True)
    start = time.time()
    db = create_db(str(gtf_path), str(db_path))
    elapsed = time.time() - start
    print(f"Built gffutils DB at {db_path}")
    print(f"Elapsed seconds: {elapsed:.2f}")
    return db


def main():
    parser = argparse.ArgumentParser(description="Build a reusable gffutils DB for the AS pipeline.")
    parser.add_argument("--gtf-path", default=DEFAULT_GTF)
    parser.add_argument("--db-path", default=DEFAULT_DB)
    args = parser.parse_args()

    gtf_path = Path(args.gtf_path)
    db_path = Path(args.db_path)

    if not gtf_path.exists():
        raise FileNotFoundError(f"Missing GTF: {gtf_path}")

    build_db(gtf_path, db_path)


if __name__ == "__main__":
    main()
