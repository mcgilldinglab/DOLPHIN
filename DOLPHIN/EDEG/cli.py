import argparse

from .call_convert import run_h5ad_rds
from .call_MAST import call_MAST_func
from .generate_EDEG import run_edeg
from .generate_JDEG import run_jdeg


def build_parser():
    parser = argparse.ArgumentParser(description="DOLPHIN exon/junction differential analysis utilities.")
    sub = parser.add_subparsers(dest="command", required=True)

    p_convert = sub.add_parser("convert", help="Convert h5ad to Seurat rds")
    p_convert.add_argument("--input-anndata", required=True)
    p_convert.add_argument("--output-rds", required=True)

    p_mast = sub.add_parser("mast", help="Run historical Seurat/MAST in one-cluster-vs-rest mode")
    p_mast.add_argument("--input-rds", required=True)
    p_mast.add_argument("--output-dir", required=True)
    p_mast.add_argument("--ident-column", default="cluster")
    p_mast.add_argument("--cluster-values", default="")
    p_mast.add_argument("--logfc-threshold", type=float, default=0.0)
    p_mast.add_argument("--skip-normalize", action="store_true")
    p_mast.add_argument("--rscript-binary", default="Rscript")

    p_edeg = sub.add_parser("edeg", help="Aggregate exon-level MAST output to gene-level EDEG")
    p_edeg.add_argument("--seurat-output", required=True)
    p_edeg.add_argument("--adata-input", required=True)
    p_edeg.add_argument("--gtf-path", required=True)
    p_edeg.add_argument("--output", required=True)

    p_jdeg = sub.add_parser("jdeg", help="Aggregate junction-level MAST output to gene-level JDEG")
    p_jdeg.add_argument("--seurat-output", required=True)
    p_jdeg.add_argument("--output", required=True)

    return parser


def main():
    args = build_parser().parse_args()
    if args.command == "convert":
        run_h5ad_rds(args.input_anndata, args.output_rds)
    elif args.command == "mast":
        cluster_values = [v for v in args.cluster_values.split(",") if v] if args.cluster_values else None
        call_MAST_func(
            input_rds=args.input_rds,
            output_dir=args.output_dir,
            ident_column=args.ident_column,
            cluster_values=cluster_values,
            logfc_threshold=args.logfc_threshold,
            normalize=not args.skip_normalize,
            rscript_binary=args.rscript_binary,
        )
    elif args.command == "edeg":
        run_edeg(
            seurat_output=args.seurat_output,
            adata_input=args.adata_input,
            gtf_path=args.gtf_path,
            output=args.output,
        )
    elif args.command == "jdeg":
        run_jdeg(seurat_output=args.seurat_output, output=args.output)


if __name__ == "__main__":
    main()
