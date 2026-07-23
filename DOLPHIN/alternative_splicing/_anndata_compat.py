def enable_nullable_string_writes():
    try:
        import anndata

        settings = getattr(anndata, "settings", None)
        if settings is not None and hasattr(settings, "allow_write_nullable_strings"):
            settings.allow_write_nullable_strings = True
    except Exception:
        pass


def write_h5ad_preserve_strings(adata, path, **kwargs):
    try:
        from anndata._io.h5ad import write_h5ad

        write_h5ad(path, adata, convert_strings_to_categoricals=False, **kwargs)
    except Exception:
        adata.write(path, **kwargs)
