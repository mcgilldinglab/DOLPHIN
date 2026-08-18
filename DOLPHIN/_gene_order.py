import re


_LEGACY_GENE_ID = re.compile(r"^gene([0-9]+)$")


def _unique_gene_ids(gene_ids):
    """Return unique gene IDs as strings while preserving first appearance."""
    unique_ids = []
    seen = set()
    for value in gene_ids:
        gene_id = str(value)
        if gene_id not in seen:
            seen.add(gene_id)
            unique_ids.append(gene_id)
    return unique_ids


def build_gene_order_map(gene_ids):
    """Build a deterministic gene order without assuming a species-specific ID format.

    Legacy DOLPHIN references use IDs such as ``gene2`` and ``gene10`` and
    historically sorted them by the numeric suffix. Preserve that ordering when
    every ID follows the legacy format. For other annotations, preserve the
    reference order represented by first appearance.
    """
    unique_ids = _unique_gene_ids(gene_ids)
    legacy_matches = [_LEGACY_GENE_ID.fullmatch(gene_id) for gene_id in unique_ids]

    if unique_ids and all(match is not None for match in legacy_matches):
        ordered_ids = [
            gene_id
            for _, gene_id in sorted(
                enumerate(unique_ids),
                key=lambda item: (
                    int(legacy_matches[item[0]].group(1)),
                    item[0],
                ),
            )
        ]
    else:
        ordered_ids = unique_ids

    return {gene_id: order for order, gene_id in enumerate(ordered_ids)}


def map_gene_orders(gene_ids, reference_gene_ids):
    """Map gene IDs to the order defined by a canonical reference sequence."""
    order_map = build_gene_order_map(reference_gene_ids)
    normalized_ids = [str(gene_id) for gene_id in gene_ids]
    missing = sorted(set(normalized_ids).difference(order_map))
    if missing:
        preview = ", ".join(missing[:5])
        raise ValueError(
            "Gene IDs were not found in the reference gene order: "
            f"{preview}"
        )
    return [order_map[gene_id] for gene_id in normalized_ids]
