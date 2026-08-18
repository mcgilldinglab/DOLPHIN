---
name: dolphin-pipeline
description: "Use for DOLPHIN repository work involving the finalized GitHub-ready pipeline: raw full-length or 10x read preprocessing, graph generation, lightweight graph-store model input, model embedding export, alternative splicing BAM or junction routes, EDEG/JDEG historical one-vs-rest MAST analysis, pipeline output contracts, validation, runtime environment setup, or preparing DOLPHIN docs/skills from the clean package."
---

# DOLPHIN Pipeline

## Core Rule

Treat the Git repository containing `.agents/skills/dolphin-pipeline` as the source of truth for the finalized engineered DOLPHIN package. When the skill is repository-scoped, resolve the package root as three directories above this skill folder. When it is installed at user scope, locate the active DOLPHIN repository from the current working directory unless the user provides another root.

Do not resurrect legacy code paths from `DOLPHIN_original_clone`, `DOLPHIN_optimized`, `AS_run`, smoke-test folders, or execution-run folders unless the user explicitly asks for historical comparison or debugging.

## First Checks

Before editing or running the pipeline:

1. Locate the DOLPHIN package root.
2. Read the relevant package docs from `docs/` and the module README before changing code.
3. Confirm the task category: preprocessing, graph/model, alternative splicing, EDEG/JDEG, validation, docs, or GitHub packaging.
4. Preserve the final output contract unless the user explicitly approves a contract change.

## Ask Before Proceeding

Ask the user before proceeding when the answer changes scientific meaning, resource use, or retained outputs:

- The modality is ambiguous: full-length vs 10x.
- The AS route is ambiguous: BAM aggregation vs direct junction aggregation.
- A run may overwrite, delete, or move existing results.
- A full run may consume substantial time, GPU, scratch, or disk space.
- The requested optimization may change output values, event definitions, aggregation semantics, or the final file contract.
- The task would expose pairwise MAST, modify the historical one-vs-rest EDEG/JDEG route, or change the EDEG input from `FeatureComp_<out_name>.h5ad`.
- The task would change Outrigger event/index/PSI semantics rather than only its execution environment.
- The durable output root or scratch/work root is not inferable from the current project context.
- The user asks to publish, package, install globally, or upload to GitHub.

Do not ask for confirmation when the missing value can be safely inferred from filenames, current working directory, existing manifests, or the user's immediately preceding instruction. State the assumption and continue.

## Reference Selection

Read [references/pipeline-contract.md](references/pipeline-contract.md) for any code, docs, or workflow task touching pipeline logic, output files, modality differences, or module boundaries.

Read [references/environment-and-validation.md](references/environment-and-validation.md) before running commands, debugging failed runs, optimizing speed, comparing BAM vs junction routes, or preparing user-facing installation instructions.

## Non-Negotiable Contracts

- Full-length and 10x should follow the same high-level logic; modality-specific code should only adapt input formats.
- Smoke tests and full runs must use the same method and same route logic. Limit cells with a subset parameter, not by switching to a different algorithm.
- The final efficient graph/model path uses grouped matrices, pooled H5ADs, and lightweight graph-store inputs. Do not make per-cell `*_fea.csv`, `*_adj.csv`, or per-cell graph H5AD fragments the retained public output.
- Prefer graph-store model input over mandatory large `.pt` serialization. Keep `.pt` only as optional backward compatibility.
- Alternative splicing keeps both BAM aggregation and direct junction aggregation routes for both full-length and 10x.
- Do not change Outrigger event/index/PSI semantics. Runtime patches may stabilize execution, but must not redefine what Outrigger calculates.
- Never fall back to repository-author input, reference, tool, result, or scratch paths. Resolve AS paths from explicit CLI arguments, documented `DOLPHIN_AS_*` environment variables, the active `PATH`, or safe system temporary-directory discovery.
- EDEG/JDEG public workflow is historical one-cluster-vs-rest MAST only. Do not expose pairwise MAST as the default public route.

## Implementation Guidance

When optimizing code, preserve output identity first and speed second. Record baseline timing, optimized timing, exact inputs, route, cell subset, and comparison metrics. For matrix-like outputs, compare shapes, cell/event overlap, exact equality where possible, and correlation where exact equality is not expected due to route differences.

When debugging failed runs, avoid opaque wrappers if they hide stderr. Re-run the failing command directly with visible stdout/stderr, then patch the DOLPHIN orchestration layer only if the fix does not change the validated scientific logic.
