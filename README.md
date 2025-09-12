<p align="center">
	<img src="https://opencobra.github.io/cobrapy/_static/img/cobrapy_logo.png" alt="COBRApy logo" width="120"/>
</p>

# COBRAxy — Metabolic analysis and visualization toolkit (Galaxy-ready)

COBRAxy (COBRApy in Galaxy) is a toolkit to compute, analyze, and visualize metabolism at the reaction level from transcriptomics and metabolomics data. It enables users to:

- derive Reaction Activity Scores (RAS) from gene expression and Reaction Propensity Scores (RPS) from metabolite abundances,
- integrate RAS into model bounds,
- perform flux sampling with either CBS (constraint-based sampling) or OPTGP,
- compute statistics (pFBA, FVA, sensitivity) and generate styled SVG/PDF metabolic maps,
- run all tools as Galaxy wrappers or via CLI on any machine.

It extends the MaREA 2 (Metabolic Reaction Enrichment Analysis) concept by adding sampling-based flux comparison and rich visualization. The repository ships both Python CLIs and Galaxy tool XMLs.

## Table of contents

- Overview and features
- Requirements
- Installation (pip/conda)
- Quick start (CLI)
- Tools and usage
	- custom_data_generator
	- ras_generator (RAS)
	- rps_generator (RPS)
	- ras_to_bounds
	- flux_simulation (CBS/OPTGP)
	- marea (enrichment + maps)
	- flux_to_map (maps from fluxes)
	- marea_cluster (clustering auxiliaries)
- Typical workflow
- Input/output formats
- Galaxy usage
- Troubleshooting
- Contributing
- License and citations
- Useful links

## Overview and features

COBRAxy builds on COBRApy to deliver end‑to‑end analysis from expression/metabolite data to flux statistics and map rendering:

- RAS and RPS computation from tabular inputs
- Bounds integration and model preparation
- Flux sampling: CBS (GLPK backend) with automatic fallback to a COBRApy interface, or OPTGP
- Flux statistics: mean/median/quantiles, pFBA, FVA, sensitivity
- Map styling/export: SVG with optional PDF/PNG export
- Ready-made Galaxy wrappers for all tools

Bundled resources in `local/` include example models (ENGRO2, Recon), gene mappings, a default medium, and SVG maps.

## Requirements

- OS: Linux, macOS, or Windows (Linux recommended; Galaxy typically runs on Linux)
- Python: 3.8.20 ≤ version < 3.12 (as per `setup.py`)
- Python packages (installed automatically by `pip install .`):
	- cobra==0.29.0, numpy==1.24.4, pandas==2.0.3, scipy==1.11, scikit-learn==1.3.2, seaborn==0.13.0
	- matplotlib==3.7.3, lxml==5.2.2, cairosvg==2.7.1, svglib==1.5.1, pyvips==2.2.3, Pillow
	- joblib==1.4.2, anndata==0.8.0, pydeseq2==0.5.1
- Optional but recommended for CBS sampling performance:
	- GLPK solver and Python bindings
		- System library: glpk (e.g., Ubuntu: `apt-get install glpk-utils libglpk40`)
		- Python: `swiglpk` (note: CBS falls back to a COBRApy interface if GLPK is unavailable)
- For pyvips: system libvips (e.g., Ubuntu: `apt-get install libvips`)

Notes:
- If you hit system-level library errors for SVG/PDF/PNG conversion or vips, install the corresponding OS packages.
- GPU is not required.

## Installation

Python virtual environment is strongly recommended.

### Install from source (pip)

1) Clone the repo and install:

```bash
git clone https://github.com/CompBtBs/COBRAxy.git
cd COBRAxy
python3 -m venv .venv && source .venv/bin/activate
pip install --upgrade pip
pip install .
```

This installs console entry points: `custom_data_generator`, `ras_generator`, `rps_generator`, `ras_to_bounds`, `flux_simulation`, `flux_to_map`, `marea`, `marea_cluster`.

### Install with conda (alternative)

```bash
conda create -n cobraxy python=3.10 -y
conda activate cobraxy
pip install .
# Optional system deps (Ubuntu): sudo apt-get install libvips libxml2 libxslt1.1 glpk-utils
# Optional Python bindings for GLPK: pip install swiglpk
```

## Quick start (CLI)

All tools provide `-h/--help` for details. Outputs are TSV/CSV and SVG/PDF files depending on the tool and flags.

Example minimal flow (using built-in ENGRO2 model and provided assets):

```bash
# 1) Generate rules/reactions/bounds/medium from a model (optional if using bundled ones)
custom_data_generator \
	-id local/models/ENGRO2.xml \
	-mn ENGRO2.xml \
	-orules out/ENGRO2_rules.tsv \
	-orxns out/ENGRO2_reactions.tsv \
	-omedium out/ENGRO2_medium.tsv \
	-obnds out/ENGRO2_bounds.tsv

# 2) Compute RAS from expression data
ras_generator \
	-td $(pwd) \
	-in my_expression.tsv \
	-ra out/ras.tsv \
	-rs ENGRO2

# 3) Integrate RAS into bounds
ras_to_bounds \
	-td $(pwd) \
	-ms ENGRO2 \
	-ir out/ras.tsv \
	-rs true \
	-idop out/ras_bounds

# 4) Flux sampling (CBS)
flux_simulation \
	-td $(pwd) \
	-ms ENGRO2 \
	-in out/ras_bounds/sample1.tsv,out/ras_bounds/sample2.tsv \
	-ni sample1,sample2 \
	-a CBS -ns 500 -sd 0 -nb 1 \
	-ot mean,median,quantiles \
	-ota pFBA,FVA,sensitivity \
	-idop out/flux

# 5) Enrichment + map styling (RAS/RPS or fluxes)
marea \
	-td $(pwd) \
	-using_RAS true -input_data out/ras.tsv \
	-comparison manyvsmany -test ks \
	-generate_svg true -generate_pdf true \
	-choice_map ENGRO2 -idop out/maps
```

## Tools and usage

Below is a high‑level summary of each CLI. Use `--help` for the full list of options.

### 1) custom_data_generator

Generate model‑derived assets.

Required inputs:
- `-id/--input`: model file (XML or JSON; gz/zip/bz2 also supported via extension)
- `-mn/--name`: the original file name including extension (Galaxy renames files; this preserves the true format)
- `-orules`, `-orxns`, `-omedium`, `-obnds`: output paths

Outputs:
- TSV with rules, reactions, exchange medium, and bounds.

### 2) ras_generator (Reaction Activity Scores)

Compute RAS from a gene expression table.

Key inputs:
- `-td/--tool_dir`: repository root path (used to locate `local/` assets)
- `-in/--input`: expression TSV (rows: genes; columns: samples)
- `-rs/--rules_selector`: model/rules choice, e.g. `ENGRO2` or `Custom` with `-rl` and `-rn`
- Optional: `-rl/--rule_list` custom rules TSV, `-rn/--rules_name` its original name/extension
- Output: `-ra/--ras_output` TSV

### 3) rps_generator (Reaction Propensity Scores)

Compute RPS from a metabolite abundance table.

Key inputs:
- `-td/--tool_dir`: repository root
- `-id/--input`: metabolite TSV (rows: metabolites; columns: samples)
- `-rc/--reaction_choice`: `default` or `custom` with `-cm/--custom` reactions TSV
- Output: `-rp/--rps_output` TSV

### 4) ras_to_bounds

Integrate RAS into reaction bounds for a given model and medium.

Key inputs:
- `-td/--tool_dir`: repository root
- `-ms/--model_selector`: one of `ENGRO2` or `Custom` with `-mo/--model` and `-mn/--model_name`
- Medium: `-mes/--medium_selector` (default `allOpen`) or `-meo/--medium` custom TSV
- RAS: `-ir/--input_ras` and `-rs/--ras_selector` (true/false)
- Output folder: `-idop/--output_path`

Outputs:
- One bounds TSV per sample in the RAS table.

### 5) flux_simulation

Flux sampling with CBS or OPTGP and downstream statistics.

Key inputs:
- `-td/--tool_dir`
- Model: `-ms/--model_selector` (ENGRO2 or Custom with `-mo`/`-mn`)
- Bounds files: `-in` (comma‑separated list) and `-ni/--names` (comma‑separated sample names)
- Algorithm: `-a CBS|OPTGP`; CBS uses GLPK if available and falls back to a COBRApy interface
- Sampling params: `-ns/--n_samples`, `-th/--thinning` (OPTGP), `-nb/--n_batches`, `-sd/--seed`
- Outputs: `-ot/--output_type` (mean,median,quantiles) and `-ota/--output_type_analysis` (pFBA,FVA,sensitivity)
- Output path: `-idop/--output_path`

Outputs:
- Per‑sample or aggregated CSV/TSV with flux samples and statistics.

### 6) marea

Statistical enrichment and map styling for RAS and/or RPS groups with optional DESeq2‑style testing via `pydeseq2`.

Key inputs:
- `-td/--tool_dir`
- Comparison: `-co manyvsmany|onevsrest|onevsmany`
- Test: `-te ks|ttest_p|ttest_ind|wilcoxon|mw|DESeq`
- Thresholds: `-pv`, `-adj` (FDR), `-fc`
- Data: RAS `-using_RAS` plus `-input_data` or multiple datasets with names; similarly for RPS with `-using_RPS`
- Map: `-choice_map HMRcore|ENGRO2|Custom` or `-custom_map` SVG
- Output: `-gs/--generate_svg`, `-gp/--generate_pdf`, output dir `-idop`

Outputs:
- Styled SVG (and optional PDF/PNG) highlighting enriched reactions by color/width per your thresholds.

### 7) flux_to_map

Like `marea`, but driven by fluxes instead of RAS/RPS. Accepts single or multiple flux datasets and produces styled maps.

### 8) marea_cluster

Convenience clustering utilities (k‑means, DBSCAN, hierarchical) for grouping samples; produces labels and optional plots.

## Typical workflow

1. Prepare a model and generate its assets (optional if using bundled assets): `custom_data_generator`
2. Compute RAS from expression: `ras_generator` (and/or compute RPS via `rps_generator`)
3. Integrate RAS into bounds: `ras_to_bounds`
4. Sample fluxes: `flux_simulation` with CBS or OPTGP
5. Analyze and visualize: `marea` or `flux_to_map` to render SVG/PDF metabolic maps
6. Optionally cluster or further analyze results: `marea_cluster`

## Input/output formats

Unless otherwise stated, inputs are tab‑separated (TSV) text files with headers.

- Expression (RAS): rows = genes (HGNC/Ensembl/symbol/Entrez supported), columns = samples
- Metabolite table (RPS): rows = metabolites, columns = samples
- Rules/Reactions: TSV with two columns: ReactionID, Rule/Reaction
- Bounds: TSV with index = reaction IDs, columns = lower_bound, upper_bound
- Medium: single‑column TSV listing exchange reactions
- Flux samples/statistics: CSV/TSV with reactions as rows and samples/statistics as columns

## Galaxy usage

Each CLI has a corresponding Galaxy tool XML in the repository (e.g., `marea.xml`, `flux_simulation.xml`). Use `shed.yml` to publish to a Galaxy toolshed. The `local/` directory provides models, mappings, and maps for out‑of‑the‑box runs inside Galaxy.

## Troubleshooting

- GLPK/CBS issues: if `swiglpk` or GLPK is missing, `flux_simulation` will attempt a COBRApy fallback. Install GLPK + `swiglpk` for best performance.
- pyvips errors: install `libvips` on your system. Reinstall the `pyvips` wheel afterward if needed.
- PDF/SVG conversions: ensure `cairosvg`, `svglib`, and system libraries (`libxml2`, `libxslt`) are installed.
- Python version: stick to Python ≥3.8.20 and <3.12.
- Memory/time: reduce `-ns` (samples) or `-nb` (batches); consider OPTGP if CBS is slow for your model.

## Contributing

Pull requests are welcome. Please:
- keep changes focused and documented,
- add concise docstrings/comments in English,
- preserve public CLI parameters and file formats.

## License and citations

This project is distributed under the MIT License. If you use COBRAxy in academic work, please cite COBRApy and MaREA, and reference this repository.

## Useful links

- COBRAxy Google Summer of Code 2024: https://summerofcode.withgoogle.com/programs/2024/projects/LSrCKfq7
- COBRApy: https://opencobra.github.io/cobrapy/
- MaREA4Galaxy: https://galaxyproject.org/use/marea4galaxy/
- Galaxy project: https://usegalaxy.org/
