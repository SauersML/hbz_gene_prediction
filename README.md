# hbz_gene_prediction

Comprehensive HBZ/HBZP1 coding-model analysis across `copy1`, `copy2`, `copy3`, the canonical HBZ hg38 locus, and canonical HBZP1 hg38, using ab initio gene prediction, splice-aware ORF checks, comparative plotting, and AlphaGenome overlays.

## Scope
Primary questions addressed:

1. Do HBZ-region copies carry coherent multi-exon coding models?
2. Are those models ORF-valid after splicing?
3. Do independent predictors agree on structure?
4. How similar are `copy1` and `copy2` at nucleotide scale?
5. How do HBZ copy models compare to HBZ/HBZP1 hg38 references?
6. Is there supportive AlphaGenome transcription/splice signal in mapped parent-sequence windows?

## Main Inputs

- `copy1.fa` (1930 nt)
- `copy2.fa` (2176 nt)
- `copy3.fa` (2182 nt)
- `HBZ_GRCh38_NC000016_152644_154503.fa`
- `HBZP1_GRCh38_NC000016_163067_165205.fa`
- AlphaGenome parent sequence: `/Users/user/Downloads/hbz/HBZ_region_UP005hap1.fasta` (hardcoded in AlphaGenome scripts)

## Dependencies

Core Python:

- Python 3.9+
- `pandas`
- `matplotlib`
- `numpy`
- `seaborn` (for `alphagenome_hbz_multi_ontology.py`)
- `biopython` (for `transfer_hbzp1_annotation_translate.py`, uses `Bio.pairwise2`)

External tools:

- `augustus`
- `glimmerhmm`
- optional: `prodigal` (legacy/side analysis files exist, but not used by current plotting scripts)

AlphaGenome:

- package exposing `alphagenome.models.dna_client`
- environment variable: `ALPHAGENOME_API_KEY`

## Reproducible Workflow (Current Scripts)
Run from repository root.

### 1) ORF and splice-aware scan

```bash
python3 hbz_orf_splice_scan.py > hbz_orf_splice_report.tsv
```

Output:

- `hbz_orf_splice_report.tsv`

### 2) Gene prediction on copy1/2/3

```bash
augustus --species=human --gff3=on copy1.fa > copy1.augustus.gff3
augustus --species=human --gff3=on copy2.fa > copy2.augustus.gff3
augustus --species=human --gff3=on copy3.fa > copy3.augustus.gff3

glimmerhmm copy1.fa /opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human -g -o copy1.glimmer.gff3
glimmerhmm copy2.fa /opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human -g -o copy2.glimmer.gff3
glimmerhmm copy3.fa /opt/homebrew/Caskroom/miniforge/base/share/glimmerhmm/trained_dir/human -g -o copy3.glimmer.gff3
```

### 3) Spliced ORF validation from AUGUSTUS GFF3

```bash
python3 extract_spliced_orf_from_augustus.py > hbz_spliced_orf_summary.tsv
```

Output:

- `hbz_spliced_orf_summary.tsv`

Behavior note:

- The script scans `*.augustus.gff3` in the repo and expects a paired FASTA named `<stem>.fa`.
- If a paired FASTA is missing (for example `HBZP1_GRCh38.augustus.gff3` -> `HBZP1_GRCh38.fa`), that file is skipped with a stderr warning instead of failing the run.

### 4) Basic summary plots for AUGUSTUS models

```bash
python3 make_hbz_figures.py
```

Outputs:

- `hbz_gene_prediction_summary.tsv`
- `hbz_gene_models.png`
- `hbz_orf_lengths.png`

### 5) Multi-predictor panel including HBZ and HBZP1 references

```bash
python3 make_multi_predictor_gene_model_figure.py
```

Outputs:

- `hbz_multi_predictor_summary_with_hbzp1.tsv`
- `hbz_anchor_alignment_100bp.tsv`
- `hbz_all_predictors_gene_models_single_panel_with_hbzp1.png`
- `hbz_all_predictors_gene_models_single_panel_with_hbzp1.svg`

### 6) AUGUSTUS protein comparison summary (console report)

```bash
python3 compare_augustus_proteins.py
```

Note: this script prints lengths and pairwise identity to stdout and does not write a file.

### 7) Dotplot (`copy1` vs `copy2`) with AUGUSTUS exon overlays

```bash
python3 dotplot_copy1_vs_copy2.py
```

Outputs:

- `copy1_vs_copy2_dotplot_augustus_exons.png`
- `copy1_vs_copy2_dotplot_augustus_exons.svg`
- `copy1_vs_copy2_dotplot_stats.txt`

## Start-Site Follow-up

Generate sampled AUGUSTUS runs first:

```bash
augustus --species=human --gff3=on --alternatives-from-sampling=true --sample=200 copy1.fa > copy1.augustus.sampled.gff3
augustus --species=human --gff3=on --alternatives-from-sampling=true --sample=200 copy2.fa > copy2.augustus.sampled.gff3
augustus --species=human --gff3=on --alternatives-from-sampling=true --sample=200 copy3.fa > copy3.augustus.sampled.gff3
```

Then:

```bash
python3 start_site_assessment.py
python3 make_start_site_figure.py
```

Outputs:

- `start_site_assessment.tsv`
- `start_site_assessment_summary.txt`
- `start_site_assessment.png`
- `start_site_assessment.svg`

## HBZP1-Specific Workflows

### Canonical HBZP1 predictor run

```bash
python3 run_hbzp1_gene_predictor_pipeline.py
```

Outputs:

- `HBZP1_GRCh38.augustus.gff3`
- `HBZP1_GRCh38.glimmer.gff3`
- `HBZP1_GRCh38.multi_predictor_summary.tsv`
- `HBZP1_GRCh38_all_predictors_gene_models.png`
- `HBZP1_GRCh38_all_predictors_gene_models.svg`

### HBZP1 annotation transfer + translation across copies

```bash
python3 transfer_hbzp1_annotation_translate.py
```

Outputs:

- `hbzp1_annotation_transfer_translation.tsv`
- `hbzp1_annotation_transfer_predicted_proteins.fa`

## AlphaGenome Workflows

All AlphaGenome scripts currently use a hardcoded parent FASTA path:

- `/Users/user/Downloads/hbz/HBZ_region_UP005hap1.fasta`

### Scripts

- `alphagenome_ontology_probe.py`: quick ontology/output availability probe (stdout).
- `run_alphagenome_hbz.py`: single-ontology fallback run, writes `alphagenome_hbz/`.
- `alphagenome_hbz_multi_ontology.py`: enrichment across multiple ontologies, writes `alphagenome_hbz_multi_ontology/`.
- `alphagenome_visualize_hbz_region.py`: multi-panel visualization for one ontology, writes `alphagenome_hbz_visualizations/`.
- `make_gene_models_with_fetal_alphagenome_tracks.py`: final combined model + fetal/placental tracks.

### Combined fetal/embryonic figure outputs

- `hbz_gene_models_with_fetal_alphagenome_tracks.png`
- `hbz_gene_models_with_fetal_alphagenome_tracks.svg`
- `hbz_splice_sites_20bp_percent.tsv`
- `hbz_copy_fullregion_mapping.tsv`

In `hbz_copy_fullregion_mapping.tsv`, exact mapped copy windows are:

- `copy1`: `1706-3635`
- `copy2`: `12467-14642`
- `copy3`: `22485-24666`

## Current Core Result Summary

- `copy1`, `copy2`, and `copy3` each support coherent multi-exon coding models under AUGUSTUS and GlimmerHMM.
- `copy2` and `copy3` are closest at predicted protein level; `copy1` is longer.
- Dotplot evidence supports forward-orientation homology between `copy1` and `copy2`.
- AlphaGenome overlays provide context-dependent transcription/splicing support, with strongest overall signal often in `copy1` windows.

## Notes on Legacy Artifacts

The repository contains historical outputs not produced by the current script set (for example older figure name variants). For reproducibility, prefer the exact commands and file names listed above.

## Non-AlphaGenome Quick Rerun

If you want to rerun the local analysis only (no AlphaGenome API), use:

```bash
python3 hbz_orf_splice_scan.py > hbz_orf_splice_report.tsv
python3 extract_spliced_orf_from_augustus.py > hbz_spliced_orf_summary.tsv
python3 make_hbz_figures.py
python3 make_multi_predictor_gene_model_figure.py
python3 compare_augustus_proteins.py
python3 dotplot_copy1_vs_copy2.py
python3 start_site_assessment.py
python3 make_start_site_figure.py
python3 run_hbzp1_gene_predictor_pipeline.py
python3 transfer_hbzp1_annotation_translate.py
```
