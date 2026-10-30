# clinicalAssociation.R — Tutorial

This app is the Shiny interface for the standalone R pipeline for computing gene signature activity scores and testing their association with clinical variables.

For larger datasets, download the app from github (link coming soon) and then run the app locally in RStudio.

---

## Installation

```r
# Install required packages if not already present
install.packages(c("pheatmap", "ggplot2"))
```

Then source the script:

```r
source("clinicalAssociation.R")
```

---

## Input File Formats

### Signature files (one file per signature)

Tab-separated. Line 1 is the signature name; subsequent lines are gene + log2FC pairs.  
Positive log2FC → up-regulated gene. Negative → down-regulated.

```
MySignature
GENE_A	2.45
GENE_B	1.10
GENE_C	-1.80
GENE_D	-3.22
```

### Expression CSV

Genes as rows, samples as columns. First column = gene symbols (row names).

```
Gene,Sample1,Sample2,Sample3
TP53,8.23,7.91,9.14
MYC,12.40,11.88,13.02
BRCA1,6.55,7.20,6.80
```

### Clinical CSV

Samples as rows, variables as columns. First column = sample IDs.  
Numeric and categorical variables can be in the same file — they are auto-detected.

```
SampleID,Age,Tumor_Stage,Gleason_Score,Subtype
Sample1,62,II,7,Luminal
Sample2,55,III,8,Basal
Sample3,71,I,6,Luminal
```

---

## Function Parameters

The underlying `run_clinical_association()` function (shown in full below) takes more parameters than the Shiny app exposes. Here's what the app's sidebar actually lets you control:

| Shiny App Input | Underlying Parameter | Description |
|---|---|---|
| Signature files (upload, multiple) | `sig_folder` | One or more signature files; each: line 1 = signature name, remaining lines = `GENE<TAB>log2FC` |
| Expression CSV (upload) | `expression` | Expression matrix, genes × samples |
| Clinical CSV (upload) | `clinical` | Clinical metadata, samples × variables; numeric vs. categorical auto-detected |
| Minimum valid samples per variable | `min_samples` | Minimum samples with valid values needed to test a clinical variable (app default: 6) |
| Heatmap size (inches) | `figure_size` | Width and height of summary heatmap PDFs in inches (app default: 20) |

The remaining parameters aren't exposed as separate controls because the app handles them automatically: `sig_pattern` matches every uploaded signature file regardless of extension; `output_root` is a temporary working location you never need to see — results are packaged into the downloadable zip when the run finishes; `sample_list` isn't needed since the function already restricts each clinical variable's z-scores and activity scores to just the samples with a valid value for that variable (see Notes below); and `support_plots`/`save_temp` are always left on, so every run includes the full set of per-signature plots and intermediate files in the download.

---


### Step-by-step Process

**1. Load signatures**  
Each file in `sig_folder` matching `sig_pattern` is read. Line 1 is the signature name. Subsequent lines split into gene symbol and log2FC. Genes with log2FC > 0 go into `up_genes`; log2FC ≤ 0 go into `down_genes`. Gene symbols are uppercased throughout.

**2. Load expression**  
The CSV is read with the first column as rownames (gene symbols, uppercased). The result is a numeric matrix of genes × samples.

**3. Auto-detect clinical variable types**  
The clinical CSV is loaded. For each variable column, the function strips NA-like strings (`"NA"`, `"[Not Available]"`, `"Unknown"`, etc.) and attempts `as.numeric()` coercion. If all non-NA values convert successfully → **numerical**; otherwise → **categorical**. Both types are held internally — no pre-splitting of the file is needed.

**4. Filter expression**  
The expression matrix is subsetted to only the genes present in at least one signature. If `sample_list` is provided, the matrix is further filtered to those samples.

**5. Per-variable loop**  
For each clinical variable in turn:

- Samples lacking a valid value for that variable are dropped.
- Variables with fewer than `min_samples` valid samples are skipped.
- Categorical variables with fewer than 2 groups are skipped.
- **Z-scores are recomputed on this reduced sample set** — the distribution is specific to the samples available for each variable.

**6. Activity scores**  
For each signature, per-sample activity scores are accumulated:

```
activity(sample) = Σ z(up-gene, sample) − Σ z(down-gene, sample)
```

A high score means the sample over-expresses up-regulated genes and under-expresses down-regulated genes, consistent with the signature.

**7. Statistical tests**

- **Numerical**: Pearson correlation between activity scores and the variable values. All r values are reported; only those with p < 0.05 appear in the "significant" output (others set to 0).
- **Categorical**: Kruskal-Wallis H-test. Benjamini-Hochberg FDR correction is applied across all signatures for each variable.

P-values are clamped at 1×10⁻¹⁰⁰ before −log10 transformation to avoid Inf.

**8. Heatmaps**  
Summary heatmaps are saved as PDFs (variables × signatures) clustered by Euclidean distance. Clustering is disabled when fewer than 2 rows or 2 columns are present.

---

## Usage Example

Generic call, using your own files:

```r
source("clinicalAssociation.R")

run_clinical_association(
  sig_folder    = "signatures/",
  sig_pattern   = "*.sig.txt",
  expression    = "data/expression.csv",
  clinical      = "data/clinical.csv",
  output_root   = "results/project_name",
  sample_list   = NULL,       # or "data/samples_of_interest.txt"
  min_samples   = 10,
  support_plots = TRUE,
  figure_size   = 20
)
```

### Current example (bundled with the Shiny app)

This is the exact call the app's **"Run Example Data"** button runs, using the dataset shipped in `Sig_ClinicalCorr_ShinyApp/data/` (2 signatures — `Lymphoid_Dysregulation`, `Myeloid_Dysregulation` — a 968-sample TB Discovery expression matrix, and a clinical file with one categorical variable `Endotypes` and one numerical variable `R/SI balance score`):

```r
source("R/clinicalAssociation.R")

run_clinical_association(
  sig_folder    = "data",
  sig_pattern   = "sig.*.txt",
  expression    = "data/TBDiscovery.07112021.csv",
  clinical      = "data/SampleClusters_res1.csv",
  output_root   = "results",
  min_samples   = 6,
  support_plots = TRUE,
  figure_size   = 20,
  save_temp     = TRUE
)
```

---

## Output Files

Outputs are organised into five locations. Running the current example above (`output_root = "results"`) produces exactly this set of files — note that the numerical clinical variable `R/SI balance score` becomes `R.SI.balance.score` in filenames, since `make.names()` is applied to build safe file paths:

```
results.categorical_variables_pvalues.euclidean.pdf
results.categorical_variables.minusLog10_pvalues.csv
results.numerical_var_correlations_significant_pcc.csv
results.numerical_var_significant_correlations.euclidean.pdf

pdf-numerical.results/
└── results.numerical.R.SI.balance.score.scatterplots.pdf

pdf-categorical.results/
└── results.categorical.Endotypes.boxplots.pdf

support-clinical.results/
├── results.Categorical_variable.Endotypes.activity_and_clinvar.csv
└── results.Numerical_variable.R.SI.balance.score.activity_and_clinvar.csv

temp-files.results/
├── results.expr.consistent.csv
├── results.Categorical_variable.Endotypes.activity_score_all.csv
├── results.Numerical_variable.R.SI.balance.score.activity_score_all.csv
├── results.numerical_var_correlations_all_pcc.csv
├── results.numerical_var_correlations_pcc_pvalue.csv
├── results.numerical_var_all_correlations.euclidean.pdf
├── results.categorical_variables.minusLog10_qvalues.csv
└── results.categorical_variables_qvalues.euclidean.pdf
```

With a different `output_root`, every `results` prefix above becomes `<base_name>` and the whole set relocates under `dirname(output_root)/` instead — see the general pattern below.

### Main results — `dirname(output_root)/`

Primary outputs for downstream interpretation. (In the current example, `output_root = "results"` has no directory component, so `dirname(output_root)` is just the working directory — these files land there directly, alongside the subfolders below, rather than inside a nested `results/` folder.)

| File | Content |
|---|---|
| `…significant_pcc.csv` | Pearson r × signature; non-significant values set to 0 |
| `…significant_correlations.euclidean.pdf` | Clustered heatmap, blue–white–red (−1 to 1) |
| `…minusLog10_pvalues.csv` | −log10(KW p-value) for every variable × signature pair |
| `…categorical_variables_pvalues.euclidean.pdf` | Clustered heatmap, grey (0) → white → red (≥4) |

### `pdf-numerical.<base>/`

| File | Content |
|---|---|
| `….numerical.<var>.scatterplots.pdf` | One page per signature: activity vs variable scatter (requires `support_plots = TRUE`) |

### `pdf-categorical.<base>/`

| File | Content |
| --- | --- |
| `….categorical.<var>.boxplots.pdf` | One page per signature: activity score boxplots by group (requires `support_plots = TRUE`) |

### `support-clinical.<base>/`

Per-variable activity score matrices joined with the clinical variable value — one file per variable. Useful for custom downstream plots or modelling.

| File | Content |
|---|---|
| `….Numerical_variable.<var>.activity_and_clinvar.csv` | Samples × (signature scores + variable value) |
| `….Categorical_variable.<var>.activity_and_clinvar.csv` | Samples × (signature scores + group label) |

### `temp-files.<base>/`

Intermediate outputs. Not primary results but useful for debugging, re-running partial analyses, or feeding other tools.

| File | Content |
|---|---|
| `….expr.consistent.csv` | Expression matrix filtered to signature genes |
| `….Numerical_variable.<var>.activity_score_all.csv` | Raw activity score matrix (signatures × samples) for that variable's sample subset |
| `….Categorical_variable.<var>.activity_score_all.csv` | Same for categorical variables |
| `….numerical_var_correlations_all_pcc.csv` | All Pearson r values (before significance filtering) |
| `….numerical_var_correlations_pcc_pvalue.csv` | Raw p-values corresponding to all PCC values |
| `….numerical_var_all_correlations.euclidean.pdf` | Heatmap of all Pearson r (duplicate of pdf-numerical copy) |
| `….categorical_variables.minusLog10_qvalues.csv` | −log10(BH q-value) for every variable × signature pair |
| `….categorical_variables_qvalues.euclidean.pdf` | Clustered heatmap of q-values, grey → white → red |

### Color palettes

| Analysis | Palette | Scale |
| --- | --- | --- |
| Numerical correlations | Blue → White → Red | −1 to 1 |
| Categorical p / q values | Grey → White → Red | 0 to 4 (−log10) |

Grey encodes exactly zero (p = 1; not significant). White marks the onset of any evidence of association. Red marks the highest significance (−log10(p) ≥ 4, i.e. p ≤ 0.0001).

---

## Notes

- **Gene symbols are matched case-insensitively** — both the expression matrix and signature files are uppercased on load.
- **Z-scores are recomputed per variable** using only the samples with a valid value. A variable covering 80 samples will produce different z-scores than one covering 200. This matches the original Python implementation.
- **Activity score is a raw sum, not a mean.** Signatures with more genes will produce larger absolute activity scores. Keep this in mind when comparing signatures of very different sizes.
- **The single clinical file is split internally.** You do not need to pre-separate numerical and categorical variables.
