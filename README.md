# 🧬 RNA-seq Differential Expression Analysis

**COVID-19 vs Healthy (Exploratory Pipeline)**

---

## 📌 What this analysis does

This script performs an **exploratory differential expression analysis** comparing
**COVID-19 samples vs Healthy samples** using raw RNA-seq count data.

It:

* reads metadata and count files
* matches samples between files
* filters low-expression genes
* computes log2 fold changes
* performs gene-wise statistical testing
* corrects for multiple testing
* extracts significant DEGs
* prepares outputs for downstream analysis (DAVID, plots)

⚠️ This pipeline is **educational / exploratory**.
For final RNA-seq inference, **DESeq2 is recommended**.

---

## 📁 Input Files

### 1️⃣ Metadata (GEO series matrix)

* Contains sample information and disease state
* Used to map samples → groups

### 2️⃣ Raw counts file

* Rows = genes (Ensembl IDs)
* Columns = samples
* Values = raw read counts (not normalized)

---

## 🧱 Objects Created & What They Store

| Object            | Type             | Meaning                         |
| ----------------- | ---------------- | ------------------------------- |
| `metadata`        | data.frame       | sample ↔ disease state          |
| `counts`          | data.frame       | raw gene counts                 |
| `covid_samples`   | character vector | column names of COVID samples   |
| `healthy_samples` | character vector | column names of Healthy samples |
| `counts_filt`     | data.frame       | filtered counts                 |
| `log2fc`          | numeric vector   | fold change per gene            |
| `pvals`           | numeric vector   | raw p-values                    |
| `results`         | data.frame       | combined gene-level stats       |
| `sig`             | data.frame       | final DEG list                  |

---

## 🧭 Step-by-Step Pipeline Flow

### STEP 1 — Load Metadata

**Purpose:** know which sample belongs to which group

* Extract sample IDs
* Extract disease state (COVID / Healthy)
* Store in a clean metadata table

---

### STEP 2 — Load Raw Counts

**Purpose:** load expression data

* Genes as rows
* Samples as columns
* Raw counts (no normalization yet)

---

### STEP 3 — Sanity Checks

**Purpose:** catch errors early

* Number of genes
* Number of samples
* Sample name preview

---

### STEP 4 — Match Metadata to Counts

**Purpose:** align sample names across files

* Count column names may not exactly match metadata IDs
* Fuzzy string matching is used
* Samples are assigned to:

  * `covid_samples`
  * `healthy_samples`

⚠️ This step is critical — wrong matching = wrong biology.

---

### STEP 5 — Filter Low-Expression Genes

**Purpose:** remove noise before DE analysis

* Compute mean expression per gene
* Keep genes with `mean > 10`

Why?

* Low-count genes inflate noise
* Improves statistical stability
* Does **not** bias results (filtering is condition-independent)

---

### STEP 6 — Calculate Log2 Fold Change

**Purpose:** measure direction + magnitude of change

For each gene:

* Compute mean COVID expression
* Compute mean Healthy expression
* Add pseudocount (`+1`)
* Calculate:

```
log2FC = log2(COVID / Healthy)
```

Interpretation:

* `log2FC > 0` → upregulated in COVID
* `log2FC < 0` → downregulated in COVID

---

### STEP 7 — Statistical Testing (t-test)

**Purpose:** test if group differences are statistically meaningful

For each gene:

* Extract COVID counts
* Extract Healthy counts
* Perform two-sample t-test
* Store one p-value per gene

⚠️ Note:

* t-tests are **not ideal** for RNA-seq
* Used here for intuition / learning
* DESeq2 should be used for final results

---

### STEP 8 — Combine Results

**Purpose:** one row = one gene

Create a results table containing:

* gene ID
* log2 fold change
* raw p-value

---

### STEP 9 — Multiple Testing Correction

**Purpose:** control false discoveries

* Apply Benjamini–Hochberg (BH) correction
* Create adjusted p-values (`padj`)
* Controls **False Discovery Rate (FDR)**

Never trust raw p-values alone in genomics.

---

### STEP 10 — Select Significant DEGs

**Purpose:** define biologically meaningful genes

Criteria:

* `|log2FC| > 1`  (≥ 2× change)
* `padj < 0.05`   (FDR controlled)

Then:

* Sort genes by absolute log2FC (strongest effects first)

Result:
➡️ final DEG table (`sig`)

---

### STEP 11 — Outputs

Generated files include:

* DEG lists for enrichment tools
* Background gene list
* Volcano plot
* Full results table

These are ready for:

* DAVID
* GO / KEGG analysis
* Reporting & visualization

---

## 🧠 Key Concepts to Remember (When You’re Rusty)

* **Counts are not normalized manually** → tools handle this
* **Filtering removes noise, not biology**
* **log2FC = effect size**
* **padj < 0.05 = significance in genomics**
* **One p-value per gene**
* **Metadata mistakes are more dangerous than code bugs**

---

## 🚨 Limitations of This Script

* Uses t-tests instead of count-based models
* Assumes approximate normality
* Not suitable for publication-grade DE analysis

✔️ Use this to **understand the workflow**
✔️ Use **DESeq2** for final inference

---

## 🧬 One-Line Mental Flow (Save This)

```
Metadata → Counts → Match samples → Filter noise → log2FC → p-values → FDR → DEGs
```

---

## 💙 Note to Future Me

If this script looks confusing again:

* Start by checking `metadata`
* Then check `covid_samples` / `healthy_samples`
* Then re-read the **Pipeline Flow** section

You already figured this out once — you can do it again 😌

