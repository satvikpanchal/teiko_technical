# Loblaw Bio — Immune Cell Population Analysis

A small but complete bioinformatics pipeline + dashboard for Bob Loblaw's
miraclib clinical trial. It loads `cell-count.csv` into a normalized SQLite
database, computes per-sample cell-population frequencies, tests whether any
population separates responders from non-responders, summarizes a baseline
subset, and serves everything through an interactive Streamlit dashboard.

---

## 1. Setup & Running

The grader runs this in **GitHub Codespaces**. Three `make` targets do
everything:

```bash
# 1. Install dependencies (pandas, numpy, scipy, matplotlib, seaborn, streamlit)
make setup

# 2. Run the full pipeline: build the DB, load data, generate all tables + plots
make pipeline

# 3. Launch the interactive dashboard
make dashboard
```

`make pipeline` runs `load_data.py` then `pipeline.py` in sequence with no
manual intervention. It creates `teiko.db` in the repo root and writes all
results to `outputs/`.

Running scripts individually works too:

```bash
python load_data.py     # creates teiko.db (idempotent — safe to re-run)
python pipeline.py       # writes outputs/ (creates the directory if missing)
streamlit run dashboard.py
```

> **Note on the input file:** `cell-count.csv` must be in the repository root
> (it is referenced by a relative path). All scripts resolve paths relative to
> their own location, so no absolute paths are hardcoded.

### Outputs produced by `make pipeline`

| File | Description |
|------|-------------|
| `teiko.db` | SQLite database (4 normalized tables) |
| `outputs/tables/cell_frequencies.csv` | Per-sample × population relative frequencies (Part 2) |
| `outputs/tables/statistical_results.csv` | Mann-Whitney U results per population (Part 3) |
| `outputs/plots/boxplot_cell_populations.png` | Responder vs non-responder boxplots (Part 3) |
| `outputs/tables/subset_analysis.csv` | Baseline subset summary (Part 4) |

---

## 2. Database Schema

The flat CSV repeats every patient's demographics on each of its samples. The
schema normalizes that into four tables following the natural biological
hierarchy **project → subject → sample → measurement**:

```
projects (project_id PK, project_name)
   │  1─to─many
subjects (subject_id PK, project_id FK, condition, age, sex, treatment, response)
   │  1─to─many
samples  (sample_id PK, subject_id FK, sample_type, time_from_treatment_start)
   │  1─to─one (per assay)
cell_counts (id PK, sample_id FK, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte)
```

**Table rationale**

- **`projects`** — one row per project; lets project metadata grow (sponsor,
  site, dates) without touching subject rows.
- **`subjects`** — demographics, treatment, and response live here because they
  are *constant for a subject across all its samples* (verified in the data).
  `response` is **nullable**: healthy / untreated subjects have no response
  label, modeled as `NULL` rather than an empty string.
- **`samples`** — the unit of measurement. Holds the per-sample attributes
  `sample_type` (PBMC vs WB) and `time_from_treatment_start`.
- **`cell_counts`** — the five immune-cell measurements, split from `samples`
  so additional assays/panels can be added later without widening one giant
  table.

This removes update anomalies (a subject's sex is stored once, not 3×),
enforces referential integrity via foreign keys, and keeps each fact in exactly
one place.

**Indexing strategy** (created in `load_data.py`):

- `idx_subjects_project` and `idx_subjects_filters (condition, treatment, response)`
  accelerate the cohort filters used throughout the analysis.
- `idx_samples_type_time (sample_type, time_from_treatment_start)` speeds the
  baseline-subset queries.
- `idx_samples_subject` and `idx_cell_counts_sample` make the
  subject→sample→counts joins index-driven.

**Scaling to hundreds of projects / thousands of samples / varied analytics**

- **Indexing** — the composite indexes above keep cohort filtering and joins
  fast as row counts grow; add covering indexes for the hottest queries.
- **`assay_type` table** — to support more than the current 5-population panel,
  pivot `cell_counts` into a long `measurements (sample_id, assay_type_id,
  population, count)` table with an `assay_types` lookup. New panels (cytokines,
  flow markers) become new rows, not schema changes.
- **Partitioning** — at very large scale, partition or shard by `project_id`
  (the natural isolation boundary), and consider moving from SQLite to
  Postgres/DuckDB. The schema ports unchanged.
- **Materialized frequency table** — `cell_frequencies` could be persisted as a
  table (or materialized view) and refreshed on load, so the dashboard reads
  pre-aggregated rows instead of recomputing.
- **Provenance columns** — `created_at` / batch IDs on `samples` enable
  time-based partitioning and reproducibility as data arrives incrementally.

---

## 3. Code Structure

Everything lives in the repo root so the scripts run with no module path setup.

| File | Purpose | Why structured this way |
|------|---------|-------------------------|
| `load_data.py` | Defines the schema and loads `cell-count.csv` into `teiko.db`. | Idempotent (drops + recreates tables each run). Parents (`projects`, `subjects`) are de-duplicated in memory, then all four tables are bulk-inserted with `executemany`. Empty cells become `NULL`. Pure stdlib (`csv` + `sqlite3`) so loading has no third-party dependency. |
| `pipeline.py` | Runs Parts 2–4 and writes every CSV/PNG to `outputs/`. | One function per part for readability; reads from the DB (Part 4 uses pure SQL filtering as required) and uses pandas/scipy/seaborn for analysis and plotting. Creates `outputs/` automatically and uses the headless matplotlib backend so it works in Codespaces/CI. |
| `dashboard.py` | Streamlit dashboard (4 tabs + sidebar filters). | Reads pre-computed `outputs/` for the static views and queries `teiko.db` **live** for the raw-data and subset-explorer tabs. Cached data loaders (`@st.cache_data`) keep it responsive. |
| `requirements.txt` | Pinned minimum versions of the analysis stack. | |
| `Makefile` | `setup`, `pipeline`, `dashboard` targets. | Single entry point for the grader. |
| `outputs/` | Generated tables and plots. | Created by `pipeline.py`. |

**Analysis decisions**

- The responder-vs-non-responder comparison is restricted to
  **melanoma + miraclib + PBMC** samples with a non-null response, exactly as
  specified. Each population is tested with a two-sided **Mann-Whitney U** test
  (non-parametric — appropriate for frequency data without a normality
  assumption); significance is flagged at *p* < 0.05 with stars on the boxplot.
- The baseline subset is **melanoma + miraclib + PBMC + time = 0**. The data
  contains both PBMC and WB sample types, so PBMC is filtered explicitly rather
  than assumed.

---

## 4. Dashboard

Run locally with `make dashboard` (or `streamlit run dashboard.py`), then open
the URL Streamlit prints (default `http://localhost:8501`). It can also be
deployed for free via **Streamlit Community Cloud** by pointing it at this repo
and `dashboard.py`.

The dashboard has four tabs:

1. **Data Overview** — summary metrics, the frequency table, and a mean-%
   bar chart per cell type.
2. **Statistical Analysis** — the results table with significant rows
   highlighted in green, the boxplot, and a plain-language summary.
3. **Subset Explorer** — the pre-computed baseline subset plus a **live**
   condition/treatment/timepoint query against the DB.
4. **Raw Data** — searchable full table from the DB with a CSV download button.

Sidebar multiselects (condition, treatment, sample type, response) filter the
relevant tabs; leaving a filter empty means "All".

---

## 5. Key Findings

Among melanoma patients on miraclib (PBMC samples), **CD4 T cells are the only
population with a statistically significant difference between responders and
non-responders** (Mann-Whitney U, *p* ≈ 0.013): responders show a higher mean
CD4 T-cell relative frequency. The other four populations (B cell, CD8 T cell,
NK cell, monocyte) do not reach significance at *p* < 0.05, though B cells trend
lower in responders (*p* ≈ 0.056). This suggests baseline/over-trial CD4 T-cell
frequency is the most promising candidate biomarker for predicting miraclib
response and warrants follow-up in a larger cohort.
