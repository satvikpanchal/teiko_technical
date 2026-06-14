"""
pipeline.py
===========
Runs the full analysis against ``teiko.db`` and writes every output table and
plot into the ``outputs/`` directory (created automatically).

Run from the repository root with::

    python pipeline.py

Parts
-----
1. Per-sample cell-population frequency table          -> outputs/cell_frequencies.csv
2. Responder vs non-responder statistics + boxplot     -> outputs/statistical_results.csv
                                                          outputs/boxplot_cell_populations.png
3. Baseline (time=0) melanoma / miraclib subset summary -> outputs/subset_analysis.csv

All paths are relative to this file so the script is location-independent.
"""

import os
import sqlite3

import matplotlib

matplotlib.use("Agg")  # headless backend so it works in Codespaces / CI

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from scipy.stats import mannwhitneyu

ROOT = os.path.dirname(os.path.abspath(__file__))
DB_PATH = os.path.join(ROOT, "teiko.db")
OUT_DIR = os.path.join(ROOT, "outputs")
TABLES_DIR = os.path.join(OUT_DIR, "tables")
PLOTS_DIR = os.path.join(OUT_DIR, "plots")

POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]
PRETTY = {
    "b_cell": "B cell",
    "cd8_t_cell": "CD8 T cell",
    "cd4_t_cell": "CD4 T cell",
    "nk_cell": "NK cell",
    "monocyte": "Monocyte",
}


# ---------------------------------------------------------------------------
# Part 1 - Frequency table
# ---------------------------------------------------------------------------
def part1_frequencies(conn):
    """Long-format relative-frequency table: one row per sample x population."""
    counts = pd.read_sql_query(
        "SELECT sample_id AS sample, b_cell, cd8_t_cell, cd4_t_cell, "
        "nk_cell, monocyte FROM cell_counts", conn
    )
    counts["total_count"] = counts[POPULATIONS].sum(axis=1)

    long = counts.melt(
        id_vars=["sample", "total_count"],
        value_vars=POPULATIONS,
        var_name="population",
        value_name="count",
    )
    long["percentage"] = (long["count"] / long["total_count"] * 100).round(4)
    long = long[["sample", "total_count", "population", "count", "percentage"]]
    long = long.sort_values(["sample", "population"]).reset_index(drop=True)

    path = os.path.join(TABLES_DIR, "cell_frequencies.csv")
    long.to_csv(path, index=False)
    print(f"[Part 1] Frequency table written: {path}  ({len(long)} rows)")
    return long


# ---------------------------------------------------------------------------
# Part 2 - Statistics: responders vs non-responders
# ---------------------------------------------------------------------------
def _melanoma_miraclib_pbmc_frequencies(conn, freq):
    """Join the frequency table to metadata, restricted to the cohort of
    interest: melanoma + miraclib + PBMC samples that have a response label."""
    meta = pd.read_sql_query(
        """
        SELECT s.sample_id AS sample, subj.response AS response
        FROM samples s
        JOIN subjects subj ON s.subject_id = subj.subject_id
        WHERE subj.condition = 'melanoma'
          AND subj.treatment = 'miraclib'
          AND s.sample_type = 'PBMC'
          AND subj.response IS NOT NULL
        """,
        conn,
    )
    return freq.merge(meta, on="sample", how="inner")


def part2_statistics(conn, freq):
    cohort = _melanoma_miraclib_pbmc_frequencies(conn, freq)

    results = []
    for pop in POPULATIONS:
        sub = cohort[cohort["population"] == pop]
        responders = sub.loc[sub["response"] == "yes", "percentage"]
        non_responders = sub.loc[sub["response"] == "no", "percentage"]

        stat, p = mannwhitneyu(
            responders, non_responders, alternative="two-sided"
        )
        results.append(
            {
                "population": pop,
                "mean_responders": round(responders.mean(), 4),
                "mean_non_responders": round(non_responders.mean(), 4),
                "p_value": p,
                "significant": bool(p < 0.05),
            }
        )

    res_df = pd.DataFrame(results)
    path = os.path.join(TABLES_DIR, "statistical_results.csv")
    res_df.to_csv(path, index=False)
    print(f"[Part 2] Statistical results written: {path}")
    print(res_df.to_string(index=False))

    _boxplot(cohort, res_df)
    return res_df, cohort


def _boxplot(cohort, res_df):
    """Five side-by-side boxplots (one per population) with significance stars."""
    sns.set_theme(style="whitegrid")
    pvals = dict(zip(res_df["population"], res_df["p_value"]))

    fig, axes = plt.subplots(1, len(POPULATIONS), figsize=(20, 6), sharey=False)
    palette = {"Responder": "#2ca58d", "Non-Responder": "#d1495b"}

    for ax, pop in zip(axes, POPULATIONS):
        sub = cohort[cohort["population"] == pop].copy()
        sub["Response"] = sub["response"].map(
            {"yes": "Responder", "no": "Non-Responder"}
        )
        sns.boxplot(
            data=sub,
            x="Response",
            y="percentage",
            hue="Response",
            order=["Responder", "Non-Responder"],
            palette=palette,
            legend=False,
            ax=ax,
            showfliers=False,
        )
        ax.set_title(PRETTY[pop])
        ax.set_xlabel("")
        ax.set_ylabel("Relative frequency (%)")

        # Significance star above the boxes when p < 0.05.
        if pvals[pop] < 0.05:
            y = sub["percentage"].max()
            pad = y * 0.06 if y else 1
            star = "***" if pvals[pop] < 0.001 else "**" if pvals[pop] < 0.01 else "*"
            ax.plot([0, 0, 1, 1], [y + pad, y + 2 * pad, y + 2 * pad, y + pad],
                    lw=1.3, c="black")
            ax.text(0.5, y + 2 * pad, f"{star}\np={pvals[pop]:.1e}",
                    ha="center", va="bottom", fontsize=9)
            ax.set_ylim(top=y + 6 * pad)

    fig.suptitle(
        "Cell Population Frequencies: Responders vs Non-Responders "
        "(Melanoma, Miraclib)",
        fontsize=15,
        fontweight="bold",
    )
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    path = os.path.join(PLOTS_DIR, "boxplot_cell_populations.png")
    fig.savefig(path, dpi=150)
    plt.close(fig)
    print(f"[Part 2] Boxplot written: {path}")


# ---------------------------------------------------------------------------
# Part 3 - Baseline subset analysis
# ---------------------------------------------------------------------------
def part3_subset(conn):
    """Melanoma PBMC samples at baseline (time=0) treated with miraclib."""
    base = pd.read_sql_query(
        """
        SELECT s.sample_id, s.subject_id, subj.project_id,
               subj.response, subj.sex,
               cc.b_cell
        FROM samples s
        JOIN subjects subj   ON s.subject_id = subj.subject_id
        JOIN cell_counts cc  ON cc.sample_id  = s.sample_id
        WHERE subj.condition = 'melanoma'
          AND subj.treatment = 'miraclib'
          AND s.sample_type  = 'PBMC'
          AND s.time_from_treatment_start = 0
        """,
        conn,
    )

    # 1. Samples per project.
    per_project = (
        base.groupby("project_id")["sample_id"].count()
        .rename("sample_count").reset_index()
    )

    # Subject-level counts (deduplicate to one row per subject).
    subj = base.drop_duplicates("subject_id")

    # 2. Responder / non-responder subject counts.
    resp_counts = subj["response"].value_counts().to_dict()
    # 3. Male / female subject counts.
    sex_counts = subj["sex"].value_counts().to_dict()

    # 4. Melanoma males: average b_cell for responders at time=0 (XX.XX).
    male_resp = base[(base["sex"] == "M") & (base["response"] == "yes")]
    avg_b_male_resp = male_resp["b_cell"].mean()
    avg_b_fmt = f"{avg_b_male_resp:.2f}" if pd.notna(avg_b_male_resp) else "NA"

    # Assemble a tidy long-format summary CSV.
    # Values are stored as strings so integer counts stay clean ("384") while
    # the average keeps its two-decimal formatting ("10401.28") in one column.
    rows = []
    for _, r in per_project.iterrows():
        rows.append({"metric": "samples_per_project",
                     "category": r["project_id"], "value": str(int(r["sample_count"]))})
    for k in ("yes", "no"):
        rows.append({"metric": "subjects_by_response",
                     "category": k, "value": str(int(resp_counts.get(k, 0)))})
    for k in ("M", "F"):
        rows.append({"metric": "subjects_by_sex",
                     "category": k, "value": str(int(sex_counts.get(k, 0)))})
    rows.append({"metric": "avg_b_cell_melanoma_male_responders_t0",
                 "category": "value", "value": avg_b_fmt})

    summary = pd.DataFrame(rows, columns=["metric", "category", "value"])
    path = os.path.join(TABLES_DIR, "subset_analysis.csv")
    summary.to_csv(path, index=False)

    # Print to stdout.
    print("\n[Part 3] Baseline subset: melanoma + miraclib + PBMC + time=0")
    print(f"  Total samples in subset: {len(base)}  |  subjects: {len(subj)}")
    print("  Samples per project:")
    for _, r in per_project.iterrows():
        print(f"    {r['project_id']}: {int(r['sample_count'])}")
    print(f"  Subjects by response  -> responders (yes): "
          f"{resp_counts.get('yes', 0)}, non-responders (no): "
          f"{resp_counts.get('no', 0)}")
    print(f"  Subjects by sex       -> male (M): {sex_counts.get('M', 0)}, "
          f"female (F): {sex_counts.get('F', 0)}")
    print(f"  Avg B cells, melanoma males, responders @ time=0: {avg_b_fmt}")
    print(f"[Part 3] Subset analysis written: {path}")
    return summary


def main():
    os.makedirs(TABLES_DIR, exist_ok=True)
    os.makedirs(PLOTS_DIR, exist_ok=True)
    if not os.path.exists(DB_PATH):
        raise FileNotFoundError(
            f"{DB_PATH} not found. Run `python load_data.py` first."
        )

    conn = sqlite3.connect(DB_PATH)
    try:
        freq = part1_frequencies(conn)
        part2_statistics(conn, freq)
        part3_subset(conn)
    finally:
        conn.close()

    print("\nPipeline complete. All outputs in ./outputs/")


if __name__ == "__main__":
    main()
