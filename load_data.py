"""
load_data.py
============
Initializes the SQLite database (``teiko.db``) with a normalized relational
schema and loads every row from ``cell-count.csv``.

Run from the repository root with::

    python load_data.py

The script is idempotent: it drops and recreates all tables on each run, so it
is always safe to re-run. All paths are relative to the directory containing
this file, so it works from GitHub Codespaces or any checkout location.
"""

import csv
import os
import sqlite3

# Resolve paths relative to this file so the script works from any CWD.
ROOT = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(ROOT, "cell-count.csv")
DB_PATH = os.path.join(ROOT, "teiko.db")

# ---------------------------------------------------------------------------
# Schema
# ---------------------------------------------------------------------------
# The data is normalized into four tables to remove the redundancy of the flat
# CSV (where every subject's demographics are repeated on each of its samples):
#
#   projects      one row per project
#   subjects      one row per patient; demographics + treatment + response live
#                 here because they are constant for a subject across samples
#   samples       one row per biological sample (the unit of measurement)
#   cell_counts   the five immune-cell measurements for each sample
#
# Foreign keys wire the hierarchy projects -> subjects -> samples -> cell_counts.
SCHEMA = """
PRAGMA foreign_keys = ON;

DROP TABLE IF EXISTS cell_counts;
DROP TABLE IF EXISTS samples;
DROP TABLE IF EXISTS subjects;
DROP TABLE IF EXISTS projects;

CREATE TABLE projects (
    project_id   TEXT PRIMARY KEY,
    project_name TEXT NOT NULL
);

CREATE TABLE subjects (
    subject_id TEXT PRIMARY KEY,
    project_id TEXT NOT NULL,
    condition  TEXT,
    age        INTEGER,
    sex        TEXT,
    treatment  TEXT,
    response   TEXT,                 -- nullable: healthy subjects have no response
    FOREIGN KEY (project_id) REFERENCES projects (project_id)
);

CREATE TABLE samples (
    sample_id                 TEXT PRIMARY KEY,
    subject_id                TEXT NOT NULL,
    sample_type               TEXT,
    time_from_treatment_start INTEGER,
    FOREIGN KEY (subject_id) REFERENCES subjects (subject_id)
);

CREATE TABLE cell_counts (
    id         INTEGER PRIMARY KEY AUTOINCREMENT,
    sample_id  TEXT NOT NULL,
    b_cell     INTEGER,
    cd8_t_cell INTEGER,
    cd4_t_cell INTEGER,
    nk_cell    INTEGER,
    monocyte   INTEGER,
    FOREIGN KEY (sample_id) REFERENCES samples (sample_id)
);

-- Indexes for the common analytical access patterns (filtering by subject
-- attributes, by sample type / timepoint, and joining counts back to samples).
CREATE INDEX idx_subjects_project   ON subjects (project_id);
CREATE INDEX idx_subjects_filters   ON subjects (condition, treatment, response);
CREATE INDEX idx_samples_subject    ON samples (subject_id);
CREATE INDEX idx_samples_type_time  ON samples (sample_type, time_from_treatment_start);
CREATE INDEX idx_cell_counts_sample ON cell_counts (sample_id);
"""


def _to_int(value):
    """Convert a CSV cell to int, treating empty strings as NULL."""
    value = (value or "").strip()
    return int(value) if value != "" else None


def _to_text(value):
    """Normalize a CSV cell to text, treating empty strings as NULL."""
    value = (value or "").strip()
    return value if value != "" else None


def main():
    if not os.path.exists(CSV_PATH):
        raise FileNotFoundError(f"Could not find {CSV_PATH}")

    conn = sqlite3.connect(DB_PATH)
    conn.execute("PRAGMA foreign_keys = ON;")
    try:
        # Recreate the schema from scratch (idempotent).
        conn.executescript(SCHEMA)

        projects = {}          # project_id -> name
        subjects = {}          # subject_id -> row tuple
        samples = []           # list of sample row tuples
        cell_counts = []       # list of count row tuples

        with open(CSV_PATH, newline="") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                pid = row["project"].strip()
                sid = row["subject"].strip()
                sample = row["sample"].strip()

                # Deduplicate parents on the fly.
                projects.setdefault(pid, pid)
                subjects.setdefault(
                    sid,
                    (
                        sid,
                        pid,
                        _to_text(row["condition"]),
                        _to_int(row["age"]),
                        _to_text(row["sex"]),
                        _to_text(row["treatment"]),
                        _to_text(row["response"]),
                    ),
                )

                samples.append(
                    (
                        sample,
                        sid,
                        _to_text(row["sample_type"]),
                        _to_int(row["time_from_treatment_start"]),
                    )
                )
                cell_counts.append(
                    (
                        sample,
                        _to_int(row["b_cell"]),
                        _to_int(row["cd8_t_cell"]),
                        _to_int(row["cd4_t_cell"]),
                        _to_int(row["nk_cell"]),
                        _to_int(row["monocyte"]),
                    )
                )

        conn.executemany(
            "INSERT INTO projects (project_id, project_name) VALUES (?, ?)",
            [(pid, name) for pid, name in projects.items()],
        )
        conn.executemany(
            "INSERT INTO subjects "
            "(subject_id, project_id, condition, age, sex, treatment, response) "
            "VALUES (?, ?, ?, ?, ?, ?, ?)",
            list(subjects.values()),
        )
        conn.executemany(
            "INSERT INTO samples "
            "(sample_id, subject_id, sample_type, time_from_treatment_start) "
            "VALUES (?, ?, ?, ?)",
            samples,
        )
        conn.executemany(
            "INSERT INTO cell_counts "
            "(sample_id, b_cell, cd8_t_cell, cd4_t_cell, nk_cell, monocyte) "
            "VALUES (?, ?, ?, ?, ?, ?)",
            cell_counts,
        )
        conn.commit()

        print("Database loaded successfully -> teiko.db")
        print(
            f"  projects: {len(projects)}  |  subjects: {len(subjects)}  |  "
            f"samples: {len(samples)}  |  cell_count rows: {len(cell_counts)}"
        )
    finally:
        conn.close()


if __name__ == "__main__":
    main()
