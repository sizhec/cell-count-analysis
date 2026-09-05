"""Create cell_counts.db from cell-count.csv: python load_data.py"""

import sqlite3
from contextlib import closing
from pathlib import Path

import pandas as pd

ROOT = Path(__file__).resolve().parent
CSV_PATH = ROOT / "cell-count.csv"
DB_PATH = ROOT / "cell_counts.db"
POPULATIONS = ["b_cell", "cd8_t_cell", "cd4_t_cell", "nk_cell", "monocyte"]

SCHEMA = """
PRAGMA foreign_keys = ON;

CREATE TABLE projects (
    project TEXT PRIMARY KEY
);
CREATE TABLE subjects (
    project TEXT NOT NULL REFERENCES projects(project),
    subject TEXT NOT NULL,
    condition TEXT NOT NULL,
    age INTEGER NOT NULL CHECK (age >= 0),
    sex TEXT NOT NULL,
    PRIMARY KEY (project, subject)
);
CREATE TABLE samples (
    sample TEXT PRIMARY KEY,
    project TEXT NOT NULL,
    subject TEXT NOT NULL,
    treatment TEXT NOT NULL,
    response TEXT CHECK (response IN ('yes', 'no')),
    sample_type TEXT NOT NULL,
    time_from_treatment_start REAL NOT NULL,
    FOREIGN KEY (project, subject) REFERENCES subjects(project, subject)
);
CREATE TABLE cell_counts (
    sample TEXT NOT NULL REFERENCES samples(sample),
    population TEXT NOT NULL,
    count INTEGER NOT NULL CHECK (count >= 0),
    PRIMARY KEY (sample, population)
);
CREATE INDEX sample_subject ON samples(project, subject);
CREATE INDEX sample_filters ON samples(treatment, sample_type, time_from_treatment_start);
"""


def main():
    data = pd.read_csv(CSV_PATH)
    projects = data[["project"]].drop_duplicates()
    subjects = data[["project", "subject", "condition", "age", "sex"]].drop_duplicates()
    samples = data[[
        "sample", "project", "subject", "treatment", "response",
        "sample_type", "time_from_treatment_start",
    ]]
    counts = data.melt(
        id_vars="sample", value_vars=POPULATIONS,
        var_name="population", value_name="count",
    )

    DB_PATH.unlink(missing_ok=True)
    with closing(sqlite3.connect(DB_PATH)) as connection:
        connection.executescript(SCHEMA)
        projects.to_sql("projects", connection, if_exists="append", index=False)
        subjects.to_sql("subjects", connection, if_exists="append", index=False)
        samples.to_sql("samples", connection, if_exists="append", index=False)
        counts.to_sql("cell_counts", connection, if_exists="append", index=False)
    print(f"Loaded {len(data):,} samples into {DB_PATH.name}")


if __name__ == "__main__":
    main()
