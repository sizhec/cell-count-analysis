"""Generate all required analysis tables and plots from cell_counts.db."""

import sqlite3
from contextlib import closing

import pandas as pd
import plotly.express as px
from scipy.stats import false_discovery_control, mannwhitneyu

from load_data import DB_PATH, POPULATIONS, ROOT

OUTPUT_DIR = ROOT / "outputs"

FREQUENCY_SQL = """
SELECT sample, total_count, population, count,
       100.0 * count / NULLIF(total_count, 0) AS percentage
FROM (
    SELECT *, SUM(count) OVER (PARTITION BY sample) AS total_count
    FROM cell_counts
)
ORDER BY sample, population
"""

ANALYSIS_SQL = f"""
SELECT s.project, s.subject, s.sample, s.time_from_treatment_start,
       s.response, f.population, f.percentage
FROM samples s
JOIN subjects u USING (project, subject)
JOIN ({FREQUENCY_SQL}) f USING (sample)
WHERE u.condition = 'melanoma'
  AND s.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND s.response IN ('yes', 'no')
"""

BASELINE_SQL = """
SELECT s.project, s.subject, s.sample, s.response, u.sex
FROM samples s
JOIN subjects u USING (project, subject)
WHERE u.condition = 'melanoma'
  AND s.treatment = 'miraclib'
  AND s.sample_type = 'PBMC'
  AND s.time_from_treatment_start = 0
ORDER BY s.project, s.subject, s.sample
"""


BASELINE_SUMMARIES = {
    "baseline_samples_by_project": f"""
        SELECT project, COUNT(*) AS sample_count
        FROM ({BASELINE_SQL}) GROUP BY project ORDER BY project
    """,
    "baseline_subjects_by_response": f"""
        SELECT response, COUNT(*) AS subject_count
        FROM (SELECT DISTINCT project, subject, response FROM ({BASELINE_SQL}))
        GROUP BY response ORDER BY response
    """,
    "baseline_subjects_by_sex": f"""
        SELECT sex, COUNT(*) AS subject_count
        FROM (SELECT DISTINCT project, subject, sex FROM ({BASELINE_SQL}))
        GROUP BY sex ORDER BY sex
    """,
}


def subject_frequencies(data):
    if data.groupby(["project", "subject"]).response.nunique().gt(1).any():
        raise ValueError("A subject cannot belong to both response groups")
    return data.groupby(
        ["project", "subject", "response", "population"], as_index=False
    )["percentage"].mean()


def statistical_results(data):
    # Each subject contributes one observation, even with repeated visits.
    subject_data = subject_frequencies(data)
    rows = []
    for population in POPULATIONS:
        group = subject_data[subject_data.population == population]
        responders = group.loc[group.response == "yes", "percentage"].dropna()
        nonresponders = group.loc[group.response == "no", "percentage"].dropna()
        statistic, p_value = float("nan"), float("nan")
        if len(responders) and len(nonresponders):
            statistic, p_value = mannwhitneyu(
                responders, nonresponders, alternative="two-sided", method="asymptotic"
            )
        rows.append({
            "population": population,
            "n_responder_subjects": len(responders),
            "n_nonresponder_subjects": len(nonresponders),
            "responder_median_percentage": responders.median(),
            "nonresponder_median_percentage": nonresponders.median(),
            "median_difference_percentage_points": responders.median() - nonresponders.median(),
            "mann_whitney_u": statistic,
            "p_value": p_value,
        })
    results = pd.DataFrame(rows)
    adjusted = false_discovery_control(results.p_value.fillna(1), method="bh")
    results["p_value_fdr_bh"] = pd.Series(adjusted).where(results.p_value.notna())
    results["significant_fdr_0_05"] = results["p_value_fdr_bh"] < 0.05
    return results


def response_plot(data, title):
    figure = px.box(
        data, x="population", y="percentage", color="response", points="outliers",
        labels={"percentage": "Relative frequency (%)", "population": "Cell population", "response": "Response"},
        title=title,
        category_orders={"response": ["yes", "no"], "population": list(POPULATIONS)},
        color_discrete_map={"yes": "#007f86", "no": "#ca5934"},
    )
    figure.update_layout(template="plotly_white", boxmode="group")
    return figure


def main():
    if not DB_PATH.exists():
        raise FileNotFoundError("Run python load_data.py before analysis.py")
    OUTPUT_DIR.mkdir(exist_ok=True)
    with closing(sqlite3.connect(DB_PATH)) as connection:
        frequencies = pd.read_sql_query(FREQUENCY_SQL, connection)
        comparison = pd.read_sql_query(ANALYSIS_SQL, connection)
        baseline = pd.read_sql_query(BASELINE_SQL, connection)
        for name, sql in BASELINE_SUMMARIES.items():
            pd.read_sql_query(sql, connection).to_csv(OUTPUT_DIR / f"{name}.csv", index=False)

    frequencies.to_csv(OUTPUT_DIR / "sample_population_frequencies.csv", index=False)
    comparison.to_csv(OUTPUT_DIR / "miraclib_response_frequencies.csv", index=False)
    statistics = statistical_results(comparison)
    statistics.to_csv(OUTPUT_DIR / "miraclib_response_statistics.csv", index=False)
    baseline.to_csv(OUTPUT_DIR / "baseline_melanoma_miraclib_pbmc_samples.csv", index=False)
    subject_data = subject_frequencies(comparison)
    subject_data.to_csv(OUTPUT_DIR / "miraclib_subject_frequencies.csv", index=False)
    response_plot(comparison, "Melanoma PBMC / miraclib: sample frequencies").write_html(
        OUTPUT_DIR / "miraclib_response_boxplot.html", include_plotlyjs=True
    )
    response_plot(subject_data, "Melanoma PBMC / miraclib: subject mean frequencies").write_html(
        OUTPUT_DIR / "miraclib_subject_boxplot.html", include_plotlyjs=True
    )

    print(statistics[["population", "p_value_fdr_bh", "significant_fdr_0_05"]].to_string(index=False))
    print(f"Generated analysis outputs in {OUTPUT_DIR.name}/")


if __name__ == "__main__":
    main()
