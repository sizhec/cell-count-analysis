"""Run with: streamlit run dashboard.py"""

import sqlite3
from contextlib import closing

import pandas as pd
import plotly.express as px
import streamlit as st

from analysis import (
    ANALYSIS_SQL, BASELINE_SQL, BASELINE_SUMMARIES, FREQUENCY_SQL,
    response_plot, statistical_results, subject_frequencies,
)
from load_data import DB_PATH

st.set_page_config(page_title="Loblaw Bio | Immune cells", layout="wide")
st.title("Immune cell analysis")
st.caption("Loblaw Bio · Melanoma · Miraclib")

if not DB_PATH.exists():
    st.error("Database not found. Run `make pipeline` first.")
    st.stop()


@st.cache_data
def query(sql: str, database_version: int) -> pd.DataFrame:
    with closing(sqlite3.connect(DB_PATH)) as connection:
        return pd.read_sql_query(sql, connection)


version = DB_PATH.stat().st_mtime_ns
frequencies = query(FREQUENCY_SQL, version)
comparison = query(ANALYSIS_SQL, version)
baseline = query(BASELINE_SQL, version)

overview, response_tab, baseline_tab = st.tabs([
    "Sample overview", "Response comparison", "Baseline cohort"
])

with overview:
    st.metric("Samples", f"{frequencies['sample'].nunique():,}")
    selected_sample = st.selectbox("Sample", frequencies["sample"].unique())
    sample_data = frequencies[frequencies["sample"] == selected_sample]
    left, right = st.columns([1, 2])
    left.metric("Total cells", f"{int(sample_data.total_count.iloc[0]):,}")
    left.dataframe(sample_data, hide_index=True, use_container_width=True)
    right.plotly_chart(
        px.bar(sample_data, x="population", y="percentage", text_auto=".1f",
               labels={"percentage": "Relative frequency (%)", "population": "Population"}),
        use_container_width=True,
    )
    st.download_button("Download all sample frequencies", frequencies.to_csv(index=False),
                       "sample_population_frequencies.csv", "text/csv")

with response_tab:
    st.caption("Melanoma PBMC samples from miraclib-treated subjects; yes = responder, no = non-responder.")
    times = sorted(comparison.time_from_treatment_start.unique())
    selected_times = st.multiselect("Days from treatment start", times, default=times)
    filtered = comparison[comparison.time_from_treatment_start.isin(selected_times)]
    if filtered.empty:
        st.info("Select at least one timepoint to see the comparison.")
    else:
        subjects = subject_frequencies(filtered)
        st.write(f"{filtered['sample'].nunique():,} samples from "
                 f"{subjects[['project', 'subject']].drop_duplicates().shape[0]:,} subjects")
        level = st.radio("Boxplot observations", ["Samples", "Subject means"], horizontal=True)
        plot_data = filtered if level == "Samples" else subjects
        st.plotly_chart(response_plot(plot_data, level), use_container_width=True)
        stats = statistical_results(filtered)
        st.caption("Two-sided Mann–Whitney U tests use each subject's mean across the selected visits. "
                   "Benjamini–Hochberg adjustment covers five populations; significance threshold is 0.05.")
        if stats.p_value.isna().any():
            st.info("Some comparisons cannot be tested because a response group has no valid observations.")
        significant = stats.loc[stats.significant_fdr_0_05, "population"].tolist()
        if significant:
            st.write("Significant populations: " + ", ".join(significant))
        elif stats.p_value.notna().all():
            st.write("No population meets the adjusted significance threshold.")
        st.dataframe(stats, hide_index=True, use_container_width=True)
        st.download_button("Download displayed statistics", stats.to_csv(index=False),
                           "response_statistics.csv", "text/csv")
        st.info("These are exploratory associations. Post-treatment visits cannot establish a baseline "
                "predictor, and significance does not measure predictive accuracy. "
                "Trying multiple time filters adds comparisons beyond the five-test correction.")

with baseline_tab:
    st.caption("Melanoma · PBMC · miraclib · day 0. Subjects are counted once within each project.")
    left, right = st.columns(2)
    left.metric("Baseline samples", len(baseline))
    right.metric("Baseline subjects", len(baseline[["project", "subject"]].drop_duplicates()))
    for column, (name, sql) in zip(st.columns(3), BASELINE_SUMMARIES.items()):
        counts = query(sql, version)
        category, count = counts.columns
        column.plotly_chart(
            px.bar(counts, x=category, y=count, text_auto=True,
                   title=name.removeprefix("baseline_").replace("_", " ").capitalize()),
            use_container_width=True,
        )
        column.dataframe(counts, hide_index=True, use_container_width=True)
    st.dataframe(baseline, hide_index=True, use_container_width=True)
    st.download_button("Download baseline samples", baseline.to_csv(index=False),
                       "baseline_samples.csv", "text/csv")
