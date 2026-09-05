# Immune cell analysis

Loads the supplied cell counts into SQLite and displays sample frequencies, treatment-response comparisons, and baseline summaries.

## Run

Use Python 3.12. In GitHub Codespaces:

```bash
make setup
make pipeline
make dashboard
```

Open port **8501** in the Codespaces **Ports** panel. Locally, open the [dashboard](http://localhost:8501).
Repository: https://github.com/sizhec/cell-count-analysis
Dashboard: https://friendly-cod-rv57w766rgp2pwq-8501.app.github.dev/

Without Make:

```bash
python -m pip install -r requirements.txt
python load_data.py
python analysis.py
python -m streamlit run dashboard.py
```

The pipeline rebuilds `cell_counts.db` directly and writes results to `outputs/`. It replaces the previous database; a failed load may leave a partial database. Correct the input and rerun the pipeline.

## Files

| File | Purpose |
| --- | --- |
| `load_data.py` | Creates the tables and loads the CSV using pandas. |
| `analysis.py` | Queries the data, calculates frequencies and statistics, and saves tables and plots. |
| `dashboard.py` | Displays the results with sample and timepoint controls. |
| `requirements.txt` | Lists the Python dependencies. |
| `Makefile` | Runs setup, the pipeline, and the dashboard. |
| `.devcontainer/devcontainer.json` | Configures Python 3.12 and port 8501 for Codespaces. |
| `cell-count.csv` | Original input data. |
| `cell_counts.db` | Generated SQLite database. |
| `outputs/` | Frequency tables, statistical results, baseline summaries, and two offline HTML boxplots. |

Loading, analysis, and display each have one Python file. The dashboard reuses the analysis queries and functions.

## Database

`projects → subjects → samples → cell_counts`

- **projects:** one row per project, keyed by the CSV project name.
- **subjects:** condition, age, and sex, keyed by `(project, subject)`.
- **samples:** treatment, response, sample type, and visit time, keyed by sample name.
- **cell_counts:** one row per sample and population, keyed by `(sample, population)`.

Foreign keys connect the tables. This avoids repeating subject details across visits and cell populations. Indexes support subject joins and cohort filters. SQLite is sufficient for hundreds of projects and thousands of samples; the same tables support new analytical queries. New populations require updating the loader's population list. A server database would be appropriate if many users needed to write concurrently.

## Analysis and results

Each frequency is `100 × population count / total cells in the sample`. Zero totals produce an undefined percentage instead of dividing by zero. The supplied data has 10,500 samples and 52,500 frequency rows.

Response comparisons include only **melanoma, miraclib, PBMC** samples with response **yes/no**. Repeated visits are averaged per subject before a two-sided Mann–Whitney U test. Benjamini–Hochberg adjustment covers the five population tests, with significance at adjusted p < 0.05. The dashboard offers sample-level and subject-mean boxplots.

**No population is significant after adjustment.** CD4 T cells have raw p = 0.01242 and adjusted p = 0.06211. These are exploratory associations, not validated predictors or causal effects. Cell fractions are dependent, confounding is not adjusted for, and trying additional time filters adds comparisons beyond the five-test correction.

Baseline means time from treatment start is **0**, with the same melanoma/miraclib/PBMC filter. Project totals count samples; response and sex totals count distinct subjects within each project.

| Baseline group | Count |
| --- | ---: |
| Samples / subjects | 656 / 656 |
| Project prj1 / prj3 samples | 384 / 272 |
| Responders / non-responders | 331 / 325 |
| Female / male subjects | 312 / 344 |

Successfully run in GitHub Codespaces using `make setup`, `make pipeline`, and `make dashboard`.
