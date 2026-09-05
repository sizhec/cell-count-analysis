PYTHON ?= python

.PHONY: setup pipeline dashboard

setup:
	$(PYTHON) -m pip install -r requirements.txt

pipeline:
	$(PYTHON) load_data.py
	$(PYTHON) analysis.py

dashboard:
	$(PYTHON) -m streamlit run dashboard.py --server.address 0.0.0.0 --server.port 8501 --server.headless true --browser.gatherUsageStats false
