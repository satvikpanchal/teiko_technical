setup:
	pip install -r requirements.txt

pipeline:
	python load_data.py
	python pipeline.py

dashboard:
	streamlit run dashboard.py

.PHONY: setup pipeline dashboard
