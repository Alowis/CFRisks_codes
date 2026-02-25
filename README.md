# Compound Hazards and Flood Losses in Europe — Supporting Code

Code repository associated with the paper:

> **Compounding hazards increase flood economic losses across Europe** https://www.researchsquare.com/article/rs-7075867/v1 (Under Review)

This repository contains the data processing scripts and analysis notebooks used to:

* Match hydrometeorological hazards to observed flood impacts
* Construct event-level compound hazard classifications
* Quantify economic losses by compound complexity
* Perform temporal window sensitivity analysis
* Generate spatial and statistical figures presented in the manuscript

The analysis is conducted at **NUTS3 level (1981–2020)**.

---

## Data Sources

### Flood Impact Data

* HANZE database
  (Filtered to exclude coastal floods)

### Hazard Data

Preprocessed event-level hazard tables for:

* High-flow (flood precursor events)
* Drought
* Heatwave
* Coldwave
* Windstorm (ERA5 gust-based threshold)

### Administrative Boundaries

* NUTS3 (2021 version)

All spatial data are projected to **EPSG:3035** for mapping and aggregation.

For reuse, please refer to the original data providers for access conditions and redistribution terms.

---

## Project Structure

```
├── data/                    # Input and intermediate datasets (if permitted)
├── output/                  # Generated figures and aggregated tables
├── trend_analysis.ipynb     # Notebook to march HANZE flood events with compound hazards, and reproduce the event-level analysis and figures
├── ml_analysis.ipynb        # Notebook to reproduce to the nuts-level machine learning analysis 
├── utils.py                  # Helper functions

```

### 1. Hazard–Flood Matching

* Extract non-coastal flood events from HANZE.
* Match modeled hazards to flood events within hazard-specific temporal windows.
* Evaluate optimal lag selection using detection diagnostics.

Output: matched compound hazard–flood event dataset.
