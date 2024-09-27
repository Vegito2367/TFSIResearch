# TFSI Crystallography Analysis

This repository contains Python code for analyzing X-ray crystallographic data of bis(trifluoromethanesulfonyl)imide (TFSI) compounds. View at:

https://vegito2367.github.io/TFSIResearch/

## Overview

Chemical data analysis is an emerging field with applications for a variety of academic and industrial pursuits. With quantum computing and artificial intelligence becoming more mainstream, it is of great interest for not only the data scientist, but the chemist as well to take advantage of these technologies to streamline data processing, analyze large data sets, and reveal new chemical insights otherwise hidden by the depth of available data. In this project, the incorporation of modern computing and fundamental chemical analysis reveal new insights into the structural behavior of the weakly coordinating anion bis(trifluoromethylsulfonyl)amide, otherwise known as TFSI. Taking a data science approach, published solid state crystal structures available in the Cambridge Crystallographic Database (CCDC) bearing the TFSI motif of interest were categorically and statistically analyzed to determine the structural characteristics attributable to TFSI. A new software program was developed using Python to parse Crystallographic Information Files (CIF) into statistically relevant information that can be compared across structural data. This research endeavor aims to highlight the applications of data science to an otherwise foreign area of research and outline the capabilities of data processing for future structural and chemical investigations.

## Features

- Import and parse CIF (Crystallographic Information File) data
- Analyze TFSI bond lengths and angles
- Visualize TFSI molecular structure
- Perform statistical analysis on crystallographic parameters

## Requirements

- Python 3.7+
- NumPy
- SciPy
- Matplotlib
- Pandas
- Plotly

Install the required packages using:

## Scripts
- `index.py` : Executes the main analysis pipeline
- `cifFileParser.py`: TFSI structural analysis functions
- `render.py`: Plotting and visualization functions
- `heatmap.py`: Interactive rendering
