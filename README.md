# TFSI Crystallography Analysis

This repository contains Python code for analyzing X-ray crystallographic data of bis(trifluoromethanesulfonyl)imide (TFSI) compounds. View at:

[Graphs Dashboard Here](https://tfsi-research.vercel.app/)

## Overview

Chemical data science is an emerging field within the broader field of chemistry that has numerous applications of high relevance to a variety of academic and industrial pursuits. With quantum computing and artificial intelligence becoming more mainstream, it is of great interest for not only the data scientist, but the chemist as well, to take advantage of these technologies to streamline data processing, analyze large data sets, and reveal new chemical insights that would be otherwise hidden by the sheer amount and depth/complexity of available data. In this project, modern computing and fundamental chemical analysis have been combined in order to pursue new insights into the structural behavior of the weakly coordinating anion bis(trifluoromethylsulfonyl)imide, otherwise known as TFSI. Taking a data science approach, published solid-state crystal structures available in the Cambridge Structural Database (CSD) including one or more TFSI species of interest were categorized and statistically analyzed using software built for this purpose. The goal of this project from a chemical perspective was to determine the structural characteristics displayed by TFSI as inferred from the structural data. The goal of this project from a data-science perspective was to develop a new software program using Python to parse Crystallographic Information Files (CIF) obtained from the CSD into statistically relevant information that could be compared across the individual structural data sets. This research endeavor aims to highlight the applications of data science to an otherwise foreign area of research (structural chemistry) and outline the capabilities of data processing for future structural and chemical investigations.  

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
