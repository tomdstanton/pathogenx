# PathoGenX Web App
###### A web-based GUI for pathogen genotype exploration

[!WARNING]
🚧 
This package is currently under construction, proceed with caution
🚧

## Introduction
This app provides a web-based GUI for the exploration of pathogen genotyping data.
This begun as a Python port of the R-Shiny app and backend R code for the Klebsiella Neonatal Sepsis Sero-epi app,
but we generalised it for the exploration of different genotypes from multiple sources.
You can read more about the `pathogenx` code in the main README.md.

## Installation

The app is an optional module under the `pathogenx` library and can therefore be installed with pip:
```shell
pip install pathogenx[app]
```

## Usage
We've exposed the app module via the PathoGenX CLI, and can be run in your browser with the following command:
```shell
pathogenx app
```