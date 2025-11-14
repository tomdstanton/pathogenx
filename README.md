# PathoGenX 🦠🧬🗺️
###### A Python library for Pathogen Genotype eXploration

[![License](https://img.shields.io/badge/license-GPLv3-blue.svg?style=flat-square&maxAge=2678400)](https://choosealicense.com/licenses/gpl-3.0/)
[![PyPI](https://img.shields.io/pypi/v/pathogenx.svg?style=flat-square&maxAge=3600&logo=PyPI)](https://pypi.org/project/pathogenx)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/pathogenx?style=flat-square&maxAge=3600&logo=anaconda)](https://anaconda.org/bioconda/pathogenx)
[![Wheel](https://img.shields.io/pypi/wheel/pathogenx.svg?style=flat-square&maxAge=3600)](https://pypi.org/project/pathogenx/#files)
[![Python Versions](https://img.shields.io/pypi/pyversions/pathogenx.svg?style=flat-square&maxAge=600&logo=python)](https://pypi.org/project/pathogenx/#files)
[![Python Implementations](https://img.shields.io/pypi/implementation/pathogenx.svg?style=flat-square&maxAge=600&label=impl)](https://pypi.org/project/pathogenx/#files)
[![Source](https://img.shields.io/badge/source-GitHub-303030.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/pathogenx/)
[![Issues](https://img.shields.io/github/issues/tomdstanton/pathogenx.svg?style=flat-square&maxAge=600)](https://github.com/tomdstanton/pathogenx/issues)
[![Docs](https://img.shields.io/readthedocs/pathogenx/latest?style=flat-square&maxAge=600)](https://pathogenx.readthedocs.io)
[![Changelog](https://img.shields.io/badge/keep%20a-changelog-8A0707.svg?maxAge=2678400&style=flat-square)](https://github.com/tomdstanton/pathogenx/blob/main/CHANGELOG.md)
[![Downloads](https://img.shields.io/pypi/dm/pathogenx?style=flat-square&color=303f9f&maxAge=86400&label=downloads)](https://pepy.tech/project/pathogenx)

## Introduction 🌐
`pathogenx` is a Python library for Pathogen Genotype eXploration.

## Installation ⚙️
**NOTE: pathogenx is not yet on PyPI or Bioconda, please install from source until it is released**

### From source:
```shell
# First clone the repo
git clone https://github.com/tomdstanton/pathogenx.git && cd pathogenx
# Then install with pip
pip install .  # -e for editable, developers only!
```

## Usage 🧑‍💻
The information below explains how to use the `pathogenx` CLI. 
For API usage, please refer to the [reference documentation](https://tomdstanton.github.io/pathogenx/reference/pathogenx/).

### Prevalence 🌎

#### Quickstart
`pathogenx prevalence kleborate.csv --stratify-by K_locus ST > prevalence.tsv`

#### Arguments
```shell
usage: pathogenx prevalence <genotype> [metadata] [distance] [options]

========================|> PathoGenX |>========================
      A Python library for Pathogen Genotype eXploration       

Inputs:

  <genotypes>           Genotype file
  <metadata>            Optional metadata file
  <distances>           Optional distance file
  --genotype-flavour    Genotype file flavour (default: pw-kleborate)
                        (choices: kleborate, pw-kleborate, kaptive)
  --metadata-flavour    Metadata file flavour (default: pw-metadata)
                        (choices: pw-metadata)
  --distance-flavour    Distance file flavour (default: pw-dist)
                        (choices: mash, ska1, ska2, pw-dist)

Calculator options:

  --stratify-by  [ ...]
                        List of columns to stratify the analysis by
  --adjust-for [ ...]   Optional list of columns for adjustment (e.g., Cluster)
  --n-distinct [ ...]   Optional list of columns to calculate distinct counts for
  --denominator         Optional column to use as the primary grouping for denominators
                        If None, the first column in stratify-by is used

Clustering options:

  --snp-distance        The maximum distance for two samples to be considered connected.
                        Only used when `method` is connected_components

Other options:

  -v, --version         Show version number and exit
  -h, --help            Show this help message and exit
```

