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

> [!WARNING]
> 🚧 This package is currently under construction, proceed with caution 🚧

## Installation ⚙️

### From source:
```shell
pip install git+https://github.com/tomdstanton/pathogenx.git
```

## Usage 🧑‍💻
The information below explains how to use the `pathogenx` CLI. 
For API usage, please refer to the [reference documentation](https://tomdstanton.github.io/pathogenx/reference/pathogenx/).

### Prevalence 🌎

#### Quickstart
To calculate prevalences from the 
[KlebNet neonatal sepsis collection](http://pathogen.watch/collection/klebnet-neonatal-sepsis),
download Kleborate, metadata and distance matrix files from Pathogenwatch, then run the following:

```shell
pathogenx prevalence pw_kleborate.csv ST > prevalence.tsv
pathogenx prevalence pw_kleborate.csv Country ST --metadata pw_metadata.csv > prevalence_by_country.tsv
pathogenx prevalence pw_kleborate.csv K_locus --distances pw_distances.csv --adjust-for Cluster > adjusted_prevalence.tsv
```

#### Arguments
```shell
usage: pathogenx prevalence <genotypes> <strata> [options]

========================|> PathoGenX |>========================
      A Python library for Pathogen Genotype eXploration       

Inputs:

  <genotypes>          Genotype file
  <strata>             List of columns to stratify the analysis by
  --metadata []        Optional metadata file
  --distances []       Optional distance file
  --genotype-flavour   Genotype file flavour (default: pw-kleborate)
                       (choices: kleborate, pw-kleborate, kaptive)
  --metadata-flavour   Metadata file flavour (default: pw-metadata)
                       (choices: pw-metadata)
  --distance-flavour   Distance file flavour (default: pw-dist)
                       (choices: mash, ska1, ska2, pw-dist)

Calculator options:

  --adjust-for [ ...]  Optional list of columns for adjustment (e.g., Cluster)
  --n-distinct [ ...]  Optional list of columns to calculate distinct counts for
  --denominator        Optional column to use as the primary grouping for denominators
                       If None, the first column in stratify-by is used

Clustering options:

  --snp-distance       The maximum distance for two samples to be considered connected.
                       Only used when `method` is connected_components (default: 20)

Other options:

  -v, --version        Show version number and exit
  -h, --help           Show this help message and exit
```

## Web-app
The PathoGenX app provides a web-based GUI for the exploration of pathogen genotyping data. It is an optional module
that can be installed like so:

```shell
pip install pathogenx[app]
```

Read more about the app [here](src/pathogenx/app/README.md)
