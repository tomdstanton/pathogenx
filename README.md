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

> [!WARNING]
> 🚧 This package is currently under construction, proceed with caution 🚧

## Introduction 🌐
`pathogenx` is a Python library for Pathogen Genotype eXploration. It started a Python port of the 
[KleborateR](https://github.com/klebgenomics/KleborateR) R code to parse 
[Kleborate](https://github.com/klebgenomics/Kleborate) and other data from [Pathogenwatch](https://pathogen.watch/),
to calculate prevalence data for sero-epidemiology, and to provide a backend for our 
[neonatal sepsis sero-epi app](https://github.com/klebgenomics/KlebNNSapp).

`pathogenx` aims to be more generalised towards other bacterial genotyping and distance calculation methods, and
will eventually force complicity with the 
[PHA4GE genotyping-specification](https://github.com/pha4ge/genotyping-specification).

`pathogenx` revolves around genotyping `Dataset` objects, which have the following attributes:

- Genotyping data - a `pandas` dataframe containing the parsed output of a genotyping tool. 
- Optional metadata - joined to genotyping results upon initialisation, possibly containing spatio-temporal data.
- Optional distances - Represented as a `scipy.sparse` matrix of pairwise distances, parsed from the outputs of tools 
such as `mash`.

We also define `Calculator`s, which calculate informative information from genotyping data, such as prevalence and 
diversity.

## Installation ⚙️

### From source:
```shell
pip install git+https://github.com/tomdstanton/pathogenx.git
```

## Usage 🧑‍💻
The information below explains how to use the `pathogenx` CLI. 
For API usage, please refer to the [reference documentation](https://tomdstanton.github.io/pathogenx/reference/pathogenx/).

### Prevalence 🌎
The `PrevalenceCalculator` object to calculate prevalence from a dataset is exposed via the command-line. For more
information about the calculator, please refer to the 
[API docs](https://tomdstanton.github.io/pathogenx/reference/pathogenx/calculators/#pathogenx.calculators.PrevalenceCalculator).

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
that can be installed and run with the following command:

```shell
pip install pathogenx[app]
pathogenx app
```

Read more about the app [here](src/pathogenx/app/README.md).
