from pathlib import Path
from typing import Generator

import pandas as pd

from pathogenx.io import InputFile

# Classes --------------------------------------------------------------------------------------------------------------
class GenotypeFile(_InputFile):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._loaders = {
            'generic': _load_table,
            'pathogenwatch': _load_pathogenwatch_genotype,
            'amrfinderplus': _load_tsv,
            'kaptive': _load_tsv,
        }


# Functions ------------------------------------------------------------------------------------------------------------
def _load_table(filepath: Path, index_col: int = 0, sep=None) -> pd.DataFrame:
    """This is a square symmetrical matrix; first line starts with Name\t"""
    return pd.read_table(filepath, index_col=index_col, sep=sep)


def _load_csv(filepath: Path, index_col: int = 0) -> pd.DataFrame:
    """This is a square symmetrical matrix; first line starts with Name\t"""
    return _load_table(filepath, index_col=index_col, sep=',')


def _load_tsv(filepath: Path, index_col: int = 0) -> pd.DataFrame:
    """This is a square symmetrical matrix; first line starts with Name\t"""
    return _load_table(filepath, index_col=index_col, sep='\t')


def _load_pathogenwatch_genotype(filepath: Path) -> pd.DataFrame:
    """This is a square symmetrical matrix; first line starts with Name\t"""
    return _load_csv(filepath, index_col=1)

