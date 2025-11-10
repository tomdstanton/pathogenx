from pathlib import Path
from typing import Union, Literal
from abc import ABC, abstractmethod

import pandas as pd
from scipy.sparse import coo_matrix


# Classes --------------------------------------------------------------------------------------------------------------
class InputFile(ABC):
    def __init__(self, filepath: Path | str, name: str = None):
        self.filepath = filepath if isinstance(filepath, Path) else Path(filepath)
        self.name = name or self.filepath.stem

    @abstractmethod
    def load(self) -> Union[pd.DataFrame, tuple[coo_matrix, list[str]]]: pass


class TextFile(InputFile):
    def __init__(self, filepath: Path, sep: str = '\t', index_col: int = 0):
        super().__init__(filepath)
        self.index_col = index_col
        self.sep = sep

    def load(self) -> pd.DataFrame:
        return pd.read_table(self.filepath, sep=self.sep, index_col=self.index_col)


class DistFile(InputFile):
    def __init__(self, filepath: Path, shape: Literal['square', 'long', 'guess'] = 'guess', sep: str = '\t',
                   usecolumns: tuple[int, int, int] = (0, 1, 2), symmetrical: bool = True):
        super().__init__(filepath)
        self.shape = shape
        self.sep = sep
        self.usecolumns = usecolumns
        self.symmetrical = symmetrical

    def load(self) -> tuple[coo_matrix, list[str]]:
        if self.shape == 'square':
            return _load_distance_square(self.filepath, self.sep, self.symmetrical)
        elif self.shape == 'long':
            return _load_distance_long(self.filepath, self.sep, self.usecolumns)
        else:
            df = pd.read_table(self.filepath, sep=self.sep, index_col=0)
            nrows, ncolumns = df.shape
            if nrows == ncolumns:
                return _load_distance_square(self.filepath, self.sep, df.index.equals(df.columns))
            else:
                return _load_distance_long(self.filepath, self.sep, self.usecolumns)


class MashFile(DistFile):
    def __init__(self, filepath: Path):
        super().__init__(filepath, usecolumns=(0, 1, 3))


class Ska1File(DistFile):
    def __init__(self, filepath: Path):
        super().__init__(filepath, shape='long', usecolumns=(0, 1, 6))


class Ska2File(DistFile):
    def __init__(self, filepath: Path):
        super().__init__(filepath, shape='long')


class PathogenwatchFile(DistFile):
    def __init__(self, filepath: Path):
        super().__init__(filepath, shape='square', symmetrical=True, sep=',')



# Functions ------------------------------------------------------------------------------------------------------------
def _load_distance_square(filepath: Path, sep: str = '\t', symmetrical: bool = True) -> tuple[coo_matrix, list[str]]:
    # Read the matrix using pandas, assuming the first column is the index
    df = pd.read_table(filepath, index_col=0, sep=sep)
    M = coo_matrix(df.values)  # Convert the pandas DataFrame to a scipy sparse COO matrix
    if not symmetrical:
        M = M.maximum(M.T)
    return M, df.columns.tolist()


def _load_distance_long(filepath: Path, sep: str = '\t', usecolumns: tuple[int, int, int] = (0, 1, 2)) -> tuple[
    coo_matrix, list[str]]:
    df = pd.read_table(filepath, usecolumns=usecolumns, sep=sep)
    # Get all unique sample names from the first two columns
    index = {sample: i for i, sample in enumerate(pd.unique(df.iloc[:, [0, 1]].values.ravel('K')))}
    num_samples = len(index)
    # Create a mapping from sample names to integer indices
    rows = df.iloc[:, 0].map(index).values
    columns = df.iloc[:, 1].map(index).values
    data = df.iloc[:, 2].values
    # Create the COO matrix
    matrix = coo_matrix((data, (rows, columns)), shape=(num_samples, num_samples))
    # The matrix should be symmetrical, so add the transpose as well
    # Ensure the matrix is square and covers all unique samples
    # matrix += coo_matrix((data, (columns, rows)), shape=(num_samples, num_samples))
    return matrix, list(index.keys())
