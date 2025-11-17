from pathlib import Path
from typing import Union, Literal, get_args, TypeVar
from abc import ABC, abstractmethod
import pandas as pd
from scipy.sparse import coo_matrix


# Constants ------------------------------------------------------------------------------------------------------------
_GENOTYPE_FLAVOURS = Literal['pw-kleborate', 'kleborate', 'kaptive']
_META_FLAVOURS = Literal['pw-metadata']
_DIST_FLAVOURS = Literal['pw-dist', 'mash', 'ska1', 'ska2']


# Classes --------------------------------------------------------------------------------------------------------------
class _InputFile(ABC):
    def __init__(self, filepath: Union[Path, str], name: str = None):
        self.filepath = filepath if isinstance(filepath, Path) else Path(filepath)
        self.name = name or self.filepath.stem

    @abstractmethod
    def load(self) -> Union[pd.DataFrame, tuple[coo_matrix, list[str]]]: pass

    def __repr__(self):
        return f"{self.__class__.__name__}('{self.name}')"


class _BaseTextFile(_InputFile):
    def __init__(self, filepath: Union[Path, str], sep: str = '\t', index_col: int = 0):
        super().__init__(filepath=filepath)
        self.index_col = index_col
        self.sep = sep

    def load(self) -> pd.DataFrame:
        return pd.read_table(self.filepath, sep=self.sep, index_col=self.index_col)
    

class GenotypeFile(_BaseTextFile):
    def __init__(self, filepath: Union[Path, str], sep: str = '\t', index_col: int = 0):
        super().__init__(filepath=filepath, sep=sep, index_col=index_col)

    @classmethod
    def from_flavour(cls, filepath: Union[Path, str], flavour: _GENOTYPE_FLAVOURS) -> 'GenotypeFile':
        if flavour in {'kleborate', 'kaptive'}:
            self = cls(filepath)
        elif flavour == 'pw-kleborate':
            self = cls(filepath, sep=',', index_col=1)
        else:
            raise ValueError(f"Unknown flavour: {flavour}")
        return self


class MetaFile(_BaseTextFile):
    def __init__(self, filepath: Union[Path, str], sep: str = '\t', index_col: int = 0):
        super().__init__(filepath=filepath, sep=sep, index_col=index_col)

    @classmethod
    def from_flavour(cls, filepath: Union[Path, str], flavour: _META_FLAVOURS) -> 'MetaFile':
        if flavour == 'pw-metadata':
            self = cls(filepath, sep=',')
        else:
            raise ValueError(f"Unknown flavour: {flavour}")
        return self


class DistFile(_InputFile):
    def __init__(self, filepath: Union[Path, str], shape: Literal['square', 'long'] = 'long',
                 sep: str = '\t', usecols: tuple[int, int, int] = (0, 1, 2), symmetrical: bool = True):
        super().__init__(filepath=filepath)
        self.shape = shape
        self.sep = sep
        self.usecols = usecols
        self.symmetrical = symmetrical

    def load(self) -> tuple[coo_matrix, list[str]]:
        """Loads distance matrix, guessing the format if no flavour was provided."""
        if self.shape == 'square':
            df = pd.read_table(self.filepath, sep=self.sep, index_col=0)
            matrix = coo_matrix(df.values)  # Convert the pandas DataFrame to a scipy sparse COO matrix
            if not self.symmetrical:
                matrix = matrix.maximum(matrix.T)
            return matrix, df.columns.tolist()
        elif self.shape == 'long':
            df = pd.read_table(self.filepath, sep=self.sep, usecols=self.usecols, header=None)
            index = {sample: i for i, sample in enumerate(pd.unique(df.iloc[:, [0, 1]].values.ravel('K')))}
            num_samples = len(index)
            rows = df.iloc[:, 0].map(index).values
            columns = df.iloc[:, 1].map(index).values
            data = df.iloc[:, 2].values
            matrix = coo_matrix((data, (rows, columns)), shape=(num_samples, num_samples))
            return matrix, list(index.keys())
        else:
            raise ValueError(f"Unknown shape: {self.shape}")

    @classmethod
    def from_flavour(cls, filepath: Union[Path, str], flavour: _DIST_FLAVOURS) -> 'DistFile':
        if flavour == 'mash':
            self = cls(filepath, shape='long', usecols=(0, 1, 3), sep='\t')
        elif flavour == 'ska1':
            self = cls(filepath, shape='long', usecols=(0, 1, 6), sep='\t')
        elif flavour == 'ska1':
            self = cls(filepath, shape='long', usecols=(0, 1, 2), sep='\t')
        elif flavour == 'pw-dist':
            self = cls(filepath, shape='square', sep=',', symmetrical=True)
        else:
            raise ValueError(f"Unknown flavour: {flavour}")
        return self

