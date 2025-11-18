"""
Module to read and parse genotype, metadata and distance data from different sources
"""
from pathlib import Path
from typing import Union, Literal
from abc import ABC, abstractmethod
import pandas as pd
from scipy.sparse import coo_matrix


# Constants ------------------------------------------------------------------------------------------------------------
_GENOTYPE_FLAVOURS = Literal['pw-kleborate', 'kleborate', 'kaptive']
_META_FLAVOURS = Literal['pw-metadata']
_DIST_FLAVOURS = Literal['pw-dist', 'mash', 'ska1', 'ska2']


# Classes --------------------------------------------------------------------------------------------------------------
class _InputFile(ABC):
    """Abstract base class for input files."""
    def __init__(self, filepath: Union[Path, str], name: str = None):
        """Initializes the InputFile.

        Args:
            filepath (Union[Path, str]): The path to the input file.
            name (str, optional): The name of the file, used for representation.
                Defaults to the file stem.
        """
        self.filepath = filepath if isinstance(filepath, Path) else Path(filepath)
        self.name = name or self.filepath.stem

    @abstractmethod
    def load(self) -> Union[pd.DataFrame, tuple[coo_matrix, list[str]]]:
        """Abstract method to load the file content."""
        pass

    def __repr__(self):
        """Return a string representation of the object."""
        return f"{self.__class__.__name__}('{self.name}')"


class _BaseTextFile(_InputFile):
    """Base class for delimited text files that can be loaded into a DataFrame."""
    def __init__(self, filepath: Union[Path, str], sep: str = '\t', index_col: int = 0):
        """Initializes the _BaseTextFile.

        Args:
            filepath (Union[Path, str]): The path to the input file.
            sep (str, optional): The separator used in the file. Defaults to '\t'.
            index_col (int, optional): The column to use as the index. Defaults to 0.
        """
        super().__init__(filepath=filepath)
        self.index_col = index_col
        self.sep = sep

    def load(self) -> pd.DataFrame:
        """Loads the text file into a pandas DataFrame.

        Returns:
            pd.DataFrame: The loaded data.
        """
        return pd.read_table(self.filepath, sep=self.sep, index_col=self.index_col)


class GenotypeFile(_BaseTextFile):
    """Represents a genotype file."""
    def __init__(self, filepath: Union[Path, str], sep: str = '\t', index_col: int = 0):
        """Initializes the GenotypeFile.

        Args:
            filepath (Union[Path, str]): The path to the genotype file.
            sep (str, optional): The separator used in the file. Defaults to '\t'.
            index_col (int, optional): The column to use as the index. Defaults to 0.
        """
        super().__init__(filepath=filepath, sep=sep, index_col=index_col)

    @classmethod
    def from_flavour(cls, filepath: Union[Path, str], flavour: _GENOTYPE_FLAVOURS) -> 'GenotypeFile':
        """Creates a GenotypeFile instance from a specific file format flavour.

        Args:
            filepath (Union[Path, str]): The path to the genotype file.
            flavour (_GENOTYPE_FLAVOURS): The format of the file.

        Returns:
            GenotypeFile: An instance of GenotypeFile with appropriate settings.

        Raises:
            ValueError: If the flavour is unknown.
        """
        if flavour in {'kleborate', 'kaptive'}:
            self = cls(filepath)
        elif flavour == 'pw-kleborate':
            self = cls(filepath, sep=',', index_col=1)
        else:
            raise ValueError(f"Unknown flavour: {flavour}")
        return self


class MetaFile(_BaseTextFile):
    """Represents a metadata file."""
    def __init__(self, filepath: Union[Path, str], sep: str = '\t', index_col: int = 0):
        """Initializes the MetaFile.

        Args:
            filepath (Union[Path, str]): The path to the metadata file.
            sep (str, optional): The separator used in the file. Defaults to '\t'.
            index_col (int, optional): The column to use as the index. Defaults to 0.
        """
        super().__init__(filepath=filepath, sep=sep, index_col=index_col)

    @classmethod
    def from_flavour(cls, filepath: Union[Path, str], flavour: _META_FLAVOURS) -> 'MetaFile':
        """Creates a MetaFile instance from a specific file format flavour.

        Args:
            filepath (Union[Path, str]): The path to the metadata file.
            flavour (_META_FLAVOURS): The format of the file.

        Returns:
            MetaFile: An instance of MetaFile with appropriate settings.

        Raises:
            ValueError: If the flavour is unknown.
        """
        if flavour == 'pw-metadata':
            self = cls(filepath, sep=',')
        else:
            raise ValueError(f"Unknown flavour: {flavour}")
        return self


class DistFile(_InputFile):
    """Represents a distance matrix file."""
    def __init__(self, filepath: Union[Path, str], shape: Literal['square', 'long'] = 'long',
                 sep: str = '\t', usecols: tuple[int, int, int] = (0, 1, 2), symmetrical: bool = True):
        """Initializes the DistFile.

        Args:
            filepath (Union[Path, str]): The path to the distance matrix file.
            shape (Literal['square', 'long'], optional): The shape of the matrix file.
                'square' is a traditional matrix, 'long' is a 3-column format. Defaults to 'long'.
            sep (str, optional): The separator used in the file. Defaults to '\t'.
            usecols (tuple[int, int, int], optional): The columns to use for long format
                (sample1, sample2, distance). Defaults to (0, 1, 2).
            symmetrical (bool, optional): Whether a square matrix is symmetrical. Defaults to True.
        """
        super().__init__(filepath=filepath)
        self.shape = shape
        self.sep = sep
        self.usecols = usecols
        self.symmetrical = symmetrical

    def load(self) -> tuple[coo_matrix, list[str]]:
        """Loads a distance matrix into a sparse matrix format.

        Returns:
            tuple[coo_matrix, list[str]]: A tuple containing the sparse distance
            matrix and a list of sample names.

        Raises:
            ValueError: If the shape is unknown.
        """
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
        """Creates a DistFile instance from a specific file format flavour.

        Args:
            filepath (Union[Path, str]): The path to the distance matrix file.
            flavour (_DIST_FLAVOURS): The format of the file.

        Returns:
            DistFile: An instance of DistFile with appropriate settings.

        Raises:
            ValueError: If the flavour is unknown.
        """
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
