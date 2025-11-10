from pathlib import Path
from typing import LiteralString, Union, Callable

import pandas as pd
from scipy.sparse import coo_matrix


# Classes --------------------------------------------------------------------------------------------------------------
class InputFileError(Exception):
    pass


class InputFile:
    def __init__(self, filepath: Path | str, name: str = None, flavour: LiteralString = 'generic'):
        self.filepath = filepath if isinstance(filepath, Path) else Path(filepath)
        self.name = name or self.filepath.stem
        self.flavour = flavour
        self._loaders: dict[LiteralString, Callable] = {}

    def __call__(self, *args, **kwargs) -> Union[pd.DataFrame, tuple[coo_matrix, list[str]]]:
        if loader := self._loaders.get(self.flavour):
            return loader(self.filepath, *args, **kwargs)
        raise InputFileError(f'Could not load {self.name} using flavour {self.flavour}')

