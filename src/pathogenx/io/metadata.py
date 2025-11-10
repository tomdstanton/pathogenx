




# Classes --------------------------------------------------------------------------------------------------------------
class MetadataFile(_InputFile):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._loaders = {
            'generic': _load_table,
            'pathogenwatch': _load_csv
        }




# Functions ------------------------------------------------------------------------------------------------------------
def _load_distance(filepath: Path, shape: Literal['square', 'long', 'guess'] = 'guess', sep: str = '\t',
                   usecolumns: tuple[int, int, int] = (0, 1, 2), symmetrical: bool = True) -> tuple[coo_matrix, list[str]]:
    if shape == 'square':
        return _load_distance_square(filepath, sep, symmetrical)
    elif shape == 'long':
        return _load_distance_long(filepath, sep, usecolumns)
    else:
        df = pd.read_table(filepath, sep=sep, index_col=0)
        nrows, ncolumns = df.shape
        if nrows == ncolumns:
            return _load_distance_square(filepath, sep, df.index.equals(df.columns))
        else:
            return _load_distance_long(filepath, sep, usecolumns)


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


def _load_mash_distance(filepath: Path) -> tuple[coo_matrix, list[str]]:
    """
    Parses a file with lines in the following format:
    genome1.fna	genome2.fna	0.0222766	0	456/1000
    We only need the first three columns.
    """
    return _load_distance_long(filepath, usecolumns=(0, 1, 3))


def _load_ska1_distance(filepath: Path) -> tuple[coo_matrix, list[str]]:
    """
    Sample 1	The name of the first sample being compared
    Sample 2	The name of the first sample being compared
    Matches	Number of split kmers found in both samples where the middle base is an A, C, G or T and matches between samples
    Mismatches	Number of split kmers found in only one of the samples
    Jaccard Index	Ratio of split kmers found in both samples to the total found in the two samples: matches/(matches+mismatches)
    Mash-like distance	A distance based on the Mash distance calculation using the Jaccard Index (j) above and the split kmer length (k): (-1/(2k+1))*ln(2j/(1+j)) for 0<j≤1 or 1 for j=0
    SNPs	Number of split kmers found in both samples where the middle base is an A, C, G or T but differs between files
    SNP distance	The ratio of SNPs to matches: SNPs/matches
    """
    return _load_distance_long(filepath, usecolumns=(0, 1, 6))


def _load_ska2_distance(filepath: Path) -> tuple[coo_matrix, list[str]]:
    return _load_distance_long(filepath)


def _load_pathogenwatch_distance(filepath: Path) -> tuple[coo_matrix, list[str]]:
    """This is a square symmetrical matrix; first line starts with Name\t"""
    return _load_distance_square(filepath, sep=',')


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


def load_pathogenwatch_datasets(path: Path) -> Generator[Dataset, None, None]:
    r = regex(r'.*pathogenwatch-(?P<species>\w+)-(?P<collection>[\w-]+)-'
              r'(?P<analysis>(kleborate|difference-matrix|metadata))\.csv')
    if not (files := [match for file in path.glob('*.csv') if (match := r.match(file.name))]):
        raise DatasetError(f'Could not find any files in {path}')
    for dataset, files in grouper(files, 2):
        files = {k: path / next(v).string for k, v in grouper(files, 3)}
        genotypes, metadata, distances = files.get('kleborate'), files.get('metadata'), files.get('difference-matrix')
        if genotypes is None:
            raise DatasetError(f'Could not find any genotypes in {path} for {dataset}')
        yield Dataset(
            _load_pathogenwatch_genotype(genotypes), _load_csv(metadata) if metadata else None,
            _load_pathogenwatch_distance(distances) if distances else None,
            name=dataset
        )

