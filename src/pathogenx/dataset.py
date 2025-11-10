from pathlib import Path
from typing import Literal, LiteralString, Union, Callable, Generator
from re import compile as regex
from warnings import warn

import pandas as pd
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components


# Classes --------------------------------------------------------------------------------------------------------------
class DatasetError(Exception):
    pass


class DatasetWarning(Warning):
    pass


class Dataset:
    def __init__(self, genotypes: pd.DataFrame, metadata: pd.DataFrame = None,
                 distances: tuple[coo_matrix, list[str]] = None, name: str = 'unknown', genotype_columns: set[str] = None,
                 qc_metrics: set[str] = None, metadata_columns: set[str] = None):
        # Add genotype data and set dataset name ----------------------------------------
        self._data = genotypes
        self._data['Dataset'] = name
        # Add metadata ------------------------------------------------------------------
        if metadata is not None:
            self._data = self._data.join(metadata, how='left')
        # Add distances ------------------------------------------------------------------
        if distances is not None:
            distances, index = distances  # Unpack tuple
            for sample in self._data.index:  # Check all samples are present in data
                if sample not in index:
                    raise DatasetError(f'Sample {sample} not in index')
            self._data = self._data.reindex(index)  # Reorder index to match distance matrix
            self.distances = distances.tocsr()  # Convert to CSR for efficient row slicing
        else:
            self.distances = None
        # Add columns -------------------------------------------------------------------
        self.metadata_columns: set[str] = metadata_columns or (set(metadata.columns) if metadata is not None else set())
        self.metadata_columns.add('Dataset')
        self.genotype_columns: set[str] = (genotype_columns or set(
            genotypes.columns)) - self.qc_metrics - self.metadata_columns

    def __repr__(self) -> str:
        return (f'Dataset({len(self._data)} samples {"with" if self.distances is not None else "without"} distances, '
                f'{len(self.genotype_columns)} genotypes, {len(self.qc_metrics)} QC metrics, '
                f'{len(self.metadata_columns)} metadata variables)')

    def __len__(self):
        return len(self._data)

    def __contains__(self, sample: str) -> bool:
        return sample in self._data.index

    def __getitem__(self, item):
        return self._data[item]

    def __iter__(self):
        return self._data.itertuples()

    @property
    def data(self) -> pd.DataFrame:
        """
        Returns a copy of the internal DataFrame to prevent unintended modification of the Dataset's state.
        """
        return self._data.copy()

    def samples(self):
        return self._data.index

    def calculate_clusters(self, method: Literal['connected_components', 'variables'] = 'connected_components',
                           group_by: list[str] = None, distance: int = 20) -> pd.Series:
        """Calculates clusters and adds/overwrites the 'Cluster' column.

        This method provides two main strategies for clustering: finding connected
        components in a distance graph or grouping by categorical variables. The
        resulting cluster labels are stored in a 'Cluster' column in the internal
        DataFrame. If the column already exists, it will be overwritten and a
        `DatasetWarning` will be issued.

        Args:
            method (Literal['connected_components', 'variables'], optional):
                The clustering strategy to use. Defaults to 'connected_components'.
                - 'connected_components': Finds clusters based on a distance
                  threshold. Requires a distance matrix.
                - 'variables': Groups samples based on the specified `group_by`
                  columns.
            group_by (list[str], optional): A list of column names to group by.
                If `None` with 'variables', each sample gets a unique cluster.
                If `None` with 'connected_components', the entire dataset is
                treated as a single group. Defaults to None.
            distance (int, optional): The maximum distance for two samples to be
                considered connected. Only used when `method` is
                'connected_components'. Defaults to 20.

        Returns:
            A pandas Series containing the calculated cluster labels for each sample.

        Raises:
            DatasetError: If `method` is 'connected_components' and no distance
                matrix is available in the Dataset.
            ValueError: If an unknown `method` is provided.
        """
        if 'Cluster' in self._data.columns:
            warn("'Cluster' column already exists and will be overwritten.", DatasetWarning)

        if method == 'variables':
            self._data['Cluster'] = 'cluster_' + (self._data.groupby(group_by).ngroup() + 1).astype(str) \
                if group_by else [f'cluster_{i + 1}' for i in range(len(self._data))]

        elif method == 'connected_components':
            if self.distances is None:
                raise DatasetError("Cannot calculate connected components: no distance matrix found in Dataset.")

            cluster_labels = pd.Series(index=self._data.index, dtype='object')
            global_cluster_counter = 1

            # Determine groups: either from group_by or a single group for the whole dataset
            groups = self._data.groupby(group_by) if group_by else [('all', self._data)]

            for _, group_df in groups:
                if group_df.empty:
                    continue

                # Get integer indices for samples in this group
                group_indices = self._data.index.get_indexer(group_df.index)

                # Efficiently subset the sparse graph and filter by distance
                subgraph = self.distances[group_indices, :][:, group_indices]
                subgraph.data[subgraph.data > distance] = 0
                subgraph.eliminate_zeros()

                _, local_labels = connected_components(csgraph=subgraph, directed=False, return_labels=True)

                if local_labels.size > 0:
                    cluster_labels.loc[group_df.index] = [f"cluster_{l + global_cluster_counter}" for l in local_labels]
                    global_cluster_counter += local_labels.max() + 1

            self._data['Cluster'] = cluster_labels

        else:
            raise ValueError(f"Unknown method: {method}")

        return self._data['Cluster']


# Functions ------------------------------------------------------------------------------------------------------------
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

