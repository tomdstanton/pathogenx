"""
Module for calculating insights from genotyping data such as prevalence.
"""

from abc import ABC, abstractmethod
from typing import List, Optional, Union
import pandas as pd
import numpy as np
from scipy.stats import norm

from pathogenx.dataset import Dataset
# from .models import ModelResult

# Constants ------------------------------------------------------------------------------------------------------------
# _ALPHA_DIVERSITY_METRICS = LiteralString['simpson', 'simpson_e', 'simpson_d']


# Classes --------------------------------------------------------------------------------------------------------------
class CalculatorResult(ABC):
    """
    Abstract base class for calculator results.

    Designed such that we can infer the parameters of the calculator
    used to generate the result, without needing the instance itself.
    """
    def __init__(self):
        """Initializes the CalculatorResult."""
        self._data = None

    @property
    def data(self) -> pd.DataFrame:
        """Returns a copy of the result data to prevent accidental modification."""
        return self._data.copy() if self._data is not None else pd.DataFrame()

    @data.setter
    def data(self, value: pd.DataFrame):
        """Sets the result data stored as a private property."""
        self._data = value

    def __len__(self):
        """Returns the number of rows in the result data.

        Returns:
            int: The length of the result data DataFrame.
        """
        return 0 if self._data is None else len(self._data)


class PrevalenceResult(CalculatorResult):
    """
    Class to store the results of a `PrevalenceCalculator`.

    Attributes:
        stratified_by (list[str]): List of columns the result is stratified by,
            in order of strata level.
        adjusted_for (list[str], optional): List of columns the prevalences were
            adjusted for (e.g., 'Cluster').
        n_distinct (list[str], optional): List of columns for which distinct counts
            (per-strata) were generated.
        denominator (str, optional): Column indicating the denominator stratum.
    """
    def __init__(self, stratified_by: List[str], adjusted_for: List[str] = None, n_distinct: List[str] = None,
                 denominator: str = None):
        """Initializes the PrevalenceResult.

        Args:
            stratified_by (List[str]): Columns the result is stratified by.
            adjusted_for (List[str], optional): Columns prevalences were adjusted for.
                Defaults to None.
            n_distinct (List[str], optional): Columns distinct counts were generated for.
                Defaults to None.
            denominator (str, optional): The denominator stratum. Defaults to None.
        """
        super().__init__()
        self.stratified_by: List[str] = stratified_by
        self.adjusted_for: Optional[List[str]] = adjusted_for
        self.n_distinct: Optional[List[str]] = n_distinct
        self.denominator: Optional[str] = denominator

    @classmethod
    def from_calculator(cls, calculator: 'PrevalenceCalculator') -> 'PrevalenceResult':
        """Initializes a PrevalenceResult from a PrevalenceCalculator instance.

        Args:
            calculator (PrevalenceCalculator): The calculator instance.

        Returns:
            PrevalenceResult: A new PrevalenceResult object configured with the
                calculator's parameters.
        """
        return cls(calculator.stratify_by, calculator.adjust_for, calculator.n_distinct, calculator.denominator)


class Calculator(ABC):
    """Abstract base class for all calculators."""
    def __init__(self, model: 'ModelResult' = None):
        """Initializes the Calculator.

        Args:
            model (ModelResult, optional): A model result to be used by the
                calculator. Defaults to None.
        """
        self.model = model

    @abstractmethod
    def calculate(self, dataset: Dataset):
        """Abstract method to perform a calculation on a dataset.

        Args:
            dataset (Dataset): The dataset to perform the calculation on.
        """
        pass


class PrevalenceCalculator(Calculator):
    """
    Calculates raw and adjusted prevalence statistics from a dataset.

    This class is designed to replicate the logic of the R 'prevalence'
    package using pandas for efficient data manipulation.

    Attributes:
        stratify_by (list[str]): List of columns to stratify the analysis by.
        adjust_for (list[str], optional): List of columns for adjustment (e.g., 'Cluster').
        n_distinct (list[str], optional): List of columns to calculate distinct counts for.
        denominator (str, optional): Column to use as the primary grouping for
            denominators. If None, the first column in stratify_by is used.
    """
    def __init__(self, stratify_by: List[str], adjust_for: List[str] = None, n_distinct: List[str] = None,
                 denominator: str = None):
        """Initializes the PrevalenceCalculator.

        Args:
            stratify_by (List[str]): Columns to stratify the analysis by.
            adjust_for (List[str], optional): Columns for adjustment. Defaults to None.
            n_distinct (List[str], optional): Columns to calculate distinct counts for.
                Defaults to None.
            denominator (str, optional): Column for primary grouping for denominators.
                Defaults to None.

        Raises:
            ValueError: If 'denominator' is provided but not in 'stratify_by'.
        """
        super().__init__()

        if denominator and denominator not in stratify_by:
            raise ValueError("If provided, 'denominator' must be in 'stratify_by'.")

        # Set denominator if not provided and stratification has more than one level
        if not denominator and len(stratify_by) > 1:
            denominator = stratify_by[0]

        self.stratify_by: List[str] = stratify_by
        self.adjust_for: Optional[List[str]] = adjust_for
        self.n_distinct: Optional[List[str]] = n_distinct
        self.denominator: Optional[str] = denominator

    def calculate(self, dataset: Union[Dataset, pd.DataFrame]) -> PrevalenceResult:
        """Calculates prevalence statistics on a dataset.

        This method computes raw and optionally adjusted counts, proportions,
        standard errors, and 95% confidence intervals (using the Wilson score
        interval method) for the specified strata. It can also calculate
        distinct counts for designated columns and rank prevalences within
        denominator groups.

        Args:
            dataset (Union[Dataset, pd.DataFrame]): The dataset on which to
                calculate prevalence.

        Returns:
            PrevalenceResult: An object containing the calculated prevalence data
                and metadata about the calculation.

        Raises:
            TypeError: If the dataset is not a pandas DataFrame or a
                pathogenx Dataset object.
        """
        if isinstance(dataset, Dataset):
            data = dataset.data  # Gets a copy from the @property
        elif isinstance(dataset, pd.DataFrame):
            data = dataset.copy()
        else:
            raise TypeError("dataset must be of type Dataset or pd.DataFrame")

        # 1. Calculate Denominators
        denominators = {}
        adj_col = self.adjust_for[0] if self.adjust_for else None
        if self.denominator:
            denom_groups = data.groupby(self.denominator)
            denominators['raw'] = denom_groups.size().rename('denominator.raw')
            if adj_col:
                denominators['adj'] = denom_groups[adj_col].nunique().rename('denominator.adj')
        else:
            denominators['raw'] = pd.Series([len(data)], name='denominator.raw')
            if adj_col:
                denominators['adj'] = pd.Series([data[adj_col].nunique()], name='denominator.adj')

        # 2. Calculate Counts within strata
        strata_groups = data.groupby(self.stratify_by)
        result_data = strata_groups.size().to_frame('count.raw')
        if adj_col:
            result_data['count.adj'] = strata_groups[adj_col].nunique()
        if self.n_distinct:
            for col in self.n_distinct:
                result_data[f'# {col}'] = strata_groups[col].nunique()

        result_data = result_data.reset_index()

        # 3. Join denominators and calculate proportions
        if self.denominator:
            result_data = result_data.merge(denominators['raw'], on=self.denominator, how='left')
            if adj_col:
                result_data = result_data.merge(denominators['adj'], on=self.denominator, how='left')
        else:
            result_data['denominator.raw'] = denominators['raw'][0]
            if adj_col:
                result_data['denominator.adj'] = denominators['adj'][0]

        # 4. Calculate Proportions, SE, and CI
        for col_type in denominators:
            count = result_data[f'count.{col_type}']
            denom = result_data[f'denominator.{col_type}']
            for x, y in zip(('prop', 'se', 'lower', 'upper'), _wilson_score_interval(count, denom)):
                result_data[f'{x}.{col_type}'] = y

        # 5. Sort data
        sort_by = ['denominator.adj', 'count.adj'] if self.adjust_for else ['denominator.raw', 'count.raw']
        result_data = result_data.sort_values(by=sort_by, ascending=False)

        # 6. Calculate Ranks
        rank_groups = result_data.groupby(self.denominator) if self.denominator else result_data
        for col_type in denominators:
            result_data[f'rank.{col_type}'] = rank_groups[f'prop.{col_type}'].rank(method='first', ascending=False)

        result = PrevalenceResult.from_calculator(self)
        result.data = result_data  # Set the result data
        return result


# Functions ------------------------------------------------------------------------------------------------------------
def _wilson_score_interval(counts: np.ndarray, denominators: np.ndarray) -> tuple[float, float, float, float]:
    """Calculates the Wilson score interval for a proportion.

    Also returns the simple proportion and standard error (using Wald method)
    for reporting alongside the more robust Wilson interval.

    Args:
        counts (np.ndarray): The number of successes (numerator).
        denominators (np.ndarray): The total number of trials (denominator).

    Returns:
        tuple[float, float, float, float]: A tuple containing the proportion,
            standard error, lower bound of the CI, and upper bound of the CI.
    """
    prop = (counts / denominators).clip(0, 1)
    # Use the more robust Wilson score interval for CI
    z = norm.ppf(1 - (0.05 / 2))  # Z-score for 95% CI
    # Calculate standard error using the Wald method for reporting
    se = np.sqrt(prop * (1 - prop) / denominators)
    # Wilson score interval calculation
    center = (counts + z ** 2 / 2) / (denominators + z ** 2)
    width = (z / (denominators + z ** 2)) * np.sqrt((prop * (1 - prop) * denominators) + (z ** 2 / 4))
    return prop, se, center - width, center + width


def _calculate_ci(samples: np.ndarray, prob: float = 0.95) -> tuple[np.ndarray, np.ndarray]:
    """Calculates quantile-based credible intervals.

    Args:
        samples (np.ndarray): An array of samples from a posterior distribution.
        prob (float, optional): The desired probability for the credible interval.
            Defaults to 0.95.

    Returns:
        tuple[np.ndarray, np.ndarray]: A tuple containing the lower and upper
            bounds of the credible interval.
    """
    lower_q = (1 - prob) / 2
    upper_q = 1 - lower_q
    # Quantiles are calculated over the flattened chain/draw dimensions
    flat_samples = samples.reshape(-1, samples.shape[-1]) if samples.ndim > 1 else samples.flatten()
    return np.quantile(flat_samples, lower_q, axis=0), np.quantile(flat_samples, upper_q, axis=0)
