"""
Module for calculating insights from genotyping data such prevalence
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
    Abstract base class for calculator results, designed such that we can infer the parameters of the calculator
    used to generate the result, without needing the instance itself.
    """
    def __init__(self):
        self._data = None

    @property
    def data(self) -> pd.DataFrame:
        """Returns a copy of the result data to prevent accidental modification."""
        return self._data.copy() if self._data is not None else pd.DataFrame()

    @data.setter
    def data(self, value: pd.DataFrame):
        """Sets the result data stored as a private property"""
        self._data = value

    def __len__(self):
        """Returns the length of the result data"""
        return 0 if self._data is None else len(self._data)


class PrevalenceResult(CalculatorResult):
    """
    Class to store the results of a `PrevalenceCalculator`.

    Parameters:
        stratified_by (list[str]): List of columns the result is stratified by - in order of strata level.
        adjusted_for (list[str]): Optional list of columns the prevalences were adjusted for (e.g., 'Cluster').
        n_distinct (list[str]): Optional list of columns distinct counts (per-strata) were generated for.
        denominator (str): Optional column indicating the denominator stratum.
    """
    def __init__(self, stratified_by: List[str], adjusted_for: List[str] = None, n_distinct: List[str] = None,
                 denominator: str = None):
        super().__init__()
        self.stratified_by: List[str] = stratified_by
        self.adjusted_for: Optional[List[str]] = adjusted_for
        self.n_distinct: Optional[List[str]] = n_distinct
        self.denominator: Optional[str] = denominator

    @classmethod
    def from_calculator(cls, calculator: 'PrevalenceCalculator') -> 'PrevalenceResult':
        """Sets up a `PrevalenceResult` instance using a `PrevalenceCalculator` instance."""
        return cls(calculator.stratify_by, calculator.adjust_for, calculator.n_distinct, calculator.denominator)
    

class Calculator(ABC):
    def __init__(self, model: 'ModelResult' = None):
        self.model = model

    @abstractmethod
    def calculate(self, dataset: Dataset): pass


class PrevalenceCalculator(Calculator):
    """
    A class to calculate raw and adjusted prevalence statistics from a dataframe.
    This class is designed to replicate the logic of the R 'prevalence' function
    using pandas for efficient data manipulation.
    
     Parameters:
        stratify_by (list[str]): List of columns to stratify the analysis by.
        adjust_for (list[str]): Optional list of columns for adjustment (e.g., 'Cluster').
        n_distinct (list[str]): Optional list of columns to calculate distinct counts for.
        denominator (str): Optional column to use as the primary grouping for denominators.
                     If None, the first column in stratify_by is used.
    """
    def __init__(self, stratify_by: List[str], adjust_for: List[str] = None, n_distinct: List[str] = None,
                 denominator: str = None):
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
        """
        Calculates prevalence statistics based on the provided dataset.

        This method computes raw and optionally adjusted counts, proportions,
        standard errors, and 95% confidence intervals (using the Wilson score
        interval method) for specified strata. It also calculates distinct
        counts for designated columns and ranks the prevalences within
        denominator groups.

        Parameters:
            dataset (Dataset | pd.Dataframe): Dataset to calculate prevalence of
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


# class DiversityCalculator(Calculator):
#     pass
#
#
# class AlphaDiversityCalculator(DiversityCalculator):
#     def __init__(self, finite: bool = False,
#                  metrics: list[_ALPHA_DIVERSITY_METRICS] = get_args(_ALPHA_DIVERSITY_METRICS)):
#         super().__init__()
#         self.finite: bool = finite
#         self.metrics: list[_ALPHA_DIVERSITY_METRICS] = metrics
#
#     def calculate(self, dataset: Dataset):
#
#
#     def simpson(self, counts, finite=False):
#         return 1 - dominance(counts, finite=finite)
#
#     def simpson_e(self, counts):
#         return 1 / (counts.size * dominance(counts))
#
#
# class BetaDiversityCalculator(DiversityCalculator):
#     pass



# Functions ------------------------------------------------------------------------------------------------------------
def _wilson_score_interval(counts: np.ndarray, denominators: np.ndarray) -> tuple[float, float, float, float]:
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
    """Calculates quantile-based credible intervals."""
    lower_q = (1 - prob) / 2
    upper_q = 1 - lower_q
    # Quantiles are calculated over the flattened chain/draw dimensions
    flat_samples = samples.reshape(-1, samples.shape[-1]) if samples.ndim > 1 else samples.flatten()
    return np.quantile(flat_samples, lower_q, axis=0), np.quantile(flat_samples, upper_q, axis=0)
