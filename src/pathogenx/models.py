"""
Object-Oriented Python port of a Bayesian mixed-effects model using Numpyro.
"""
from abc import ABC, abstractmethod
from typing import Literal
from pathlib import Path
import pickle

import pandas as pd
import numpy as np

import jax
import jax.numpy as jnp
from jax.scipy.special import expit
import numpyro
import numpyro.distributions as dist
from numpyro.infer import MCMC, NUTS, SVI, ELBO, Predictive
from numpyro.infer.autoguide import AutoNormal
import numpyro.optim as optim

from pathogenx.io import Dataset


# Classes --------------------------------------------------------------------------------------------------------------
class ModelResult(ABC):
    def __init__(self): self._data: pd.DataFrame = None

    @abstractmethod
    @property
    def data(self) -> pd.DataFrame: return self._data.copy()

    @abstractmethod
    @data.setter
    def data(self, value: pd.DataFrame): self._data = value

    @abstractmethod
    def save(self, filepath: Path): pass

    @classmethod
    @abstractmethod
    def from_file(cls, filepath: Path): pass


class BayesianMixedModelResult(ModelResult):
    def __init__(self):
        super().__init__()


class Model(ABC):
    def __init__(self): self.results = None

    @abstractmethod
    def fit(self, dataset: Dataset) -> ModelResult: pass


class BayesianMixedModel(Model):
    """
    A wrapper for a Bayesian mixed-effects logistic regression model using Numpyro.

    This class is designed to be data-agnostic. You provide the column
    names for the phenotype (dependent variable) and the grouping factor
    (random effect) during initialization.
    """

    def __init__(self, phenotype_col: str, group_col: str, seed: int = 0):
        super().__init__()
        self.phenotype_col = phenotype_col
        self.group_col = group_col
        self.rng_key = jax.random.PRNGKey(seed)

    @classmethod
    def from_file(cls, filepath: Path):
        """
        Loads pre-computed model results from a pickle file.
        """
        if not filepath.exists() or not filepath.is_file():
            raise FileNotFoundError(f"No file found at {filepath}")
        return pickle.loads(filepath.read_bytes())

    def _fit_mcmc(self, draws: int = 2000, chains: int = 4, warmup: int = 1000, **kwargs):
        """
        Fits the model using full MCMC sampling (NUTS sampler).

        Args:
            draws (int): Number of samples to draw (per chain).
            chains (int): Number of chains to run.
            warmup (int): Number of warmup steps (per chain).
            **kwargs: Additional arguments passed to `numpyro.infer.MCMC()`.
        """
        self.rng_key, fit_key = jax.random.split(self.rng_key)
        kernel = NUTS(model)
        mcmc = MCMC(kernel, num_warmup=warmup, num_samples=draws, num_chains=chains, **kwargs)
        mcmc.run(fit_key, group_idx=group_codes_jax, y_obs=phenotype_data_jax)
        # Get samples with chains as the first dimension
        # This returns a dict: {var_name: jnp.array[chain, draw, ...]}
        self.results = mcmc.get_samples(group_by_chain=True)
        return self.results

    def _fit_vi(self, n_samples: int = 1000, n_steps: int = 20000, **kwargs):
        """
        Fits the model using Variational Inference (ADVI).

        Args:
            n_samples (int): How many samples to draw from the approximation.
            n_steps (int): Number of optimization steps.
            **kwargs: Additional arguments passed to `numpyro.infer.SVI()`.
        """
        self.rng_key, fit_key, sample_key = jax.random.split(self.rng_key, 3)
        guide = AutoNormal(model)
        optimizer = optim.Adam(step_size=0.01)
        svi = SVI(model, guide, optimizer, loss=ELBO())
        result = svi.run(fit_key, n_steps, group_idx=group_codes_jax, y_obs=phenotype_data_jax, **kwargs)

        # Sample from the approximate posterior
        predictive = Predictive(guide, params=result.params, num_samples=n_samples)
        posterior_samples = predictive(sample_key, group_idx=group_codes_jax)

        # Standardize results to {var: [chain, draw, ...]} format (1 chain)
        # This makes it compatible with our plotting/summary functions
        results = {
            k: v[np.newaxis, ...] for k, v in posterior_samples.items()
        }
        return results

    def fit(self, dataset: Dataset, method: Literal['mcmc', 'vi'] = 'mcmc', **kwargs) -> BayesianMixedModelResult:
        data = dataset.data
        # Factorize data for Numpyro model
        group_codes, group_labels = pd.factorize(data[self.group_col])
        n_groups = len(group_labels)
        # Convert data to JAX arrays for the model
        phenotype_data_jax = jnp.array(data[self.phenotype_col].values)
        group_codes_jax = jnp.array(group_codes)
        # --- Priors ---
        # Global intercept
        intercept = numpyro.sample("Intercept", dist.Normal(0, 1.0))
        # Priors for the random effects (non-centered parameterization)
        sigma_group = numpyro.sample("sigma_group", dist.HalfNormal(1.0))
        # Numpyro needs sample_shape tuple
        offset_group = numpyro.sample("offset_group", dist.Normal(0, 1.0), sample_shape=(n_groups,))
        effect_group = numpyro.deterministic("effect_group", offset_group * sigma_group)
        # --- Linear Model (Logit-Link) ---
        # Map the random effects to the corresponding observations
        mu = intercept + effect_group[group_idx]
        # --- Likelihood ---
        # Use a plate for vectorized operations
        with numpyro.plate("data", len(group_idx)):
            numpyro.sample("y_obs", dist.Bernoulli(logits=mu), obs=y_obs)

        if method == 'mcmc':
            return self._fit_mcmc(**kwargs)
        elif method == 'vi':
            return self._fit_vi(**kwargs)
        else:
            raise NotImplementedError(f"Method {method} not implemented.")


# def _model_pipeline(self):
#     if self.results is None:
#         raise BayesianMixedModelError("Model has not been fitted. Call a fit method first.")
#
#         # --- Adjusted (Group-Specific) Prevalence ---
#         # Add intercept to each group effect. Shape: (chains, draws, n_groups)
#         # Intercept shape: (chains, draws) -> add new axis for broadcasting
#     logits_adj = self.results["Intercept"][..., np.newaxis] + self.results["effect_group"]
#     prop_adj = expit(logits_adj)  # Apply inverse-logit
#
#     # --- Raw (Global) Prevalence ---
#     logits_raw = self.results["Intercept"]
#     prop_raw = expit(logits_raw)
#
#     # --- Summarize and Format Results ---
#     # Calculate mean and credible intervals for adjusted prevalences
#     mean_adj = np.mean(prop_adj, axis=(0, 1))
#     lower_adj, upper_adj = _calculate_ci(prop_adj, ci_prob)
#
#     # Create the results DataFrame
#     result_df = pd.DataFrame({
#         self.group_col: self.group_labels,
#         'prop.adj': mean_adj,
#         'lower.adj': lower_adj,
#         'upper.adj': upper_adj,
#     })
#
#     # Add raw prevalence to each row for consistency
#     result_df['prop.raw'] = np.mean(prop_raw)
#     lower_raw, upper_raw = _calculate_ci(prop_raw, ci_prob)
#     result_df['lower.raw'] = lower_raw
#     result_df['upper.raw'] = upper_raw
#
#     # Create a PrevalenceResult object to hold the data
#     result_obj = PrevalenceResult(
#         stratified_by=[self.group_col],
#         adjusted_for=[self.phenotype_col]  # Use phenotype as the "adjustment" concept
#     )
#     result_obj.data = result_df
#     return result_obj
