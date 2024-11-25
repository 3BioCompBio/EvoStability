
"""
- Some basic statistical functions (such as Pearson Correlatoin Coefficient to compare 2 arrays).

"""

# Imports ----------------------------------------------------------------------
import numpy as np
from scipy import stats as st


# 2D-Metrics -------------------------------------------------------------------
def get_MSE(arr1, arr2):
    """
    Mean Squared Error
    """
    return np.mean( (np.array(arr1) - np.array(arr2))**2 )

def get_MAE(arr1, arr2):
    """
    Mean Absolute Error
    """
    return np.mean(np.absolute( np.array(arr1) - np.array(arr2) ))

def get_MSigE(arr1, arr2):
    """
    Mean Signed Error
    """
    return np.mean(np.array(arr1) - np.array(arr2))

def get_RMSE(arr1, arr2):
    """
    Root Mean Squared Error
    """
    return np.sqrt(np.mean( (np.array(arr1) - np.array(arr2))**2 ))

def get_pearson(arr1, arr2):
    """
    Pearson Correlation: cov(arr1, arr2) / sqrt(var(arr1)*var(arr2))
    """
    return np.corrcoef(np.array(arr1), np.array(arr2))[0, 1]

def get_spearman(arr1, arr2):
    """
    Spearman correlation of the ranks
    Wiki: https://en.wikipedia.org/wiki/Spearman%27s_rank_correlation_coefficient
    """
    return st.spearmanr(arr1, arr2)[0]

def get_kendall(arr1, arr2):
    """
    Kendall Correlation
    (|ordered pairs| - |disordered pairs|) / |pairs|
    Wiki: https://en.wikipedia.org/wiki/Kendall_rank_correlation_coefficient
    """
    return st.kendalltau(arr1, arr2)[0]

def get_rank_deviation(arr1, arr2):
    """Return Rank Deviation (absolute difference in ranks normalized by the length of the array)."""
    l = len(arr1)
    rank1 = st.rankdata(arr1)
    rank2 = st.rankdata(arr2)
    rank_deviation = np.absolute(rank1 - rank2) / l
    return rank_deviation

def get_spearman_impact(arr1, arr2):
    """Return the impact of one point on the global Spearman Correlation (normalized by the standard deviation of the impact)"""
    l = len(arr1)
    sp  = get_spearman(arr1, arr2)
    sp_impact = []
    for i in range(l):
        arr1_without = arr1[:i] + arr1[i+1:]
        arr2_without = arr2[:i] + arr2[i+1:]
        sp_without = get_spearman(arr1_without, arr2_without)
        sp_impact.append(sp - sp_without)
    sp_impact = np.array(sp_impact)
    sp_impact = sp_impact / np.std(sp_impact)
    return sp_impact

def get_pearson_impact(arr1, arr2):
    """Return the impact of one point on the global Pearson Correlation (normalized by the standard deviation of the impact)"""
    l = len(arr1)
    pr  = get_pearson(arr1, arr2)
    pr_impact = []
    for i in range(l):
        arr1_without = arr1[:i] + arr1[i+1:]
        arr2_without = arr2[:i] + arr2[i+1:]
        pr_without = get_spearman(arr1_without, arr2_without)
        pr_impact.append(pr - pr_without)
    pr_impact = np.array(pr_impact)
    pr_impact = pr_impact / np.std(pr_impact)
    return pr_impact

# 1D-Metrics -------------------------------------------------------------------
def get_mean(arr):
    return np.mean(np.array(arr))

def get_var(arr):
    return np.var(np.array(arr))

def get_std(arr):
    return np.std(np.array(arr))

def get_absmean(arr):
    return np.mean(np.abs(np.array(arr)))

def get_rsmean(arr):
    return np.sqrt(np.mean(np.array(arr)**2))