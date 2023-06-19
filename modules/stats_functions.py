import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
from scipy import stats
import scipy


def polynomial(x, *coeffs):
    return np.polyval(coeffs, x)


def cosine_similarity(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return dot_product / (norm_v1 * norm_v2)


def bootstrap_t_test(array1, array2, n_bootstrap=1000):
    # Combine the two arrays
    val = np.concatenate((array1, array2))
    grp = np.array([1] * len(array1) + [2] * len(array2))

    # Calculate the observed t statistic
    ## The t-statistic is calculated as the difference between the group means divided by the standard error of the difference.

    observed_t_stat = stats.ttest_ind(
        val[grp == 1], val[grp == 2], equal_var=True
    ).statistic

    # Perform bootstrapping
    t_values = np.zeros(n_bootstrap)
    n1 = len(array1)
    n2 = len(array2)
    for j in range(n_bootstrap):
        sample = np.random.choice(val, size=n1 + n2, replace=True)
        group1 = sample[:n1]
        group2 = sample[n1 : n1 + n2]
        if np.std(group1) == 0 or np.std(group2) == 0:
            t_values[j] = np.nan
        else:
            t_values[j] = stats.ttest_ind(group1, group2, equal_var=True).statistic

    # Calculate p-value
    p_value = np.nanmean(np.abs(t_values) >= np.abs(observed_t_stat))

    # Confidence interval info
    dfcol1 = val
    dfcol2 = grp
    diffs = []
    for i in range(n_bootstrap):
        dfcol1, dfcol2 = equalize_series_lengths(dfcol1, dfcol2)
        diff = np.array(dfcol1) - np.array(dfcol2)
        diffs.append((diff))
    diffs_mean = np.mean(diffs, axis=0)
    bs_result = bs.bootstrap(diffs_mean, stat_func=bs_stats.mean)

    confidence_interval = (
        round(bs_result.lower_bound, 3),
        round(bs_result.upper_bound, 3),
    )

    return p_value, confidence_interval


def equalize_series_lengths(series1, series2):
    if len(series1) < len(series2):
        series1 = series1.sample(len(series2), replace=True)
    elif len(series2) < len(series1):
        series2 = series2.sample(len(series1), replace=True)
    return series1, series2


def mean_confidence_interval(data, confidence=0.95):
    a = 1.0 * np.array(data)
    n = len(a)
    mean, se = np.mean(a), scipy.stats.sem(a)
    h = se * scipy.stats.t.ppf((1 + confidence) / 2.0, n - 1)
    return mean, (mean - h, mean + h)


def mean_standard_error(data):
    a = 1.0 * np.array(data)
    n = len(a)
    mean, se = np.mean(a), scipy.stats.sem(a)
    return mean, se

