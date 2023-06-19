TAS = analysis_module([cyto, speck])
from scipy import stats
import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.compare_functions as bs_compare
import bootstrapped.stats_functions as bs_stats

self.modules["TS_Cyto"].Max_Normalized_Measurement
module = TAS.modules["TS_Cyto"]
module = TAS.modules["TS_Speck"]
module.data["Treatment"].unique()
module.data[(module.data["Treatment"] == "ATP") & (module.data["Time (hrs)"] == 0)]

from modules.plots_module import plotting_module

plotting_module.plot_lineplot(module.data, "Measurement")

plotting_module.plot_pointplot(ratio_subset, "Ratio_IL1b_IL18")
ratio_info = TAS.modules["TS_Cyto"].ratio_data["IL18:IL1b"]
ratio_subset = ratio_info[ratio_info["Treatment"] == "ATP"]
ratio_time_slice = ratio_subset[ratio_subset["Time (hrs)"] == 10]

scipy.stats.sem(ratio_time_slice["Ratio_IL18_IL1b"])
ratio_name = "IL18:IL1b"
module = TAS.modules["TS_Speck"]
TAS.aggregate_time_comparisons(module)
TAS.plot_ratio(module)
array1 = array1["Measurement"]

array2 = array2["Measurement"]
dfcol2 = array2
equalize_series_lengths(value, mvalue)

array1 = module.data[
    (module.data["Treatment"] == "ATP") & (module.data["Time (hrs)"] == 0)
]["Measurement"]
array2 = module.data[
    (module.data["Treatment"] == "ATP") & (module.data["Time (hrs)"] == 1)
]["Measurement"]

bootstrap_t_test(array1, array2)

val = array1
grp = array2



def equalize_series_lengths(series1, series2):
    if len(series1) < len(series2):
        series1 = series1.sample(len(series2), replace=True)
    elif len(series2) < len(series1):
        series2 = series2.sample(len(series1), replace=True)
    return series1, series2


def bootstrap_t_test(array1, array2, n_bootstrap=1000):
    # Combine the two arrays
    val = np.concatenate((array1, array2))
    grp = np.array([1] * len(array1) + [2] * len(array2))

    # Calculate the observed t statistic
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


def aggregate_time_comparisons(
    self, module, output_directory: Path = Path("reports"), time1: float = 0
):
    # Get a list of all treatments, including modified treatments
    treatments = [treatment for treatment in module.data["Treatment"].unique()]
    analytes = [analyte for analyte in module.data["Analyte"].unique()]
    # Initialize an empty dictionary to store the significant times for each treatment
    significant_times = {}
    filename = output_directory / f"{module.name}_significant_times.xlsx"

    for analyte in analytes:
        data = module.data[module.data["Analyte"] == analyte]
        # Iterate over each treatment
        for treatment in treatments:
            # Calculate the time range with a significant difference from zero count
            table = self.time_point_comparison(treatment, analyte, module)
            significant_range = (
                table[table["P-Value"] < 0.05].Time.min(),
                table[table["P-Value"] < 0.05].Time.max(),
            )

            # Store the significant times in the dictionary
            significant_times[treatment] = significant_range

            # Write the significant times to the report
            print(
                f"{treatment} is significant between {significant_range[0]} and {significant_range[1]}."
            )

            sheet_name = treatment
            if len(analytes) > 1:
                sheet_name = str(analyte + "_" + treatment)
            # Write the significant times to an Excel file
            # If filename doesn't exist, create it
            if not filename.exists():
                table.to_excel(filename, sheet_name=sheet_name, index=True)
            with pd.ExcelWriter(
                filename, mode="a", if_sheet_exists="replace"
            ) as writer:
                table.to_excel(writer, sheet_name=sheet_name, index=True)


def calculate_ratio(self, module, measurement_type="Measurement"):
    if len(module.data["Analyte"].unique()) < 2:
        raise ValueError("Model does not have two analytes.")
    df2, df1 = self.split_analytes(module).values()
    analyte1, analyte2 = self.split_analytes(module).keys()
    merged_data = pd.merge(
        df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
    )
    merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
        merged_data[f"{measurement_type}_x"] / merged_data[f"{measurement_type}_y"]
    )
    merged_data[f"Inv_Ratio_{analyte1}_{analyte2}"] = (
        merged_data[f"{measurement_type}_y"] / merged_data[f"{measurement_type}_x"]
    )
    ATP_Sample = merged_data[
        (merged_data["Treatment"] == "ATP") & (merged_data["Time (hrs)"] == 23)
    ]
    MSU_Sample = merged_data[
        (merged_data["Treatment"] == "MSU") & (merged_data["Time (hrs)"] == 23)
    ]

    bootstrap_t_test(
        ATP_Sample[f"Ratio_{analyte1}_{analyte2}"],
        MSU_Sample[f"Ratio_{analyte1}_{analyte2}"],
    )
    bootstrap_t_test(
        ATP_Sample[f"Inv_Ratio_{analyte1}_{analyte2}"],
        MSU_Sample[f"Inv_Ratio_{analyte1}_{analyte2}"],
    )


dfcol1 = ATP_Sample["Measurement_x"]
dfcol2 = MSU_Sample["Measurement_x"]

from scipy.stats import f_oneway, bootstrap
import numpy as np
import bootstrapped.bootstrap as bs
import bootstrapped.compare_functions as bs_compare
import bootstrapped.stats_functions as bs_stats

# Your two lists of numbers
list1 = np.array([1, 2, 3, 4, 5])
list2 = np.array([2, 3, 4, 5, 6])


def new_bootsrap(dfcol1, dfcol2):
    # merge arrays
    diff = np.array(dfcol1) - np.array(dfcol2)
    bs_result = bs.bootstrap(diff, stat_func=bs_stats.mean)

    # Print the bootstrap estimate of the mean difference and the 95% confidence interval
    print(f"Bootstrap Mean Difference: {bs_result.value:.3f}")
    print(
        f"95% Confidence Interval: ({bs_result.lower_bound:.3f}, {bs_result.upper_bound:.3f})"
    )

    concat = np.concatenate((dfcol1, dfcol2))
    num_greater = 0.0
    num_permutations = 10000
    diff_observed = np.mean(dfcol2) - np.mean(dfcol1)

    for _ in range(num_permutations):
        perm = np.random.permutation(concat)
        perm_list1 = perm[: len(list1)]
        perm_list2 = perm[len(list1) :]
        diff_perm = np.mean(perm_list2) - np.mean(perm_list1)
        if diff_perm > diff_observed:
            num_greater += 1

    p_value = num_greater / num_permutations
    confidence_interval = (
        round(bs_result.lower_bound, 3),
        round(bs_result.upper_bound, 3),
    )
    return p_value, confidence_interval
