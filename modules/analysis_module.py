from typing import Type, List, Dict, Any
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from pathlib import Path
import scipy.stats
from itertools import combinations
from itertools import product
from scipy.stats import f_oneway
from statsmodels.stats.multicomp import pairwise_tukeyhsd


def polynomial(x, *coeffs):
    return np.polyval(coeffs, x)


def cosine_similarity(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return dot_product / (norm_v1 * norm_v2)


def bootstrap_t_test(array1, array2, n_bootstrap=1000, ci=0.95):
    # Calculate the observed mean difference
    observed_mean_diff = np.mean(array1) - np.mean(array2)

    # Combine the two arrays
    combined_array = np.concatenate((array1, array2))

    # Perform bootstrapping
    bootstrap_mean_diffs = []
    for _ in range(n_bootstrap):
        np.random.shuffle(combined_array)
        sample1 = combined_array[: len(array1)]
        sample2 = combined_array[len(array1) :]
        bootstrap_mean_diff = np.mean(sample1) - np.mean(sample2)
        bootstrap_mean_diffs.append(bootstrap_mean_diff)

    # Calculate p-value
    p_value = (
        np.sum(np.abs(bootstrap_mean_diffs) >= np.abs(observed_mean_diff)) / n_bootstrap
    )

    # Calculate confidence interval
    lower_bound = np.percentile(bootstrap_mean_diffs, (1 - ci) / 2 * 100)
    upper_bound = np.percentile(bootstrap_mean_diffs, (1 + ci) / 2 * 100)
    confidence_interval = (round(lower_bound, 3), round(upper_bound, 3))

    return p_value, confidence_interval


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


class arbitrary_max:
    def __init__(self):
        self.real_max = None
        self.real_time = None
        self.window = None
        self.time = None
        self.max = None


class placeholder:
    def __init__(self):
        pass


class analysis_module:
    def __init__(self, component_list: List = []):
        self.modules = {module.name: module for module in component_list}
        self.time_compare()

    def time_compare(self):
        """
        Compare time-dependent factors across modules.

        This method compares the "Normalized_Measurement" and "Change_Rate"
        factors across all modules in the instance. For each factor, it
        calculates the maximum value for each group of (Analyte, Treatment,
        Experimental_Replicate) in each module and adds a new DataFrame
        attribute to the corresponding module object with the maximum
        value for that factor.

        Args:
            self (TimeComparison): The instance to compare.

        Returns:
            None.
        """
        # Loop over all modules in the instance
        for module in self.modules:
            # Loop over both factors to compare
            for factor in ["Normalized_Measurement", "Change_Rate"]:
                # Group the module data by (Analyte, Treatment, Experimental_Replicate)
                # and calculate the maximum value of the factor for each group
                max_factor = (
                    self.modules[module]
                    .data.groupby(["Analyte", "Treatment", "Experimental_Replicate"])
                    .agg(Max_Factor=(factor, "max"))
                    .reset_index()
                )
                # Merge the maximum factor values back into the module data to identify
                # the rows with the maximum factor values for each group
                max_factor = self.modules[module].data.merge(
                    max_factor,
                    left_on=["Analyte", "Treatment", "Experimental_Replicate", factor],
                    right_on=[
                        "Analyte",
                        "Treatment",
                        "Experimental_Replicate",
                        "Max_Factor",
                    ],
                )
                # Sort the resulting DataFrame by (Analyte, Treatment, Experimental_Replicate)
                sorted_data = max_factor.sort_values(
                    by=["Analyte", "Treatment", "Experimental_Replicate"]
                )
                # Set the attribute of the module object to the sorted DataFrame
                setattr(self.modules[module], f"Max_{factor}", sorted_data)

    def check_modifier_impact(self, module=None):
        if module is None:
            module = self.modules[next(iter(self.modules))]
        # Check if a modifier exists, if not, return None
        if module.comp.modifier is None:
            return None
        # Check that the modifier is a string, if not, return None
        if not isinstance(module.comp.modifier, str):
            return None
        results = []
        for treatment in module.comp.treatments:
            result = {"Treatment": treatment, "Modifier": module.comp.modifier}
            for factor in ["Normalized_Measurement", "Change_Rate"]:
                # Get the data for the current treatment and factor
                data = getattr(module, f"Max_{factor}")
                value = data[data["Treatment"] == treatment][factor]
                mvalue = data[
                    data["Treatment"] == str(result["Modifier"] + "_" + treatment)
                ][factor]
                time = data[data["Treatment"] == treatment]["Time (hrs)"]
                mtime = data[
                    data["Treatment"] == str(result["Modifier"] + "_" + treatment)
                ]["Time (hrs)"]
                p, ci = bootstrap_t_test(value, mvalue)
                tp, tci = bootstrap_t_test(time, mtime)
                result[f"base_peak_{factor}"] = value.mean()
                result[f"mod_peak_{factor}"] = mvalue.mean()
                result[f"peak_{factor}_p"] = p
                result[f"peak_{factor}_ci"] = ci
                result[f"base_peak_{factor}_time"] = time.mean()
                result[f"mod_peak_{factor}_time"] = mtime.mean()
                result[f"peak_{factor}_time_p"] = tp
                result[f"peak_{factor}_time_ci"] = tci
            results.append(result)
        results_df = pd.DataFrame(results)
        return results_df

    def time_point_comparison(
        self,
        Treatment: str = None,
        Analyte: str = None,
        module=None,
        time1: float = 0,
    ):
        """
        Compare the mean speck formation of a treatment at time1 with the mean speck formation at all other time points.

        Args:
            treatment (str): The treatment to analyze.
            time1 (float, optional): The initial time point for comparison. Defaults to 0.

        Returns:
            pd.DataFrame: A DataFrame containing time points, mean speck formation, differences, and p-values.
        """
        if Treatment is None:
            Treatment = module.data["Treatment"].unique()[0]
        if Analyte is None:
            Analyte = module.data["Analyte"].unique()[0]
        data = module.data[
            (module.data["Treatment"] == Treatment)
            & (module.data["Analyte"] == Analyte)
        ]
        unique_timepoints = data["Time (hrs)"].unique()
        unique_timepoints = np.delete(
            unique_timepoints, np.where(unique_timepoints == time1)
        )
        results = pd.DataFrame(columns=["Time", "Mean", "Difference", "P-Value", "CI"])
        time1_data = data[data["Time (hrs)"] == time1]
        for time2 in unique_timepoints:
            time2_data = data[(data["Time (hrs)"] == time2)]
            p, ci = bootstrap_t_test(
                time1_data["Measurement"], time2_data["Measurement"]
            )
            new_row = pd.DataFrame(
                {
                    "Time": float(time2),
                    "Mean": time2_data["Measurement"].mean(),
                    "Difference": -(
                        time1_data["Measurement"].mean()
                        - time2_data["Measurement"].mean()
                    ),
                    "P-Value": p,
                    "CI": str(str(ci[0]) + ":" + str(ci[1])),
                },
                index=[0],
            )
            results = pd.concat([results, new_row], ignore_index=True)
            # results.set_index("Time", inplace=True)
        return results

    def create_summary_table(self, module=None):
        """
        Creates a summary table with treatment statistics including max measurement, max change rate, and time to max values.

        Returns:
            pd.DataFrame: A summary table with treatment statistics.
        """
        if module is None:
            module = self.modules[next(iter(self.modules))]

        columns = [
            "Treatment",
            "Max Measurement",
            "Max Measurement (Normalized)",
            "Max Change Rate",
            "Time to Max Measurement",
            "Time to Max Change Rate",
        ]
        # If there are more than 1 analytes, add analyte column to start of columns list
        if len(module.data["Analyte"].unique()) > 1:
            columns.insert(0, "Analyte")
        summary = pd.DataFrame(columns=columns)
        all_treatments = module.data["Treatment"].unique()
        for analyte in module.data["Analyte"].unique():
            for treatment in all_treatments:
                set = module.Max_Normalized_Measurement[
                    (module.Max_Normalized_Measurement["Treatment"] == treatment)
                    & (module.Max_Normalized_Measurement["Analyte"] == analyte)
                ]
                max_measurement = set["Measurement"].mean()
                max_measurement_normalized = set["Normalized_Measurement"].mean()
                time_to_max_measurement = set["Time (hrs)"].mean()
                tset = module.Max_Change_Rate[
                    (module.Max_Change_Rate["Treatment"] == treatment)
                    & (module.Max_Change_Rate["Analyte"] == analyte)
                ]
                max_change_rate = tset["Change_Rate"].mean()
                max_change_rate_normalized = tset["Normalized_Change_Rate"].mean()
                time_to_max_change_rate = tset["Time (hrs)"].mean()
                new_row = pd.DataFrame(
                    {
                        "Analyte": analyte,
                        "Treatment": treatment,
                        "Max Measurement": max_measurement,
                        "Max Measurement (Normalized)": max_measurement_normalized,
                        "Max Change Rate": max_change_rate,
                        "Max Change Rate (Normalized)": max_change_rate_normalized,
                        "Time to Max Measurement": time_to_max_measurement,
                        "Time to Max Change Rate": time_to_max_change_rate,
                    },
                    index=[0],
                )
                if len(module.data["Analyte"].unique()) == 1:
                    new_row.drop(columns="Analyte", inplace=True)
                summary = pd.concat([summary, new_row], ignore_index=True)
        return summary

    def plot_ratio(
        self,
        model,
        measurement_type="Measurement",
        normalize_start=False,
    ):
        ##If the model's data does not have two analytes, return an error
        if len(model.data["Analyte"].unique()) < 2:
            raise ValueError("Model does not have two analytes.")
        df1, df2 = self.split_analytes(model).values()
        analyte1, analyte2 = self.split_analytes(model).keys()
        merged_data = pd.merge(
            df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
        )
        merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
            merged_data[f"{measurement_type}_x"] / merged_data[f"{measurement_type}_y"]
        )
        # If normalize_start is True, normalize the ratio to the first time point for each treatment and experimental replicate
        if normalize_start:
            merged_data[f"Ratio_{analyte1}_{analyte2}"] = merged_data.groupby(
                ["Treatment", "Experimental_Replicate"]
            )[f"Ratio_{analyte1}_{analyte2}"].transform(lambda x: x / x.iloc[0])

        self.plot_lineplot(merged_data, y=f"Ratio_{analyte1}_{analyte2}", errorbar="se")
        return merged_data

    def split_analytes(self, module):
        """
        Split the data into separate dataframes for each analyte.

        Args:
            model (time_sequence): The time_sequence object to split.

        Returns:
            dict: A dictionary of dataframes for each analyte.
        """
        analytes = module.data["Analyte"].unique()
        split_data = {}
        for analyte in analytes:
            split_data[analyte] = module.data[model.data["Analyte"] == analyte]
        return split_data

    def plot_lineplot(
        self,
        dataframe: pd.DataFrame,
        y: str,
        x: str = "Time (hrs)",
        hue: str = "Treatment",
        errorbar: str = "se",
        filepath: str = None,
        show_p=False,
    ):
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.lineplot(data=dataframe, x=x, y=y, hue=hue, errorbar=errorbar, ax=ax)
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f"{y} Over Time")
        ax.legend()
        if show_p == True:
            ax = self.show_p(data=dataframe, y=y, separator=hue, ax=ax)
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()

    def plot_multimodel(
        self,
        modelA,
        modelB,
        measurement_type="Normalized_Measurement",
        output_directory=None,
    ):
        ##Confirm that ModelA and ModelB have the same treatments
        if modelA.comp.treatments != modelB.comp.treatments:
            raise ValueError("Models do not have the same treatments.")
        for treatment in modelA.comp.treatments:
            dfA = modelA.data[modelA.data["Treatment"] == treatment]
            dfB = modelB.data[modelB.data["Treatment"] == treatment]
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.set_ylabel("Normalized Cytokine Measurement")
            ax2 = ax.twinx()
            ax2.set_ylabel("Normalized Speck Formation")
            for model in [dfA, dfB]:
                ax_to_use = ax if model is dfA else ax2
                if len(model["Analyte"].unique()) > 1:
                    sns.lineplot(
                        data=model,
                        x="Time (hrs)",
                        y=measurement_type,
                        hue="Analyte",
                        ax=ax_to_use,
                        errorbar="se",
                    )
                else:
                    sns.lineplot(
                        data=model,
                        x="Time (hrs)",
                        y=measurement_type,
                        color="green",
                        ax=ax_to_use,
                        errorbar="se",
                    )
            plt.xlabel("Time (hrs)")
            # plt.ylabel("Normalized Value")
            plt.title(
                str(
                    "Normalized Measurements for IL1b, IL18, and Normalized Speck Formation - "
                    + treatment
                )
            )
            if output_directory is not None:
                filepath = output_directory / f"{treatment}_multimodel.png"
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()

    def plot_speck_lineplots(
        self,
        module,
        treatments: List = None,
        output_path: Path = None,
        show_p=False,
    ):
        """
        Creates a line plot of the normalized speck formation over time, separated by treatment type.

        """
        if treatments is None:
            treatments = module.comp.treatments

        data = module.data[
            module.data["Treatment"].isin(treatments)
        ]  # Filter the data to only include the specified treatments

        # Check if there is only one treatment and multiple analytes
        if len(treatments) == 1 and len(data.Analyte.unique()) > 1:
            hue = "Analyte"
        # Check if there is only one analyte and multiple treatments
        elif len(data.Analyte.unique()) == 1 and len(treatments) >= 1:
            hue = "Treatment"
        # Otherwise, raise an error
        else:
            raise ValueError(
                "Must specify either one treatment and multiple analytes, or one analyte and multiple treatments."
            )

        if len(data[hue].unique()) == 1:
            y = "Measurement"
        else:
            y = "Normalized_Measurement"

        if output_path is not None:
            filename = output_path / f"{module.name}_{'_'.join(treatments)}.png"
            self.plot_lineplot(data, y=y, hue=hue, filepath=filename, show_p=show_p)
        else:
            self.plot_lineplot(data, y=y, hue=hue, show_p=show_p)

    def show_p(self, data, y, separator, ax):
        # Get the current ylim values
        ymin, ymax = ax.get_ylim()
        # Set the new ylim with an increased upper limit
        ax.set_ylim(ymin, ymax * 1.25)
        # Get the treatments in the data
        treatments = data[separator].unique()
        p_values_text = ""
        for treatment in treatments:
            treatment_data = data[data[separator] == treatment]
            max_time = treatment_data.loc[treatment_data[y].idxmax(), "Time (hrs)"]

            # max_time_per_replicate = (
            #     treatment_data.groupby("Experimental_Replicate")[y]
            #     .idxmax()
            #     .reset_index()
            # )
            # # Merge the original treatment_data with max_time_per_replicate to get the values for each replicate at their max times
            # max_values_per_replicate = treatment_data.loc[max_time_per_replicate[y]]
            # # Filter the max_values_per_replicate DataFrame for the current treatment
            # treatment_max_values = max_values_per_replicate[
            #     max_values_per_replicate[separator] == treatment
            # ]

            time_zero_data = np.array(
                treatment_data[treatment_data["Time (hrs)"] == 0][y]
            )
            max_time_data = np.array(
                treatment_data[treatment_data["Time (hrs)"] == max_time][y]
            )
            # print(time_zero_data, max_time_data)
            p, ci = bootstrap_t_test(max_time_data, time_zero_data)
            # print(p)
            p_values_text += f"{treatment} (Max = {max_time} hrs):\n p = {p:.4f}\n"

        ax.text(
            0.02,
            0.98,
            p_values_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            horizontalalignment="left",
        )
        return ax

    def find_arbitrary_max(self, module):
        self.time_compare()
        max_values = module.Max_Normalized_Measurement
        max_dict = {}
        module.arbitrary_maximums = placeholder()
        for analyte in max_values["Analyte"].unique():
            analyte_data = max_values[max_values["Analyte"] == analyte]
            mean_max_times = analyte_data.groupby(["Treatment"]).mean()["Time (hrs)"]
            # Round the mean max times to the nearest measurement time
            mean_max_times = mean_max_times.apply(
                lambda x: analyte_data["Time (hrs)"].tolist()[
                    np.abs(np.array(analyte_data["Time (hrs)"].tolist()) - x).argmin()
                ]
            )
            for treatment, time in mean_max_times.items():
                print(treatment)
                time_table = self.time_point_comparison(
                    treatment, module=module, time1=time
                )
                max_window = time_table[time_table["P-Value"] > 0.05]
                arbitrary_max_time = time_table[time_table["P-Value"] > 0.05][
                    "Time"
                ].min()
                arb_obj = arbitrary_max()
                arb_obj.real_time = time
                arb_obj.real_max = analyte_data.groupby(["Treatment"]).mean()[
                    "Normalized_Measurement"
                ][treatment]
                arb_obj.time = arbitrary_max_time
                arb_obj.window = max_window
                arb_obj.max = module.data[
                    (module.data["Time (hrs)"] == arbitrary_max_time)
                    & (module.data["Analyte"] == analyte)
                    & (module.data["Treatment"] == treatment)
                ]["Normalized_Measurement"].mean()
                setattr(module.arbitrary_maximums, treatment, arb_obj)

    def compare_max_time_distances(self, moduleA, moduleB):
        # Get tlist of all treatments which appear in both moduleA and moduleB data
        def get_times_by_treatment_analyte(df):
            result = {}
            analytes = df["Analyte"].unique()

            for analyte in analytes:
                analyte_dict = {}
                for treatment in df["Treatment"].unique():
                    times = df[
                        (df["Treatment"] == treatment) & (df["Analyte"] == analyte)
                    ]["Time (hrs)"].tolist()
                    analyte_dict[treatment] = times
                result[analyte] = analyte_dict

            return result

        def compare_differences(results, treatment1, treatment2):
            mean_diff1 = results[treatment1]["mean_diff"]
            mean_diff2 = results[treatment2]["mean_diff"]
            ci_diff1 = results[treatment1]["ci_diff"]
            ci_diff2 = results[treatment2]["ci_diff"]

            mean_diff_diff = mean_diff1 - mean_diff2

            se_diff1 = (ci_diff1[1] - ci_diff1[0]) / (
                2
                * scipy.stats.t.ppf(
                    (1 + 0.95) / 2,
                    len(times_A[treatment]) + len(times_B[treatment]) - 2,
                )
            )
            se_diff2 = (ci_diff2[1] - ci_diff2[0]) / (
                2
                * scipy.stats.t.ppf(
                    (1 + 0.95) / 2,
                    len(times_A[treatment]) + len(times_B[treatment]) - 2,
                )
            )

            se_diff_diff = np.sqrt(se_diff1**2 + se_diff2**2)
            ci_diff_diff = (ci_diff1[0] - ci_diff2[1], ci_diff1[1] - ci_diff2[0])
            return mean_diff_diff, ci_diff_diff

        max_A = get_times_by_treatment_analyte(moduleA.Max_Normalized_Measurement)
        max_B = get_times_by_treatment_analyte(moduleB.Max_Normalized_Measurement)
        treatments = list(
            set(moduleA.data["Treatment"].unique()).intersection(
                set(moduleB.data["Treatment"].unique())
            )
        )
        results = {}
        if len(max_A) == 1:
            max_A = {moduleA.name: max_A[next(iter(max_A.keys()))]}
        if len(max_B) == 1:
            max_B = {moduleB.name: max_B[next(iter(max_B.keys()))]}
        for treatment in treatments:
            for analyte_A, times_A in max_A.items():
                name_A = analyte_A
                for analyte_B, times_B in max_B.items():
                    name_B = analyte_B
                    key = f"{name_A}-{name_B}_{treatment}"
                    mean_A, se_A = mean_standard_error(times_A[treatment])
                    mean_B, se_B = mean_standard_error(times_B[treatment])

                    mean_diff = mean_A - mean_B
                    se_diff = np.sqrt(se_A**2 + se_B**2)
                    ci_diff = scipy.stats.t.interval(
                        0.95,
                        len(times_A[treatment]) + len(times_B[treatment]) - 2,
                        loc=mean_diff,
                        scale=se_diff,
                    )
                    results[key] = {"mean_diff": mean_diff, "ci_diff": ci_diff}

        comparison_results = {}
        pairwise_combinations = list(product(list(max_A.keys()), list(max_B.keys())))
        for pair in pairwise_combinations:
            prefix = f"{pair[0]}-{pair[1]}"
            for treatment_combination in combinations(treatments, 2):
                option1 = f"{prefix}_{treatment_combination[0]}"
                option2 = f"{prefix}_{treatment_combination[1]}"
                mean_diff_diff, ci_diff_diff = compare_differences(
                    results, option1, option2
                )
                comparison_results[f"{option1}:{option2}"] = {
                    "mean_diff_diff": mean_diff_diff,
                    "ci_diff_diff": ci_diff_diff,
                    "significance": not (ci_diff_diff[0] <= 0 <= ci_diff_diff[1]),
                }
        return comparison_results

    # def compare_max_time_distances(self, moduleA, moduleB):
    #     # Get tlist of all treatments which appear in both moduleA and moduleB data
    #     def get_times_by_treatment_analyte(df):
    #         result = {}
    #         analytes = df["Analyte"].unique()

    #         for analyte in analytes:
    #             analyte_dict = {}
    #             for treatment in df["Treatment"].unique():
    #                 times = df[
    #                     (df["Treatment"] == treatment) & (df["Analyte"] == analyte)
    #                 ]["Time (hrs)"].tolist()
    #                 analyte_dict[treatment] = times
    #             result[analyte] = analyte_dict

    #         return result

    #     max_A = get_times_by_treatment_analyte(moduleA.Max_Normalized_Measurement)
    #     max_B = get_times_by_treatment_analyte(moduleB.Max_Normalized_Measurement)

    #     treatments = list(
    #         set(moduleA.data["Treatment"].unique()).intersection(
    #             set(moduleB.data["Treatment"].unique())
    #         )
    #     )
    #     results = {}
    #     if len(max_A) == 1:
    #         max_A = {moduleA.name: max_A[next(iter(max_A.keys()))]}
    #     if len(max_B) == 1:
    #         max_B = {moduleB.name: max_B[next(iter(max_B.keys()))]}
    #     for treatment in treatments:
    #         for analyte_A, times_A in max_A.items():
    #             name_A = analyte_A
    #             for analyte_B, times_B in max_B.items():
    #                 name_B = analyte_B
    #                 key = f"{name_A}-{name_B}_{treatment}"

    #                 # T-test
    #                 t_stat, p_val = scipy.stats.ttest_ind(
    #                     times_A[treatment], times_B[treatment]
    #                 )
    #                 # Mean difference with CI
    #                 mean_A, se_A = mean_standard_error(times_A[treatment])
    #                 mean_B, se_B = mean_standard_error(times_B[treatment])

    #                 mean_diff = mean_A - mean_B
    #                 se_diff = np.sqrt(se_A**2 + se_B**2)
    #                 ci_diff = scipy.stats.t.interval(
    #                     0.95,
    #                     len(times_A[treatment]) + len(times_B[treatment]) - 2,
    #                     loc=mean_diff,
    #                     scale=se_diff,
    #                 )
    #                 results[key] = {
    #                     "t_stat": t_stat,
    #                     "p_val": p_val,
    #                     "mean_diff": mean_diff,
    #                     "ci_diff": ci_diff,
    #                 }

    #     comparison_results = {}
    #     pairwise_combinations = list(product(list(max_A.keys()), list(max_B.keys())))
    #     for pair in pairwise_combinations:
    #         prefix = f"{pair[0]}-{pair[1]}"
    #         for treatment_combination in combinations(treatments, 2):
    #             option1 = f"{prefix}_{treatment_combination[0]}"
    #             option2 = f"{prefix}_{treatment_combination[1]}"
    #             t_stat_diff = results[option1]["t_stat"] - results[option2]["t_stat"]
    #             mean_diff_diff, ci_diff_diff = compare_differences(
    #                 results, option1, option2
    #             )
    #             comparison_results[f"{option1}:{option2}"] = {
    #                 "t_stat_diff": t_stat_diff,
    #                 "t_significance": min(
    #                     results[option1]["p_val"], results[option2]["p_val"]
    #                 )
    #                 < 0.05,  # change this to use adjusted p-value if multiple comparison correction is applied
    #                 "mean_diff_diff": mean_diff_diff,
    #                 "ci_diff_diff": ci_diff_diff,
    #                 "significance": not (ci_diff_diff[0] <= 0 <= ci_diff_diff[1]),
    #             }
    #     return comparison_results

    def compare_max_normalized_measurement_anova(self, df):
        anova_results = {}

        for analyte in df["Analyte"].unique():
            df_analyte = df[df["Analyte"] == analyte]
            groups = [
                df_analyte["Normalized_Measurement"][
                    df_analyte["Treatment"] == treatment
                ]
                for treatment in df_analyte["Treatment"].unique()
            ]
            f_val, p_val = f_oneway(*groups)
            anova_results[analyte] = {"f_val": f_val, "p_val": p_val}
            if p_val < 0.05:
                post_hoc = pairwise_tukeyhsd(
                    df_analyte["Normalized_Measurement"], df_analyte["Treatment"]
                )
                anova_results[analyte]["post_hoc"] = post_hoc.summary()

        return anova_results


# self.time_compare()
# module = self.modules["TS_Speck"]
# module.arbitrary_maximums.ATP.window
# moduleA=TAS.modules["TS_Speck"]
# moduleB=TAS.modules["TS_Cyto"]
