import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from openpyxl import Workbook, load_workbook
import seaborn as sns
import scipy.stats as stats
from pathlib import Path
from typing import Type, List, Dict
from scipy.optimize import curve_fit
from scipy.integrate import quad, simps
import itertools
from scipy.stats import pearsonr
import bootstrapped.bootstrap as bs
import bootstrapped.stats_functions as bs_stats
from scipy.stats import t
import numpy as np
from scipy.optimize import curve_fit
from scipy.integrate import quad
from sklearn.metrics.pairwise import cosine_similarity
import matplotlib.pyplot as plt


def polynomial(x, *coeffs):
    return np.polyval(coeffs, x)


def cosine_similarity(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return dot_product / (norm_v1 * norm_v2)


class ts_analysis_module:
    """
    A class to analyze time sequence data from the time_sequence_module.

    Attributes:
        comp (obj): An instance of the time_sequence_module class.
        data (pd.DataFrame): A DataFrame containing the normalized data.
        max_growth_rate (pd.DataFrame): A DataFrame containing the maximum growth rate data.
        max_speck_time (pd.DataFrame): A DataFrame containing the maximum speck time data.
    """

    def __init__(self, module):
        """
        Initializes the ts_analysis_module with a time_sequence_module instance.

        Args:
            module (obj): An instance of the time_sequence_module class.
        """
        self.comp = module
        self.data = module.data
        self.max_growth_rate = self.compare_growth_rate()
        self.max_speck_time = self.compare_max_speck_time()

    def line_plot(self):
        """
        Creates a line plot of the normalized speck formation over time, separated by treatment type.

        Returns:
            tuple: A tuple containing the figure and axes objects of the plot.
        """
        fig, ax = plt.subplots()
        sns.lineplot(
            ax=ax,
            data=self.data,
            x="Time (hrs)",
            y="NormalizedSpeckFormation",
            hue="Treatment",
        )
        ax.set_title("Speck Formation by Treatment Type")
        ax.set_ylabel("Normalized Speck Formation")
        ax.set_xlabel("Time")
        return (fig, ax)

    def line_plot_selected(self, treatments: List[str], show_p=False):
        """
        Creates a line plot of the normalized speck formation over time, separated by treatment type.

        Returns:
            tuple: A tuple containing the figure and axes objects of the plot.
        """
        data = self.data[self.data["Treatment"].isin(treatments)]
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.lineplot(
            ax=ax,
            data=data,
            x="Time (hrs)",
            y="NormalizedSpeckFormation",
            hue="Treatment",
        )
        ax.set_title("Speck Formation by Treatment Type")
        ax.set_ylabel("Normalized Speck Formation")
        ax.set_xlabel("Time")

        # Get the current ylim values
        ymin, ymax = ax.get_ylim()
        # Set the new ylim with an increased upper limit
        ax.set_ylim(ymin, ymax * 1.25)

        if show_p:
            p_values_text = ""
            for treatment in treatments:
                treatment_data = data[data["Treatment"] == treatment]
                max_time = treatment_data.loc[
                    treatment_data["NormalizedSpeckFormation"].idxmax(), "Time (hrs)"
                ]

                max_time_per_replicate = (
                    treatment_data.groupby("Treatment_Replicate")[
                        "NormalizedSpeckFormation"
                    ]
                    .idxmax()
                    .reset_index()
                )
                # Merge the original treatment_data with max_time_per_replicate to get the values for each replicate at their max times
                max_values_per_replicate = treatment_data.loc[
                    max_time_per_replicate["NormalizedSpeckFormation"]
                ]
                # Filter the max_values_per_replicate DataFrame for the current treatment
                treatment_max_values = max_values_per_replicate[
                    max_values_per_replicate["Treatment"] == treatment
                ]

                time_zero_data = np.array(
                    treatment_data[treatment_data["Time (hrs)"] == 0][
                        "NormalizedSpeckFormation"
                    ]
                )
                max_time_data = np.array(
                    treatment_max_values["NormalizedSpeckFormation"]
                )
                p, ci = self.bootstrap_t_test(max_time_data, time_zero_data)
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

        return (fig, ax)

    def bootstrap_t_test(self, array1, array2, n_bootstrap=1000, ci=0.95):
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
            np.sum(np.abs(bootstrap_mean_diffs) >= np.abs(observed_mean_diff))
            / n_bootstrap
        )

        # Calculate confidence interval
        lower_bound = np.percentile(bootstrap_mean_diffs, (1 - ci) / 2 * 100)
        upper_bound = np.percentile(bootstrap_mean_diffs, (1 + ci) / 2 * 100)
        confidence_interval = (round(lower_bound, 3), round(upper_bound, 3))

        return p_value, confidence_interval

    def compare_max_speck_time(self):
        """
        Compares the time it takes for each treatment to reach maximum speck formation using the Mann-Whitney U test.

        Returns:
            pd.DataFrame: A DataFrame containing the maximum speck formation times for each treatment.
        """
        max_speck_time = (
            self.data.groupby(["Treatment", "Treatment_Replicate"])
            .agg(max_speck=("NormalizedSpeckFormation", "max"))
            .reset_index()
        )
        max_speck_time = self.data.merge(
            max_speck_time,
            left_on=["Treatment", "Treatment_Replicate", "NormalizedSpeckFormation"],
            right_on=["Treatment", "Treatment_Replicate", "max_speck"],
        )
        # self.max_speck_time = max_speck_time
        # Perform Mann-Whitney U tests to compare the time to reach the maximum speck formation between treatments
        print("Time to reach the maximum speck formation:")
        for i, t1 in enumerate(self.comp.comp.treatments):
            for t2 in self.comp.comp.treatments[i + 1 :]:
                self.compare_treatments(t1, t2, max_speck_time, "Time (hrs)")
            if self.comp.comp.modifier:
                self.compare_treatments(
                    t1,
                    str(self.comp.comp.modifier + "_" + t1),
                    max_speck_time,
                    "Time (hrs)",
                )
        return max_speck_time

    def compare_growth_rate(self):
        """
        Compares the maximum growth rate between treatments using the Mann-Whitney U test.

        Returns:
            pd.DataFrame: A DataFrame containing the maximum growth rates for each treatment.
        """

        # Get the maximum GrowthRate for each Treatment_Replicate
        max_growth_rate = (
            self.data.groupby(["Treatment", "Treatment_Replicate"])
            .agg(max_growth=("GrowthRate", "max"))
            .reset_index()
        )
        max_growth_rate = self.data.merge(
            max_growth_rate,
            left_on=["Treatment", "Treatment_Replicate", "GrowthRate"],
            right_on=["Treatment", "Treatment_Replicate", "max_growth"],
        )
        # Perform Mann-Whitney U tests to compare the GrowthRate between treatments
        print("GrowthRate comparison:")
        for i, t1 in enumerate(self.comp.comp.treatments):
            for t2 in self.comp.comp.treatments[i + 1 :]:
                self.compare_treatments(t1, t2, max_growth_rate, "GrowthRate")
            if self.comp.comp.modifier:
                self.compare_treatments(
                    t1,
                    str(self.comp.comp.modifier + "_" + t1),
                    max_growth_rate,
                    "GrowthRate",
                )
        return max_growth_rate

    def compare_treatments(self, treatment1, treatment2, df, metric):
        """
        Compares two treatments using the Mann-Whitney U test.

        Args:
            treatment1 (str): The first treatment to compare.
            treatment2 (str): The second treatment to compare.
            df (pd.DataFrame): A DataFrame containing the data to compare.
            metric (str): The metric to compare between treatments (e.g. 'GrowthRate', 'Time (hrs)').
        """
        data1 = df[df["Treatment"] == treatment1][metric]
        data2 = df[df["Treatment"] == treatment2][metric]

        u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative="two-sided")
        print(f"{treatment1} vs {treatment2} ({metric}): p-value = {p_value}")

    def check_modifier_impact(self):
        """
        Checks the impact of the modifier on treatments by comparing the maximum speck time and maximum growth rate.
        """
        if not self.comp.comp.modifier:
            print("No modifier specified.")
            return

        # Calculate maximum speck time and maximum growth rate
        max_speck_time = self.compare_max_speck_time()
        max_growth_rate = self.compare_growth_rate()

        results = []

        # Loop through the treatments and compare them with their modified counterparts
        for t in self.comp.comp.treatments:
            result = {"Treatment": t, "Modifier": self.comp.comp.modifier}

            # Max growth rate
            base_growth = max_growth_rate[max_growth_rate["Treatment"] == t][
                "GrowthRate"
            ]
            mod_growth = max_growth_rate[
                max_growth_rate["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["GrowthRate"]

            p, ci = self.bootstrap_t_test(base_growth, mod_growth)
            result["Max_growth_rate_p_value"] = p
            result["Max_growth_rate_CI"] = ci
            result["Max_growth_rate_power_difference"] = np.mean(mod_growth) / np.mean(
                base_growth
            )
            result["Max_growth_rate_mean_(base)"] = np.mean(base_growth)
            result[f"Max_growth_rate_mean_({self.comp.comp.modifier})"] = np.mean(
                mod_growth
            )

            # Time to max growth rate
            base_time_growth = max_growth_rate[max_growth_rate["Treatment"] == t][
                "Time (hrs)"
            ]
            mod_time_growth = max_growth_rate[
                max_growth_rate["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["Time (hrs)"]

            p, ci = self.bootstrap_t_test(base_time_growth, mod_time_growth)
            result["Time_max_growth_rate_p_value"] = p
            result["Time_max_growth_rate_CI"] = ci
            result["Time_max_growth_rate_power_difference"] = np.mean(
                mod_time_growth
            ) / np.mean(base_time_growth)
            result["Time_max_growth_rate_mean_(base)"] = np.mean(base_time_growth)
            result[f"Time_max_growth_rate_mean_({self.comp.comp.modifier})"] = np.mean(
                mod_time_growth
            )

            # Max speck time
            base_time = max_speck_time[max_speck_time["Treatment"] == t]["Time (hrs)"]
            mod_time = max_speck_time[
                max_speck_time["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["Time (hrs)"]

            p, ci = self.bootstrap_t_test(base_time, mod_time)
            result["Max_speck_time_p_value"] = p
            result["Max_speck_time_CI"] = ci
            result["Max_speck_time_power_difference"] = np.mean(mod_time) / np.mean(
                base_time
            )
            result["Max_speck_time_mean_(base)"] = np.mean(base_time)
            result[f"Max_speck_time_mean_({self.comp.comp.modifier})"] = np.mean(
                mod_time
            )

            # Maximum amount of specks
            base_specks = self.data[self.data["Treatment"] == t][
                "NormalizedSpeckFormation"
            ]
            mod_specks = self.data[
                self.data["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["NormalizedSpeckFormation"]

            p, ci = self.bootstrap_t_test(base_specks, mod_specks)
            result["Max_specks_p_value"] = p
            result["Max_specks_CI"] = ci
            result["Max_specks_power_difference"] = np.mean(mod_specks) / np.mean(
                base_specks
            )
            result["Max_specks_mean_(base)"] = np.mean(base_specks)
            result[f"Max_specks_mean_({self.comp.comp.modifier})"] = np.mean(mod_specks)

            results.append(result)

            # Convert the results list to a DataFrame
        results_df = pd.DataFrame(results)
        return results_df

    def create_summary_table(self):
        """
        Creates a summary table with treatment statistics including max specks, max growth rate, and time to max values.

        Returns:
            pd.DataFrame: A summary table with treatment statistics.
        """
        summary = pd.DataFrame(
            columns=[
                "Treatment",
                "Max Specks (Absolute)",
                "Max Specks (Normalized)",
                "Max Growth Rate",
                "Time to Max Specks",
                "Time to Max Growth Rate",
            ]
        )
        all_treatments = list(self.comp.comp.treatments) + [
            f"{self.comp.comp.modifier}_{treatment}"
            for treatment in self.comp.comp.treatments
        ]
        for treatment in all_treatments:
            max_specks_t_abs = self.max_speck_time[
                self.max_speck_time["Treatment"] == treatment
            ]["SpeckFormation"].mean()
            max_specks_t_norm = self.max_speck_time[
                self.max_speck_time["Treatment"] == treatment
            ]["NormalizedSpeckFormation"].mean()
            max_growth_rate_t = self.max_growth_rate[
                self.max_growth_rate["Treatment"] == treatment
            ]["GrowthRate"].mean()
            time_max_specks_t = self.max_speck_time[
                self.max_speck_time["Treatment"] == treatment
            ]["Time (hrs)"].mean()
            time_max_growth_rate_t = self.max_growth_rate[
                self.max_growth_rate["Treatment"] == treatment
            ]["Time (hrs)"].mean()

            summary = summary.append(
                {
                    "Treatment": treatment,
                    "Max Specks (Absolute)": max_specks_t_abs,
                    "Max Specks (Normalized)": max_specks_t_norm,
                    "Max Growth Rate": max_growth_rate_t,
                    "Time to Max Specks": time_max_specks_t,
                    "Time to Max Growth Rate": time_max_growth_rate_t,
                },
                ignore_index=True,
            )
        return summary

    # def aggregate_ctsc0t(self):
    #     all_treatments = list(self.comp.comp.treatments) + [
    #         f"{self.comp.comp.modifier}_{treatment}"
    #         for treatment in self.comp.comp.treatments
    #     ]
    #     for treatment in all_treatments:
    #         table = self.compare_time_sequence_count_to_zero_time(treatment)
    #         table[table["P-value"] < 0.05].index.min()
    #         table[table["P-value"] < 0.05].index.max()

    def aggregate_ctsc0t(self):
        # Get a list of all treatments, including modified treatments
        all_treatments = list(self.comp.comp.treatments) + [
            f"{self.comp.comp.modifier}_{treatment}"
            for treatment in self.comp.comp.treatments
        ]

        # Initialize an empty dictionary to store the significant times for each treatment
        significant_times = {}
        filename = "significant_times_by_treatment.xlsx"
        # Iterate over each treatment
        for treatment in all_treatments:
            # Calculate the time range with a significant difference from zero count
            table = self.compare_time_sequence_count_to_zero_time(treatment)
            significant_range = (
                table[table["P-value"] < 0.05].index.min(),
                table[table["P-value"] < 0.05].index.max(),
            )

            # Store the significant times in the dictionary
            significant_times[treatment] = significant_range

            # Write the significant times to the report
            print(
                f"{treatment} is significant between {significant_range[0]} and {significant_range[1]}."
            )

            # Write the significant times to an Excel file
            with pd.ExcelWriter(
                filename, mode="a", if_sheet_exists="replace"
            ) as writer:
                table[table["P-value"] < 0.05].to_excel(
                    writer, sheet_name=treatment, index=True
                )

    def compare_time_sequence_count_to_zero_time(self, treatment, time1=0):
        """
        Compare the mean speck formation of a treatment at time1 with the mean speck formation at all other time points.

        Args:
            treatment (str): The treatment to analyze.
            time1 (float, optional): The initial time point for comparison. Defaults to 0.

        Returns:
            pd.DataFrame: A DataFrame containing time points, mean speck formation, differences, and p-values.
        """
        # Get all unique time points after time1
        unique_timepoints = self.data[self.data["Time (hrs)"] > time1][
            "Time (hrs)"
        ].unique()

        # Create an empty DataFrame to store the results
        results = pd.DataFrame(columns=["Time", "Mean", "P-value"])

        # Filter the data for the treatment and time1
        data_time1 = self.data[
            (self.data["Treatment"] == treatment) & (self.data["Time (hrs)"] == time1)
        ]["NormalizedSpeckFormation"]
        mean_time1 = data_time1.mean()

        # Iterate through the unique time points and compare the means
        for time2 in unique_timepoints:
            data_time2 = self.data[
                (self.data["Treatment"] == treatment)
                & (self.data["Time (hrs)"] == time2)
            ]["NormalizedSpeckFormation"]

            # Perform the Bootstrap test to compare the means
            p, ci = self.bootstrap_t_test(data_time1, data_time2)
            mean_time2 = data_time2.mean()

            # Calculate the difference between mean at time2 and mean at time1
            difference = mean_time2 - mean_time1

            # Add the results to the DataFrame
            new_row = pd.DataFrame(
                {
                    "Time": [time2],
                    "Mean": [mean_time2],
                    "P-value": [p],
                }
            )

            results = pd.concat([results, new_row], ignore_index=True)

        # Set the DataFrame index to the Time column
        results.set_index("Time", inplace=True)

        return results

    def plot_time_course_analysis(
        self, treatment1=None, treatment2=None, degree=3, points_only=False
    ):
        """
        Plot the time course analysis for one or two treatments.

        Args:
            treatment1 (str, optional): The first treatment to analyze. Defaults to None.
            treatment2 (str, optional): The second treatment to analyze. Defaults to None.
            degree (int, optional): Degree of the polynomial fit. Defaults to 3.
            points_only (bool, optional): If True, plot only the data points without the fitted curve. Defaults to False.
        """
        grouped_data = (
            self.data.groupby(["Treatment", "Time (hrs)"]).mean().reset_index()
        )

        if treatment1 is not None and treatment2 is None:
            treatment_list = [treatment1]
            if self.comp.comp.modifier:
                treatment_list.append(str(self.comp.comp.modifier + "_" + treatment1))
        elif treatment1 is not None and treatment2 is not None:
            treatment_list = [treatment1, treatment2]
        else:
            treatment_list = grouped_data["Treatment"].unique()

        fit_params = {}

        for treatment in treatment_list:
            treatment_data = grouped_data[grouped_data["Treatment"] == treatment]
            x_data = treatment_data["Time (hrs)"].values
            y_data = treatment_data["NormalizedSpeckFormation"].values

            coeffs, _ = curve_fit(polynomial, x_data, y_data, p0=np.ones(degree + 1))
            fit_params[treatment] = coeffs

            x_plot = np.linspace(min(x_data), max(x_data), 1000)
            y_plot = polynomial(x_plot, *coeffs)

            plt.plot(x_data, y_data, "o", label=treatment)
            if not points_only:
                plt.plot(x_plot, y_plot, label=f"{treatment} fit")

        plt.xlabel("Time (hrs)")
        plt.ylabel("Normalized Speck Formation")
        plt.legend()
        plt.show()

    def quantify_curve_difference(self, treatment1, treatment2, degree=3):
        """
        Quantify the difference between the time course curves of two treatments.

        Args:
            treatment1 (str): The first treatment to analyze.
            treatment2 (str): The second treatment to analyze.
            degree (int, optional): Degree of the polynomial fit. Defaults to 3.
        """
        grouped_data = (
            self.data.groupby(["Treatment", "Time (hrs)"]).mean().reset_index()
        )

        def fitted_curve(treatment):
            treatment_data = grouped_data[grouped_data["Treatment"] == treatment]
            x_data = treatment_data["Time (hrs)"].values
            y_data = treatment_data["NormalizedSpeckFormation"].values

            coeffs, _ = curve_fit(polynomial, x_data, y_data, p0=np.ones(degree + 1))

            return lambda x: polynomial(x, *coeffs)

        curve1 = fitted_curve(treatment1)
        curve2 = fitted_curve(treatment2)

        min_time = min(grouped_data["Time (hrs)"])
        max_time = max(grouped_data["Time (hrs)"])

        # Area between the curves
        area, _ = quad(lambda x: np.abs(curve1(x) - curve2(x)), min_time, max_time)
        print(f"Area between the curves: {area}")

        # RMSE
        x_common = np.linspace(min_time, max_time, 1000)
        y1_common = curve1(x_common)
        y2_common = curve2(x_common)

        rmse = np.sqrt(np.mean((y1_common - y2_common) ** 2))
        print(f"Root Mean Squared Error: {rmse}")
        # Normalize the curves
        y1_mean = np.mean(y1_common)
        y1_std = np.std(y1_common)
        y1_normalized = (y1_common - y1_mean) / y1_std

        y2_mean = np.mean(y2_common)
        y2_std = np.std(y2_common)
        y2_normalized = (y2_common - y2_mean) / y2_std

        # Calculate cosine similarity
        cos_sim = cosine_similarity(y1_normalized, y2_normalized)
        print(f"Cosine similarity: {cos_sim}")

    def plot_normalized_speck_formation(df, treatment1, treatment2):
        # Filter the dataframe for the specified treatments
        df = self.data
        df_filtered = df[df["Treatment"].isin([treatment1, treatment2])]

        # Group by treatment and time
        grouped_data = (
            df_filtered.groupby(["Treatment", "Time (hrs)", "ExperimentalReplicate"])
            .mean()
            .reset_index()
        )
        x_max = grouped_data[grouped_data["Treatment"] == treatment1].max()
        # Normalize the speck formation
        for exp_rep in grouped_data["ExperimentalReplicate"].unique():
            grouped_data.loc[
                grouped_data["ExperimentalReplicate"] == exp_rep,
                "NormalizedSpeckFormation",
            ] = (
                grouped_data.loc[
                    grouped_data["ExperimentalReplicate"] == exp_rep, "SpeckFormation"
                ]
                .groupby("Treatment")["SpeckFormation"]
                .transform(lambda x: x / x.max())
            )
        # grouped_data['NormalizedSpeckFormation'] = grouped_data.groupby('Treatment')['SpeckFormation'].transform(lambda x: x / x.max())

        # Plot the normalized speck formation curves
        # for treatment in [treatment1, treatment2]:
        #     treatment_data = grouped_data[grouped_data['Treatment'] == treatment]
        # plt.plot(treatment_data['Time (hrs)'], treatment_data['NormalizedSpeckFormation'], label=treatment)
        # sns.lineplot(
        #     x="Time (hrs)",
        #     y="NormalizedSpeckFormation",
        #     hue="Treatment",
        #     data=grouped_data,
        #     errorbar="se",
        # )
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.lineplot(
            ax=ax,
            data=grouped_data,
            x="Time (hrs)",
            y="NormalizedSpeckFormation",
            hue="Treatment",
            errorbar="se",
        )
        ax.set_title("Speck Formation by Treatment Type")
        ax.set_ylabel("Normalized Speck Formation")
        ax.set_xlabel("Time")

        plt.show()

    # treatment1 = "MSU"
    # treatment2 = "MCC950_MSU"

    def plot_curve_difference(self, treatment1, treatment2, degree=3):
        """
        Quantify the difference between the time course curves of two treatments.

        Args:
            treatment1 (str): The first treatment to analyze.
            treatment2 (str): The second treatment to analyze.
            degree (int, optional): Degree of the polynomial fit. Defaults to 3.
        """
        grouped_data = (
            self.data.groupby(["Treatment", "Time (hrs)"]).mean().reset_index()
        )

        def fitted_curve(treatment):
            treatment_data = grouped_data[grouped_data["Treatment"] == treatment]
            x_data = treatment_data["Time (hrs)"].values
            y_data = treatment_data["SpeckFormation"].values

            coeffs, _ = curve_fit(polynomial, x_data, y_data, p0=np.ones(degree + 1))

            return lambda x: polynomial(x, *coeffs)

        curve1 = fitted_curve(treatment1)
        curve2 = fitted_curve(treatment2)

        min_time = min(grouped_data["Time (hrs)"])
        max_time = max(grouped_data["Time (hrs)"])

        # Area between the curves
        area, _ = quad(lambda x: np.abs(curve1(x) - curve2(x)), min_time, max_time)
        print(f"Area between the curves: {area}")

        # RMSE
        x_common = np.linspace(min_time, max_time, 1000)
        y1_common = curve1(x_common)
        y2_common = curve2(x_common)

        rmse = np.sqrt(np.mean((y1_common - y2_common) ** 2))
        print(f"Root Mean Squared Error: {rmse}")
        # Normalize the curves
        y1_mean = np.mean(y1_common)
        y1_std = np.std(y1_common)
        y1_normalized = (y1_common - y1_mean) / y1_std

        y2_mean = np.mean(y2_common)
        y2_std = np.std(y2_common)
        y2_normalized = (y2_common - y2_mean) / y2_std
        # Plot the standardized curves
        plt.plot(x_common, y1_normalized, label=treatment1)
        plt.plot(x_common, y2_normalized, label=treatment2)

        plt.xlabel("Time (hrs)")
        plt.ylabel("Standardized Normalized Speck Formation")
        plt.title("Standardized Normalized Speck Formation Curves")

        plt.legend(title="Treatment")
        plt.show()

        # Calculate cosine similarity
        cos_sim = cosine_similarity(y1_normalized, y2_normalized)
        print(f"Cosine similarity: {cos_sim}")

    def compare_replicate_curves(self, treatment):
        treatment_data = self.data[self.data["Treatment"] == treatment]
        replicates = treatment_data["Treatment_Replicate"].unique()

        # Plot the curves for each replicate
        for rep in replicates:
            rep_data = treatment_data[treatment_data["Treatment_Replicate"] == rep]
            plt.plot(
                rep_data["Time (hrs)"], rep_data["NormalizedSpeckFormation"], label=rep
            )

        plt.xlabel("Time (hrs)")
        plt.ylabel("Normalized Speck Formation")
        plt.legend()
        plt.title(f"Curves for Each Replicate of {treatment}")
        plt.show()

        # Calculate the ABC, RMSE, and Cosine Similarity for each pair of replicates
        comparison_results = []
        for rep1, rep2 in itertools.combinations(replicates, 2):
            rep1_data = treatment_data[treatment_data["Treatment_Replicate"] == rep1]
            rep2_data = treatment_data[treatment_data["Treatment_Replicate"] == rep2]

            # Calculate ABC
            y1 = rep1_data["NormalizedSpeckFormation"].values
            y2 = rep2_data["NormalizedSpeckFormation"].values
            x = rep1_data["Time (hrs)"].values
            abc = simps(np.abs(y1 - y2), x)

            # Calculate RMSE
            rmse = np.sqrt(np.mean((y1 - y2) ** 2))

            # Normalize the curves
            y1_mean = np.mean(y1)
            y1_std = np.std(y1)
            y1_normalized = (y1 - y1_mean) / y1_std

            y2_mean = np.mean(y2)
            y2_std = np.std(y2)
            y2_normalized = (y2 - y2_mean) / y2_std

            # Calculate cosine similarity
            cos_sim = cosine_similarity(y1_normalized, y2_normalized)

            comparison_results.append((rep1, rep2, abc, rmse, cos_sim))

        # Build the table with the comparison results
        results_table = pd.DataFrame(
            comparison_results,
            columns=["Replicate 1", "Replicate 2", "ABC", "RMSE", "Cosine Similarity"],
        )

        # Calculate the average ABC, RMSE, and Cosine Similarity
        avg_abc = results_table["ABC"].mean()
        avg_rmse = results_table["RMSE"].mean()
        avg_cos_sim = results_table["Cosine Similarity"].mean()

        print("Comparison results:")
        print(results_table)
        print(f"Average ABC: {avg_abc:.4f}")
        print(f"Average RMSE: {avg_rmse:.4f}")
        print(f"Average Cosine Similarity: {avg_cos_sim:.4f}")

    def compare_standardized_replicate_curves(self, treatment, shift_x=False):
        treatment_data = self.data[self.data["Treatment"] == treatment]
        replicates = treatment_data["Treatment_Replicate"].unique()
        adjusted_treatment_data = pd.DataFrame()

        def calc_max_time(rep_data):
            max_std = rep_data[
                rep_data["NormalizedSpeckFormation"]
                == np.max(rep_data["NormalizedSpeckFormation"])
            ]
            max_time = max_std["Time (hrs)"].values[0]
            return max_time

        max_times = []
        if shift_x == True:
            # Plot the standardized curves for each replicate
            for rep in replicates:
                rep_data = treatment_data[treatment_data["Treatment_Replicate"] == rep]
                max_times.append(calc_max_time(rep_data))
        for rep in replicates:
            rep_data = treatment_data[treatment_data["Treatment_Replicate"] == rep]
            if shift_x == True:
                print(calc_max_time(rep_data))
                max_time_diff = calc_max_time(rep_data) - min(max_times)
                print(max_time_diff)
                rep_data["Time (hrs)"] = rep_data["Time (hrs)"] - max_time_diff
            # print(rep_data)
            y = rep_data["NormalizedSpeckFormation"].values / np.max(
                rep_data["NormalizedSpeckFormation"].values
            )
            adjusted_treatment_data = pd.concat([adjusted_treatment_data, rep_data])
            plt.plot(rep_data["Time (hrs)"], y, label=rep)

        plt.xlabel("Time (hrs)")
        plt.ylabel("Standardized Normalized Speck Formation")
        plt.legend()
        plt.title(f"Standardized Curves for Each Replicate of {treatment}")
        # plt.show()
        # Calculate the ABC, RMSE, Cosine Similarity, and Estimated Difference for each pair of replicates
        comparison_results = []
        for rep1, rep2 in itertools.combinations(replicates, 2):
            rep1_data = adjusted_treatment_data[
                adjusted_treatment_data["Treatment_Replicate"] == rep1
            ]
            rep2_data = adjusted_treatment_data[
                adjusted_treatment_data["Treatment_Replicate"] == rep2
            ]

            x1 = rep1_data["Time (hrs)"].values
            x2 = rep2_data["Time (hrs)"].values
            y1_standardized = rep1_data["NormalizedSpeckFormation"].values / np.max(
                rep1_data["NormalizedSpeckFormation"].values
            )
            y2_standardized = rep2_data["NormalizedSpeckFormation"].values / np.max(
                rep2_data["NormalizedSpeckFormation"].values
            )

            x_common = np.union1d(x1, x2)
            y1_common = np.interp(x_common, x1, y1_standardized)
            y2_common = np.interp(x_common, x2, y2_standardized)

            # Calculate ABC
            abc = simps(np.abs(y1_common - y2_common), x_common)

            # Calculate RMSE
            rmse = np.sqrt(np.mean((y1_common - y2_common) ** 2))

            # Calculate Pearson correlation
            pearson_corr, _ = pearsonr(y1_common, y2_common)

            # Normalize the curves
            y1_mean = np.mean(y1_common)
            y1_std = np.std(y1_common)
            y1_normalized = (y1_common - y1_mean) / y1_std

            y2_mean = np.mean(y2_common)
            y2_std = np.std(y2_common)
            y2_normalized = (y2_common - y2_mean) / y2_std

            # Calculate cosine similarity
            cos_sim = np.dot(y1_normalized, y2_normalized) / (
                np.linalg.norm(y1_normalized) * np.linalg.norm(y2_normalized)
            )
            # Calculate Estimated Difference
            est_diff = rmse * (1 - cos_sim)

            comparison_results.append(
                (rep1, rep2, abc, rmse, cos_sim, pearson_corr, est_diff)
            )

        # Build the table with the comparison results
        results_table = pd.DataFrame(
            comparison_results,
            columns=[
                "Replicate 1",
                "Replicate 2",
                "ABC",
                "RMSE",
                "Cosine Similarity",
                "Pearson Correlation",
                "Estimated Difference",
            ],
        )

        # Calculate the average ABC, RMSE, Cosine Similarity, and Estimated Difference
        avg_abc = results_table["ABC"].mean()
        avg_rmse = results_table["RMSE"].mean()
        avg_cos_sim = results_table["Cosine Similarity"].mean()
        avg_pearson_corr = results_table["Pearson Correlation"].mean()
        avg_est_diff = results_table["Estimated Difference"].mean()
        print(results_table)
        return (plt, results_table)
