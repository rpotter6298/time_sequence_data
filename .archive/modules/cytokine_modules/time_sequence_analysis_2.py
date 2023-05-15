import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from pathlib import Path
from typing import Type, List, Dict
from scipy.optimize import curve_fit
from scipy.integrate import quad, simps
import itertools
from scipy.stats import pearsonr
from matplotlib.collections import LineCollection
import os


def polynomial(x, *coeffs):
    return np.polyval(coeffs, x)


def cosine_similarity(v1, v2):
    dot_product = np.dot(v1, v2)
    norm_v1 = np.linalg.norm(v1)
    norm_v2 = np.linalg.norm(v2)
    return dot_product / (norm_v1 * norm_v2)


class cyt_ts_analysis_module:
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
        self.data_split = self.split_analytes()
        self.groupby_cols = [
            n
            for n in ["Analyte", "Treatment", "Experimental_Replicate"]
            if n in self.data.columns
        ]
        self.max_rate = self.max_rate_of_change_table()
        self.measurement_max = self.max_measurement_table()

    def max_measurement_table(self):
        """
        Compares the time it takes for each treatment to reach maximum speck formation using the Mann-Whitney U test.

        Returns:
            pd.DataFrame: A DataFrame containing the maximum speck formation times for each treatment.
        """
        max_measurement = (
            self.data.groupby(self.groupby_cols)
            .agg(Max_Measurement=("Normalized_Measurement", "max"))
            .reset_index()
        )
        max_measurement = self.data.merge(
            max_measurement,
            left_on=self.groupby_cols + ["Normalized_Measurement"],
            right_on=self.groupby_cols + ["Max_Measurement"],
        )
        return max_measurement

    def max_rate_of_change_table(self):
        """
        Compares the growth rate between treatments using the Mann-Whitney U test.

        Returns:
            pd.DataFrame: A DataFrame containing the maximum growth rates for each treatment.
        """
        max_rate = (
            self.data.groupby(self.groupby_cols)
            .agg(Max_Rate=("Rate_of_Change", "max"))
            .reset_index()
        )
        max_rate = self.data.merge(
            max_rate,
            left_on=self.groupby_cols + ["Rate_of_Change"],
            right_on=self.groupby_cols + ["Max_Rate"],
        )
        return max_rate

    def split_analytes(self):
        """
        Splits the data by analyte.

        Returns:
            dict: A dictionary containing the dataframes split by analyte.
        """
        analytes = self.data["Analyte"].unique()
        analyte_dict = {}
        for analyte in analytes:
            analyte_dict[analyte] = self.data[self.data["Analyte"] == analyte]
        return analyte_dict

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

    def ratio_plot(self):
        df1, df2 = self.split_analytes().values()
        analyte1, analyte2 = self.split_analytes().keys()
        merged_data = pd.merge(
            df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
        )
        merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
            merged_data["Concentration_x"] / merged_data["Concentration_y"]
        )
        # Normalize ratio per replicate and treatment so that starting ratio is 1
        merged_data[f"Normalized_Ratio_{analyte1}_{analyte2}"] = merged_data[
            f"Ratio_{analyte1}_{analyte2}"
        ] / merged_data.groupby(["Treatment", "Experimental_Replicate"])[
            f"Ratio_{analyte1}_{analyte2}"
        ].transform(
            "first"
        )
        # Save the normalized ratio data to self.ratio_data
        self.ratio_data = merged_data
        # Plot the ratio for each treatment with ci determined by replicates
        fig, ax = plt.subplots(figsize=(10, 6))
        sns.lineplot(
            data=merged_data,
            x="Time (hrs)",
            y=f"Normalized_Ratio_{analyte1}_{analyte2}",
            hue="Treatment",
            ci="sd",
            ax=ax,
        )
        ax.set_xlabel("Time (hrs)")
        ax.set_ylabel(f"Normalized Ratio {analyte1} / {analyte2}")
        ax.set_title(f"Ratio of {analyte1} to {analyte2} for Each Treatment")
        ax.legend()
        plt.show()

    def split_line_plot(self):
        df1, df2 = self.split_analytes().values()

        # Set up colors for Analyte and Treatment
        treatment_colors = {
            "ATP": "green",
            "MSU": "purple",
            "Nigericin": "blue",
            "None": "black",
        }
        fig, ax = plt.subplots(figsize=(10, 6))

        # Plot lines for each analyte with a thicker line width
        sns.lineplot(
            data=df1,
            x="Time (hrs)",
            y="Normalized_Measurement",
            hue="Treatment",
            palette=treatment_colors,
            linewidth=4,
            legend=False,
        )
        sns.lineplot(
            data=df1,
            x="Time (hrs)",
            y="Normalized_Measurement",
            hue="Treatment",
            palette=["white"],
            linewidth=1,
            legend=False,
            ci=None,
        )

        # Plot lines for each treatment with a thinner line width, overlaying the analyte lines
        sns.lineplot(
            data=df2,
            x="Time (hrs)",
            y="Normalized_Measurement",
            hue="Treatment",
            palette=treatment_colors,
            linewidth=4,
        )
        sns.lineplot(
            data=df2,
            x="Time (hrs)",
            y="Normalized_Measurement",
            hue="Treatment",
            palette=["black"],
            linewidth=1,
            legend=False,
            ci=None,
        )
        ax.set_xlabel("Time (hrs)")
        ax.set_ylabel("Normalized Measurement")
        ax.set_title(
            "Normalized Measurement for Each Analyte/Treatment Combination with Confidence Interval Shading"
        )

        # Add a legend
        handles, labels = ax.get_legend_handles_labels()

        # Add white and black lines for IL18 and IL1b
        from matplotlib.lines import Line2D

        handles.append(Line2D([0], [0], color="white", linewidth=3))
        labels.append("IL18")
        handles.append(Line2D([0], [0], color="black", linewidth=3))
        labels.append("IL1b")

        ax.legend(handles=handles, labels=labels)

        plt.show()

    # def line_plot(self):
    #     """
    #     Creates a line plot of the normalized speck formation over time, separated by treatment type.

    #     Returns:
    #         tuple: A tuple containing the figure and axes objects of the plot.
    #     """
    #     fig, ax = plt.subplots()
    #     sns.lineplot(
    #         ax=ax,
    #         data=self.data,
    #         x="Time (hrs)",
    #         y="NormalizedSpeckFormation",
    #         hue="Treatment",
    #     )
    #     ax.set_title("Speck Formation by Treatment Type")
    #     ax.set_ylabel("Normalized Speck Formation")
    #     ax.set_xlabel("Time")
    #     return (fig, ax)

    def line_plot_selected(
        self,
        treatments: List[str],
        show_p=False,
        measurement_type="Normalized_Measurement",
    ):
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
            y="Normalized_Measurement",
            hue="Treatment",
        )
        ax.set_title(str(measurement_type + " by Treatment Type"))
        ax.set_ylabel(measurement_type)
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
                    treatment_data["Normalized_Measurement"].idxmax(), "Time (hrs)"
                ]
                time_zero_data = treatment_data[treatment_data["Time (hrs)"] == 0][
                    "Normalized_Measurement"
                ]
                max_time_data = treatment_data[
                    treatment_data["Time (hrs)"] == max_time
                ]["Normalized_Measurement"]
                u_stat, p_value = stats.ttest_rel(time_zero_data, max_time_data)
                p_values_text += (
                    f"{treatment} (Max = {max_time} hrs): p = {p_value:.4f}\n"
                )

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

    # def check_modifier_impact(self):
    #     """
    #     Checks the impact of the modifier on treatments by comparing the maximum speck time and maximum growth rate.
    #     """
    #     if not self.comp.comp.modifier:
    #         print("No modifier specified.")
    #         return

    #     # Calculate maximum speck time and maximum growth rate
    #     max_speck_time = self.compare_max_speck_time()
    #     max_growth_rate = self.compare_growth_rate()

    #     # Loop through the treatments and compare them with their modified counterparts
    #     for t in self.comp.comp.treatments:
    #         print(f"Comparing treatment {t} with {self.comp.comp.modifier}_{t}:")

    #         # Max growth rate
    #         base_growth = max_growth_rate[max_growth_rate["Treatment"] == t][
    #             "GrowthRate"
    #         ]
    #         mod_growth = max_growth_rate[
    #             max_growth_rate["Treatment"] == f"{self.comp.comp.modifier}_{t}"
    #         ]["GrowthRate"]

    #         u_stat, p_value = stats.mannwhitneyu(base_growth, mod_growth)
    #         print(
    #             f"Max growth rate - p-value: {p_value:.4f}, power difference: {np.mean(mod_growth) / np.mean(base_growth):.4f}, means: {np.mean(base_growth):.4f} vs {np.mean(mod_growth):.4f}"
    #         )

    #         # Time to max growth rate
    #         base_time_growth = max_growth_rate[max_growth_rate["Treatment"] == t][
    #             "Time (hrs)"
    #         ]
    #         mod_time_growth = max_growth_rate[
    #             max_growth_rate["Treatment"] == f"{self.comp.comp.modifier}_{t}"
    #         ]["Time (hrs)"]

    #         u_stat, p_value = stats.mannwhitneyu(base_time_growth, mod_time_growth)
    #         print(
    #             f"Time to max growth rate - p-value: {p_value:.4f}, power difference: {np.mean(mod_time_growth) / np.mean(base_time_growth):.4f}, means: {np.mean(base_time_growth):.4f} vs {np.mean(mod_time_growth):.4f}"
    #         )

    #         # Max speck time
    #         base_time = max_speck_time[max_speck_time["Treatment"] == t]["Time (hrs)"]
    #         mod_time = max_speck_time[
    #             max_speck_time["Treatment"] == f"{self.comp.comp.modifier}_{t}"
    #         ]["Time (hrs)"]

    #         u_stat, p_value = stats.mannwhitneyu(base_time, mod_time)
    #         print(
    #             f"Max speck time - p-value: {p_value:.4f}, power difference: {np.mean(mod_time) / np.mean(base_time):.4f}, means: {np.mean(base_time):.4f} vs {np.mean(mod_time):.4f}"
    #         )

    #         # Maximum amount of specks
    #         base_specks = self.data[self.data["Treatment"] == t][
    #             "NormalizedSpeckFormation"
    #         ]
    #         mod_specks = self.data[
    #             self.data["Treatment"] == f"{self.comp.comp.modifier}_{t}"
    #         ]["NormalizedSpeckFormation"]

    #         u_stat, p_value = stats.mannwhitneyu(base_specks, mod_specks)
    #         print(
    #             f"Max specks - p-value: {p_value:.4f}, power difference: {np.mean(mod_specks) / np.mean(base_specks):.4f}, means: {np.mean(base_specks):.4f} vs {np.mean(mod_specks):.4f}"
    #         )

    #         print()

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

            u_stat, p_value = stats.mannwhitneyu(base_growth, mod_growth)
            result["Max_growth_rate_p_value"] = p_value
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

            u_stat, p_value = stats.mannwhitneyu(base_time_growth, mod_time_growth)
            result["Time_max_growth_rate_p_value"] = p_value
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

            u_stat, p_value = stats.mannwhitneyu(base_time, mod_time)
            result["Max_speck_time_p_value"] = p_value
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

            u_stat, p_value = stats.mannwhitneyu(base_specks, mod_specks)
            result["Max_specks_p_value"] = p_value
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
            max_growth_rate_t = (
                self.max_growth_rate[self.max_growth_rate["Treatment"] == treatment][
                    "GrowthRate"
                ].mean()
            ) / 0.5
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

    def compare_mean_replicates_table(self, treatment, time1=0):
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
        results = pd.DataFrame(columns=["Time", "Mean", "Difference", "P-value"])

        # Filter the data for the treatment and time1
        data_time1 = self.data[
            (self.data["Treatment"] == treatment) 
            & (self.data["Time (hrs)"] == time1)
        ]["NormalizedSpeckFormation"]
        mean_time1 = data_time1.mean()

        # Iterate through the unique time points and compare the means
        for time2 in unique_timepoints:
            data_time2 = self.data[
                (self.data["Treatment"] == treatment)
                & (self.data["Time (hrs)"] == time2)
            ]["NormalizedSpeckFormation"]

            # Perform the Mann-Whitney U test to compare the means
            u_stat, p_value = stats.mannwhitneyu(data_time1, data_time2)
            mean_time2 = data_time2.mean()

            # Calculate the difference between mean at time2 and mean at time1
            difference = mean_time2 - mean_time1

            # Add the results to the DataFrame
            results = results.append(
                {
                    "Time": time2,
                    "Mean": mean_time2,
                    "Difference": difference,
                    "P-value": p_value,
                },
                ignore_index=True,
            )

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
        unique_timepoints = self.data[self.data["Time (hrs)"].astype(float) > time1][
            "Time (hrs)"
        ].unique()

        # Create an empty DataFrame to store the results
        results = pd.DataFrame(columns=["Analyte", "Time", "Mean", "P-value"])

        for analyte in self.data["Analyte"].unique():
            print(analyte)
            # Filter the data for the treatment and time1
            data_time1 = self.data[
                (self.data["Treatment"] == treatment)
                & (self.data["Time (hrs)"].astype(float) == float(time1))
                & (self.data["Analyte"] == analyte)
            ]["Normalized_Measurement"]
            mean_time1 = data_time1.mean()

            # Iterate through the unique time points and compare the means
            for time2 in unique_timepoints:
                data_time2 = self.data[
                    (self.data["Treatment"] == treatment)
                    & (self.data["Time (hrs)"].astype(float) == float(time2))
                    & (self.data["Analyte"] == analyte)
                ]["Normalized_Measurement"]

                # Perform the Bootstrap test to compare the means
                p, ci = self.bootstrap_t_test(data_time1, data_time2)
                mean_time2 = data_time2.mean()

                # Calculate the difference between mean at time2 and mean at time1
                difference = mean_time2 - mean_time1

                # Add the results to the DataFrame
                new_row = pd.DataFrame(
                    {
                        "Analyte": [analyte],
                        "Time": [time2],
                        "Mean": [mean_time2],
                        "P-value": [p],
                    }
                )

                results = pd.concat([results, new_row], ignore_index=True)

        # Set the DataFrame index to the Time column
        # results.set_index("Time", inplace=True)

        return results

    def aggregate_ctsc0t(self):
        # Get a list of all treatments, including modified treatments
        all_treatments = self.data["Treatment"].unique()

        # Initialize an empty dictionary to store the significant times for each treatment
        significant_times = {}
        filename = "cytokine_significant_times_by_treatment.xlsx"
        # Iterate over each treatment
        for treatment in all_treatments:
            # Calculate the time range with a significant difference from zero count
            table = self.compare_time_sequence_count_to_zero_time(treatment)
            significant_range = (
                table[table["P-value"] < 0.05].Time.min(),
                table[table["P-value"] < 0.05].Time.max(),
            )

            # Store the significant times in the dictionary
            significant_times[treatment] = significant_range

            # Write the significant times to the report
            print(
                f"{treatment} is significant between {significant_range[0]} and {significant_range[1]}."
            )
            # Check if file exists, otherwise create it
            if not os.path.isfile(filename):
                table[table["P-value"] < 0.05].to_excel(
                    filename, sheet_name=treatment, index=True
                )
            # Write the significant times to an Excel file
            with pd.ExcelWriter(
                filename, mode="a", if_sheet_exists="replace"
            ) as writer:
                table[table["P-value"] < 0.05].to_excel(
                    writer, sheet_name=treatment, index=True
                )
