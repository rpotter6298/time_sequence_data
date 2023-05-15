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



def polynomial(x, *coeffs):
    return np.polyval(coeffs, x)


class time_sequence_module:
    def __init__(self, data_class, time_limits: dict = None):
        self.comp = data_class
        self.raw_data = data_class.data
        self.time_limits = time_limits
        self.data = self.normalize()

    def normalize(self, robust=True):
        tech_merge = input("Merge Technical Replicates? (y/n): ")
        if tech_merge == "y":
            data = self.merge_technical_replicates()
        elif tech_merge == "n":
            data = self.preserve_technical_replicates()
        else:
            print("Invalid input.")
            return
        initial_values = data[data["Time (hrs)"] == 0].set_index("Treatment_Replicate")[
            "SpeckFormation"
        ]
        data["NormalizedSpeckFormation"] = data.apply(
            lambda row: row["SpeckFormation"]
            / initial_values[row["Treatment_Replicate"]],
            axis=1,
        )
        if self.time_limits is not None:
            data = self.limit_data(data, self.time_limits)
        if robust == True:
            data = self.calc_growth_rate(data)
        return data

    def limit_data(self, data, time_limits):
        limited_data = data
        if time_limits is None:
            print("No time limits specified. Returning all data.")
            return limited_data
        else:
            for treatment in time_limits.keys():
                if treatment == "ALL":
                    limited_data = limited_data[
                        limited_data["Time (hrs)"] <= time_limits[treatment]
                    ]
                else:
                    limited_data = limited_data[
                        (limited_data["Treatment"].str.contains(treatment, case=False))
                        & (limited_data["Time (hrs)"] <= time_limits[treatment])
                        | (
                            ~limited_data["Treatment"].str.contains(
                                treatment, case=False
                            )
                        )
                    ]
            return limited_data

    def merge_technical_replicates(self):
        data = self.raw_data.copy()
        data["Treatment_Replicate"] = (
            data["Treatment"] + "_" + data["ExperimentalReplicate"]
        )
        merged_data = (
            data.groupby(["Treatment_Replicate", "Time (hrs)"])["SpeckFormation"]
            .mean()
            .reset_index()
        )
        # Merge the grouped data back with the original data, dropping the "SpeckFormation" column
        data = data.drop("SpeckFormation", axis=1).merge(
            merged_data, on=["Treatment_Replicate", "Time (hrs)"]
        )
        data["TechnicalReplicate"] = pd.to_numeric(
            data["TechnicalReplicate"], errors="coerce"
        )
        data = data[data["TechnicalReplicate"] <= 1]
        # Set the remaining "TechnicalReplicate" values to 0
        data["TechnicalReplicate"] = 0
        # Assign the merged data back to self.raw_data
        return data

    def preserve_technical_replicates(self):
        data = self.raw_data.copy()
        data["Treatment_Replicate"] = (
            data["Treatment"]
            + "_"
            + data["ExperimentalReplicate"].astype(str)
            + "_"
            + data["TechnicalReplicate"].astype(str)
        )
        return data

    def calc_growth_rate(self, data):
        growth_rate = data.copy()
        growth_rate["GrowthRate"] = growth_rate.groupby(
            ["Treatment", "Treatment_Replicate"]
        )["NormalizedSpeckFormation"].diff()
        data["GrowthRate"] = growth_rate["GrowthRate"]
        return data

class analysis_module:
    def __init__(self, module):
        self.comp = module
        self.data = module.data
        self.max_growth_rate = self.compare_growth_rate()
        self.max_speck_time = self.compare_max_speck_time()

    def line_plot(self):
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

    def compare_max_speck_time(self):
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
        data1 = df[df["Treatment"] == treatment1][metric]
        data2 = df[df["Treatment"] == treatment2][metric]

        u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative="two-sided")
        print(f"{treatment1} vs {treatment2} ({metric}): p-value = {p_value}")

    def check_modifier_impact(self):
        if not self.comp.comp.modifier:
            print("No modifier specified.")
            return

        # Calculate maximum speck time and maximum growth rate
        max_speck_time = self.compare_max_speck_time()
        max_growth_rate = self.compare_growth_rate()

        # Loop through the treatments and compare them with their modified counterparts
        for t in self.comp.comp.treatments:
            print(f"Comparing treatment {t} with {self.comp.comp.modifier}_{t}:")

            # Max growth rate
            base_growth = max_growth_rate[max_growth_rate["Treatment"] == t][
                "GrowthRate"
            ]
            mod_growth = max_growth_rate[
                max_growth_rate["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["GrowthRate"]

            u_stat, p_value = stats.mannwhitneyu(base_growth, mod_growth)
            print(
                f"Max growth rate - p-value: {p_value:.4f}, power difference: {np.mean(mod_growth) / np.mean(base_growth):.4f}, means: {np.mean(base_growth):.4f} vs {np.mean(mod_growth):.4f}"
            )

            # Time to max growth rate
            base_time_growth = max_growth_rate[max_growth_rate["Treatment"] == t][
                "Time (hrs)"
            ]
            mod_time_growth = max_growth_rate[
                max_growth_rate["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["Time (hrs)"]

            u_stat, p_value = stats.mannwhitneyu(base_time_growth, mod_time_growth)
            print(
                f"Time to max growth rate - p-value: {p_value:.4f}, power difference: {np.mean(mod_time_growth) / np.mean(base_time_growth):.4f}, means: {np.mean(base_time_growth):.4f} vs {np.mean(mod_time_growth):.4f}"
            )

            # Max speck time
            base_time = max_speck_time[max_speck_time["Treatment"] == t]["Time (hrs)"]
            mod_time = max_speck_time[
                max_speck_time["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["Time (hrs)"]

            u_stat, p_value = stats.mannwhitneyu(base_time, mod_time)
            print(
                f"Max speck time - p-value: {p_value:.4f}, power difference: {np.mean(mod_time) / np.mean(base_time):.4f}, means: {np.mean(base_time):.4f} vs {np.mean(mod_time):.4f}"
            )

            # Maximum amount of specks
            base_specks = self.data[self.data["Treatment"] == t][
                "NormalizedSpeckFormation"
            ]
            mod_specks = self.data[
                self.data["Treatment"] == f"{self.comp.comp.modifier}_{t}"
            ]["NormalizedSpeckFormation"]

            u_stat, p_value = stats.mannwhitneyu(base_specks, mod_specks)
            print(
                f"Max specks - p-value: {p_value:.4f}, power difference: {np.mean(mod_specks) / np.mean(base_specks):.4f}, means: {np.mean(base_specks):.4f} vs {np.mean(mod_specks):.4f}"
            )

            print()

    def create_summary_table(self):
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
        # Get all unique time points after time1
        unique_timepoints = self.data[self.data["Time (hrs)"] > time1][
            "Time (hrs)"
        ].unique()

        # Create an empty DataFrame to store the results
        results = pd.DataFrame(columns=["Time", "Mean", "Difference", "P-value"])

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

    def plot_time_course_analysis(self, treatment1=None, treatment2=None, degree=3, points_only=False):
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

        # Calculate the ABC and RMSE for each pair of replicates
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

            comparison_results.append((rep1, rep2, abc, rmse))

        # Build the table with the comparison results
        results_table = pd.DataFrame(
            comparison_results, columns=["Replicate 1", "Replicate 2", "ABC", "RMSE"]
        )

        # Calculate the average ABC and RMSE
        avg_abc = results_table["ABC"].mean()
        avg_rmse = results_table["RMSE"].mean()

        print("Comparison results:")
        print(results_table)
        print(f"Average ABC: {avg_abc:.4f}")
        print(f"Average RMSE: {avg_rmse:.4f}")


testbox = time_sequence_module(structure, time_limits={"ALL": 21, "Nigericin": 12})
testbox.data

pd.set_option("display.max_rows", 20)

a1 = analysis_module(testbox)
a1.compare_replicate_curves("Nigericin")
a1.quantify_curve_difference("Nigericin", "MCC950_Nigericin")
a1.plot_time_course_analysis("MSU")
a1.plot_time_course_analysis("ATP", points_only=True)
a1.compare_mean_replicates_table("ATP", 0)
a1.create_summary_table()
a1.data
a1.comp.comp.treatments
a1.check_modifier_impact()
a1.compare_max_speck_time()
a1.compare_growth_rate()
a1.max_growth_rate
a1.max_growth_rate[a1.max_growth_rate["Treatment"] == "Nigericin"]
atp_numbers = a1.max_growth_rate[a1.max_growth_rate["Treatment"] == "ATP"]["max_growth"]
mcc950_atp_numbers = a1.max_growth_rate[
    a1.max_growth_rate["Treatment"] == "MCC950_ATP"
]["max_growth"]
atp_numbers.mean()
mcc950_atp_numbers.mean()


a1.max_growth_rate[a1.max_growth_rate["Treatment"] == "MCC950_Nigericin"]
a1.max_speck_time[a1.max_speck_time["Treatment"] == "Nigericin"]
a1.data[a1.data["Treatment"] == "Nigericin"]
fig, ax = a1.line_plot()
fig.show()
testbox.normalize()
