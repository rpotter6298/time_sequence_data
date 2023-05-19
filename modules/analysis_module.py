from typing import Type, List, Dict, Any
from modules.stats_functions import *
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
                max_factor = max_factor.drop_duplicates(
                    subset=["Analyte", "Treatment", "Experimental_Replicate"],
                    keep="first",
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
                time_0_set = module.data[
                    (module.data["Treatment"] == treatment)
                    & (module.data["Analyte"] == analyte)
                    & (module.data["Time (hrs)"] == 0)
                ]
                max_measurement = set["Measurement"].mean()

                max_measuremnt_p, ci = bootstrap_t_test(
                    time_0_set["Measurement"], set["Measurement"]
                )
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
                        "Max Measurement P-Value": max_measuremnt_p,
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
            split_data[analyte] = module.data[module.data["Analyte"] == analyte]
        return split_data

    def prepare_lineplot(
        self,
        module,
        treatments: List = None,
        output_path: Path = None,
        measurement_type=None,
        analyte_selection=None,
    ):
        if treatments is None:
            treatments = module.comp.treatments
        if treatments == "all":
            treatments = module.data["Treatment"].unique()
        data = module.data[
            module.data["Treatment"].isin(treatments)
        ]  # Filter the data to only include the specified treatments
        if analyte_selection is not None:
            data = data[data["Analyte"].isin(analyte_selection)]
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
        if measurement_type:
            y = measurement_type
        elif len(data[hue].unique()) == 1:
            y = "Measurement"
        else:
            y = "Normalized_Measurement"

        # Specify the order of treatments to control the colors
        hue_order = None  # Initialize hue_order as None

        if hue == "Treatment":
            hue_order = treatments  # Set hue_order only when hue is "Treatment"

        if output_path is not None:
            filename = output_path / f"{module.name}_{'_'.join(treatments)}.png"
            lineplot_dict = {
                "dataframe": data,
                "y": y,
                "hue": hue,
                "hue_order": hue_order,
                "filepath": filename,
            }
        else:
            lineplot_dict = {
                "dataframe": data,
                "y": y,
                "hue": hue,
                "hue_order": hue_order,
            }
        return lineplot_dict

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

    def compute_ratio(
        self,
        module_name,
        measurement_type="Measurement",
        normalize_start=False,
    ):
        module = self.modules[module_name]
        ##If the module's data does not have two analytes, return an error
        if len(module.data["Analyte"].unique()) < 2:
            raise ValueError("Model does not have two analytes.")
        self.modules[module_name].ratio_data = {}
        for analyte_pair in list(combinations(module.data["Analyte"].unique(), 2)):
            split_data = self.split_analytes(module_name, analyte_pair)
            df1, df2 = split_data.values()
            analyte1, analyte2 = split_data.keys()
            merged_data = pd.merge(
                df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
            )
            merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
                merged_data[f"{measurement_type}_x"]
                / merged_data[f"{measurement_type}_y"]
            )
            merged_data[f"Ratio_{analyte2}_{analyte1}"] = (
                merged_data[f"{measurement_type}_y"]
                / merged_data[f"{measurement_type}_x"]
            )

            # If normalize_start is True, normalize the ratio to the first time point for each treatment and experimental replicate
            if normalize_start:
                merged_data[f"Ratio_{analyte1}_{analyte2}"] = merged_data.groupby(
                    ["Treatment", "Experimental_Replicate"]
                )[f"Ratio_{analyte1}_{analyte2}"].transform(lambda x: x / x.iloc[0])
            self.modules[module_name].ratio_data[
                analyte1 + ":" + analyte2
            ] = merged_data

    def split_analytes(self, module_name, analyte_pair):
        """
        Split the data into separate dataframes for each analyte.

        Args:
            model (time_sequence): The time_sequence object to split.

        Returns:
            dict: A dictionary of dataframes for each analyte.
        """
        module = self.modules[module_name]
        split_data = {}
        for analyte in analyte_pair:
            split_data[analyte] = module.data[module.data["Analyte"] == analyte]
        return split_data
