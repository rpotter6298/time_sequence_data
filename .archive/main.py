import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats
from pathlib import Path
from typing import Type, List, Dict
import re




class time_sequence:
    def __init__(self, file_path, limits: dict = None, mode=0):
        self.raw_data = pd.read_excel(file_path, sheet_name="Sheet1")
        self.limits = limits
        if mode == 0:
            self.data = self.normalize()
        else:
            self.data = self.raw_data
        self.treatments = self.data["Treatment"].unique()

    def individual_counts_normalize(self):
        file = pd.read_excel('count_all_data.xlsx', sheet_name='ATP rep 1.1')

    def normalize(self):
        data = self.raw_data.set_index("Time").stack().reset_index()
        data.columns = ["Time", "Treatment_Replicate", "SpeckFormation"]
        data["Treatment"] = data["Treatment_Replicate"].str.rsplit("_", n=1).str[0]
        data["Replicate"] = data["Treatment_Replicate"].str.rsplit("_", n=1).str[1]
        # Get the values of SpeckFormation at Time 0 for each Treatment_Replicate
        initial_values = data[data["Time"] == 0].set_index("Treatment_Replicate")[
            "SpeckFormation"
        ]

        # Divide the SpeckFormation values by their corresponding value at Time 0
        data["NormalizedSpeckFormation"] = data.apply(
            lambda row: row["SpeckFormation"]
            / initial_values[row["Treatment_Replicate"]],
            axis=1,
        )
        if self.limits is not None:
            data = self.limit_data(data, self.limits)
        return data

    def line_plot(self):
        fig, ax = plt.subplots()
        sns.lineplot(
            ax=ax,
            data=self.data,
            x="Time",
            y="NormalizedSpeckFormation",
            hue="Treatment",
        )
        ax.set_title("Speck Formation by Treatment Type")
        ax.set_ylabel("Normalized Speck Formation")
        ax.set_xlabel("Time")
        return (fig, ax)

    def limit_data(self, data: pd.DataFrame, limits: dict = None):
        limited_data = data
        if limits is None:
            print("No limits given, returning all data")
            return data
        else:
            for treatment in limits.keys():
                # print(treatment)
                if treatment == "All":
                    limited_data = limited_data[data["Time"] <= limits[treatment]]
                else:
                    limited_data = limited_data[
                        (data["Treatment"].str.contains(treatment, case=False))
                        & (data["Time"] <= limits[treatment])
                        | (~data["Treatment"].str.contains(treatment, case=False))
                    ]
            return limited_data

    def define_treatment_modifier(self, modifier: str):
        self.treatment_modifier = modifier
        self.treatments = (
            self.data["Treatment"].str.replace(modifier + "_", "").unique()
        )

    def calc_growth_rate(self):
        growth_rate = self.data.copy()
        growth_rate["GrowthRate"] = growth_rate.groupby(
            ["Treatment", "Treatment_Replicate"]
        )["NormalizedSpeckFormation"].diff()
        self.growth_rate = growth_rate
        self.data["GrowthRate"] = growth_rate["GrowthRate"]

    def compare_max_speck_time(self):
        # Calculate the time to reach the maximum speck formation for each treatment replicate
        max_speck_time = self.data.loc[
            self.data.groupby(["Treatment", "Treatment_Replicate"])[
                "NormalizedSpeckFormation"
            ].idxmax()
        ]
        self.max_speck_time = max_speck_time
        # Perform Mann-Whitney U tests to compare the time to reach the maximum speck formation between treatments
        print("Time to reach the maximum speck formation:")
        for i, t1 in enumerate(self.treatments):
            for t2 in self.treatments[i + 1 :]:
                self.compare_treatments_original(t1, t2, max_speck_time, "Time")
            if self.treatment_modifier:
                self.compare_treatments_original(
                    t1, str(self.treatment_modifier + "_" + t1), max_speck_time, "Time"
                )

    def compare_treatments_original(self, treatment1, treatment2, df, metric):
        data1 = df[df["Treatment"] == treatment1][metric]
        data2 = df[df["Treatment"] == treatment2][metric]

        u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative="two-sided")
        print(f"{treatment1} vs {treatment2} ({metric}): p-value = {p_value}")

    def compare_max_growth_rate(self):
        # Calculate the maximum growth rate for each treatment replicate
        max_growth_rate = self.growth_rate.loc[
            self.growth_rate.groupby(["Treatment", "Treatment_Replicate"])[
                "GrowthRate"
            ].idxmax()
        ]
        self.max_growth_rate = max_growth_rate
        # Perform Mann-Whitney U tests to compare the maximum growth rate between treatments
        print("Maximum growth rate comparison:")
        for i, t1 in enumerate(self.treatments):
            for t2 in self.treatments[i + 1 :]:
                self.compare_treatments_original(t1, t2, max_growth_rate, "GrowthRate")
            if self.treatment_modifier:
                self.compare_treatments_original(
                    t1,
                    str(self.treatment_modifier + "_" + t1),
                    max_growth_rate,
                    "GrowthRate",
                )

    def mannwhitney_pvalue(self, treatment1, treatment2, metric):
        data1 = self.data[self.data["Treatment"] == treatment1][metric]
        data2 = self.data[self.data["Treatment"] == treatment2][metric]

        # Perform the Mann-Whitney U test
        result = stats.mannwhitneyu(data1, data2, alternative="two-sided")

        # Return the p-value
        return result.pvalue

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
        all_treatments = list(self.treatments) + [
            f"{self.treatment_modifier}_{treatment}" for treatment in self.treatments
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
            ]["Time"].mean()
            time_max_growth_rate_t = self.max_growth_rate[
                self.max_growth_rate["Treatment"] == treatment
            ]["Time"].mean()

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



speck_class = time_sequence(Path("speck_data.xlsx"), limits={"All": 21})
speck_class.define_treatment_modifier("MCC950")
speck_class.calc_growth_rate()
speck_class.compare_max_speck_time()
speck_class.compare_max_growth_rate()
speck_class.create_summary_table()
speck_class.compare_max_speck_time()
speck_class.compare_max_growth_rate()
speck_class.calc_growth_rate()
speck_class.growth_rate
max_speck_time = speck_class.data.loc[
    speck_class.data.groupby(["Treatment", "Treatment_Replicate"])[
        "NormalizedSpeckFormation"
    ].idxmax()
]
df = max_speck_time
treatment1 = "ATP"
treatment2 = "Nigericin"
metric = "Time"

speck_class.normalize()
speck_class.data["Treatment"].unique()
fig, ax = speck_class.line_plot()
fig.show()
speck_class.limit_data(speck_class.data, limits={"All": 21, "nigericin": 3})

xls = pd.read_excel('count_all_data.xlsx', sheet_name=None)
treatment_options = ["ATP", "MSU", "Nigericin", "MCC950_ATP", "MCC950_MSU", "MCC950_Nigericin"]
def match_treatment(sheet_name, treatment_options):
    for treatment in treatment_options:
        if treatment.lower().replace("_", "").replace(" ", "") in sheet_name.lower().replace("_", "").replace(" ", ""):
            return treatment
all_sheets = []
for sheet_name, df in xls.items():
    # Remove the top row and set the second row as column names
    df = df.drop(0).reset_index(drop=True)
    df.columns = df.iloc[0]
    
    # Add a column with the sheet name
    df["SheetName"] = sheet_name
    
    # Extract experimental and technical replicates
    decimal_match = re.search(r"(\d+)\.?(\d+)?", sheet_name)
    if decimal_match:
        df["ExperimentalReplicate"] = decimal_match.group(1)
        df["TechnicalReplicate"] = decimal_match.group(2) if decimal_match.group(2) else 1

    # Remove "rep" from the sheet name if it exists
    modified_sheet_name = sheet_name.lower().replace("rep", "")

    # Create a column for treatment
    df["Treatment"] = match_treatment(modified_sheet_name, treatment_options)
    
    all_sheets.append(df)

# Combine all sheets into a single DataFrame
combined_data = pd.concat(all_sheets, ignore_index=True)

