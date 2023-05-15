# import os
# os.chdir('time_sequence_data')
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats as stats


# read in the Excel file
data = pd.read_excel("speck_data.xlsx", sheet_name="Sheet1")
data = data.set_index("Time").stack().reset_index()
data.columns = ["Time", "Treatment_Replicate", "SpeckFormation"]
data["Treatment"] = data["Treatment_Replicate"].str.rsplit("_", n=1).str[0]

# Get the values of SpeckFormation at Time 0 for each Treatment_Replicate
initial_values = data[data["Time"] == 0].set_index("Treatment_Replicate")[
    "SpeckFormation"
]

# Divide the SpeckFormation values by their corresponding value at Time 0
data["NormalizedSpeckFormation"] = data.apply(
    lambda row: row["SpeckFormation"] / initial_values[row["Treatment_Replicate"]],
    axis=1,
)

# Create Line Plot
sns.lineplot(data=limited_data, x="Time", y="NormalizedSpeckFormation", hue="Treatment")
plt.title("Speck Formation by Treatment Type")
plt.show()

mean_speck_by_treatment = data.groupby(["Treatment", "Time"]).mean().reset_index()
mean_speck_by_treatment["Difference"] = mean_speck_by_treatment.groupby("Treatment")[
    "NormalizedSpeckFormation"
].diff()

sns.lineplot(data=mean_speck_by_treatment, x="Time", y="Difference", hue="Treatment")
plt.title("Rate of Speck Formation by Treatment Type")
plt.show()

treatments = ["ATP", "MSU", "Nigericin"]
mean_normalized_speck_by_treatment = (
    data.groupby(["Treatment", "Time"]).mean("NormalizedSpeckFormation").reset_index()
)

for treatment in treatments:
    treatment_data = mean_normalized_speck_by_treatment[
        mean_normalized_speck_by_treatment["Treatment"].isin(
            [treatment, f"MCC950_{treatment}"]
        )
    ]
    plt.figure()
    sns.lineplot(
        data=treatment_data, x="Time", y="NormalizedSpeckFormation", hue="Treatment"
    )
    plt.title(f"Normalized Speck Formation: {treatment} vs MCC950_{treatment}")
    plt.show()

nigericin_limit = 3
other_limit = 21

# Limit the data based on the treatment type
limited_data = mean_normalized_speck_by_treatment[
    (
        (mean_normalized_speck_by_treatment["Treatment"].str.contains("Nigericin"))
        & (mean_normalized_speck_by_treatment["Time"] <= nigericin_limit)
    )
    | (
        (~mean_normalized_speck_by_treatment["Treatment"].str.contains("Nigericin"))
        & (mean_normalized_speck_by_treatment["Time"] <= other_limit)
    )
]

# Find the time to reach the maximum speck formation for each treatment
max_speck_time = limited_data.loc[
    limited_data.groupby("Treatment")["NormalizedSpeckFormation"].idxmax()
]

print(max_speck_time[["Treatment", "Time"]])


def compare_treatments(treatment1, treatment2, df):
    data1 = df[df["Treatment"] == treatment1]["NormalizedSpeckFormation"]
    data2 = df[df["Treatment"] == treatment2]["NormalizedSpeckFormation"]

    u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative="two-sided")
    print(f"{treatment1} vs {treatment2}: p-value = {p_value}")


for i, t1 in enumerate(treatments):
    for t2 in treatments[i + 1 :]:
        compare_treatments(t1, t2, mean_normalized_speck_by_treatment)

for treatment in treatments:
    compare_treatments(
        treatment, f"MCC950_{treatment}", mean_normalized_speck_by_treatment
    )


nigericin_limit = 3
other_limit = 21

# Limit the data based on the treatment type
limited_data = mean_normalized_speck_by_treatment[
    (
        (mean_normalized_speck_by_treatment["Treatment"].str.contains("Nigericin"))
        & (mean_normalized_speck_by_treatment["Time"] <= nigericin_limit)
    )
    | (
        (~mean_normalized_speck_by_treatment["Treatment"].str.contains("Nigericin"))
        & (mean_normalized_speck_by_treatment["Time"] <= other_limit)
    )
]

# Calculate the growth rate for the limited data
growth_rate = limited_data.copy()
growth_rate["GrowthRate"] = growth_rate.groupby("Treatment")[
    "NormalizedSpeckFormation"
].diff()

# Find the time to reach the maximum growth rate for each treatment
max_growth_rate_time = growth_rate.loc[
    growth_rate.groupby("Treatment")["GrowthRate"].idxmax()
]

print("Time to reach the maximum growth rate:")
print(max_growth_rate_time[["Treatment", "Time"]])


original_growth_rate = data.copy()
original_growth_rate["GrowthRate"] = original_growth_rate.groupby(
    ["Treatment", "Treatment_Replicate"]
)["NormalizedSpeckFormation"].diff()

max_speck_time_original = original_growth_rate.loc[
    original_growth_rate.groupby(["Treatment", "Treatment_Replicate"])[
        "NormalizedSpeckFormation"
    ].idxmax()
]
max_growth_rate_time_original = original_growth_rate.loc[
    original_growth_rate.groupby(["Treatment", "Treatment_Replicate"])[
        "GrowthRate"
    ].idxmax()
]


def compare_treatments_original(treatment1, treatment2, df, metric):
    data1 = df[df["Treatment"] == treatment1][metric]
    data2 = df[df["Treatment"] == treatment2][metric]

    u_stat, p_value = stats.mannwhitneyu(data1, data2, alternative="two-sided")
    print(f"{treatment1} vs {treatment2} ({metric}): p-value = {p_value}")


print("Time to reach the maximum speck formation:")
for i, t1 in enumerate(treatments):
    for t2 in treatments[i + 1 :]:
        compare_treatments_original(t1, t2, max_speck_time_original, "Time")

print("\nTime to reach the maximum growth rate:")
for i, t1 in enumerate(treatments):
    for t2 in treatments[i + 1 :]:
        compare_treatments_original(t1, t2, max_growth_rate_time_original, "Time")
