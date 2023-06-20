from modules import *
from pathlib import Path
import json
from pprint import pprint

base_treatments = ["ATP", "MSU", "Nigericin"]
modifier = "MCC950"
short_name_dict = {"nig": "Nigericin"}
time_limits = {"ALL": 21}

cyto = time_sequence_module(
    Cytokine_Adapter("Cyto", base_treatments, short_name_dict=short_name_dict)
)
# cyto.data = cyto.limit_data(cyto.data, time_limits=time_limits)

speck = time_sequence_module(
    Count_Aggregator(
        "Speck", base_treatments, modifier=modifier, short_name_dict=short_name_dict
    )
)
speck.structure_data("SpeckFormation")
speck.data = speck.limit_data(speck.data, time_limits=time_limits)
# Removal of outlier data
speck.data = speck.data[
    ~(
        (speck.data["Treatment"].str.contains("MCC950"))
        & (speck.data["Experimental_Replicate"] == "1")
    )
]
speck.data = speck.data[
    ~(
        (speck.data["Treatment"] == "MSU")
        & (speck.data["Experimental_Replicate"] == "1")
    )
]
speck.data = speck.data[
    ~(
        (speck.data["Treatment"] == "MSU")
        & (speck.data["Experimental_Replicate"] == "2")
    )
]
# count the number of replicates within each treatment group in speck.data and cyto.data
for module in [speck, cyto]:
    for treatment in module.data.Treatment.unique():
        subset = module.data[module.data["Treatment"] == treatment]
        count_replicates = subset.Experimental_Replicate.unique()
        print(treatment, len(count_replicates))

TAS = analysis_module([cyto, speck])
TAS.time_compare()
TAS.compute_ratio("TS_Cyto")

##REPORTS SECTION
## Create Summaries
for module in TAS.modules:
    (TAS.create_summary_table(TAS.modules[module])).to_excel(
        Path("reports", "".join(["summary_", module, ".xlsx"])), index=False
    )
## Create Modifier Impact Reports
for module in [module for module in TAS.modules if TAS.modules[module].comp.modifier]:
    TAS.check_modifier_impact(TAS.modules[module]).T.to_excel(
        Path("reports", "".join(["modifier_impact_", module, ".xlsx"])), header=False
    )
## Create Time Point Comparisons
for module in TAS.modules:
    TAS.aggregate_time_comparisons(TAS.modules[module])


###PLOTS SECTION
##Line Plots of Speck Data With and Without Modifiers
module = TAS.modules["TS_Speck"]
for treatment in module.comp.treatments:

    def ax_modification(ax):
        ax.set_title(str(f"Normalized ASC-Speck Count with MCC950 - {treatment}"))
        ax.set_ylabel("Normalized ASC-Speck Count")
        # ymin, ymax = ax.get_ylim()
        # ax.set_ylim(ymin, ymax * 1.2)

    plotting_module.plot_lineplot(
        **TAS.prepare_lineplot(
            TAS.modules["TS_Speck"],
            treatments=[treatment, f"MCC950_{treatment}"],
            measurement_type="Normalized_Measurement",
        ),
        manual_ax_modification=ax_modification,
    )


def ax_modification(ax):
    ax.set_title(str("Inflammasome Speck Count over Time"))
    ax.set_ylabel("ASC-Speck Count")


## SPECK CURVES IN THE SAME PLOT WITH RAW COUNTS
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments="all",
        measurement_type="Measurement",
    ),
    manual_ax_modification=ax_modification,
)
## NON-MCC SPECK CURVES IN THE SAME PLOT WITH RAW COUNTS
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments=["ATP", "MSU", "Nigericin"],
        measurement_type="Measurement",
    ),
    manual_ax_modification=ax_modification,
)

## SEPARATE ATP SPECK CURVE WITH RAW COUNTS
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"], treatments=["ATP"], measurement_type="Measurement"
    ),
    manual_ax_modification=ax_modification,
)

## CYTOKINES PLOTS, ALL TREATMENTS, SINGLE ANALYTE
for analyte in cyto.data["Analyte"].unique():

    def ax_modification(ax):
        if analyte == "IL1b":
            aname = "IL-1β"
        elif analyte == "IL18":
            aname = "IL-18"
        ax.set_title(str(aname + " Concentration over Time"))
        ax.set_ylabel("Concentration (pg/mL)")
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.2)

    plotting_module.plot_lineplot(
        **TAS.prepare_lineplot(
            TAS.modules["TS_Cyto"],
            treatments="all",
            measurement_type="Measurement",
            analyte_selection=[analyte],
        ),
        manual_ax_modification=ax_modification,
    )


## SPECK COUNT, CYTO RATIO and ALL TREATMENTS
def ax_modification(
    ax,
):
    ax.set_title("Ratio of IL-1β:IL-18 Concentration over Time")
    ax.set_ylabel("IL-1β:IL-18 Concentration Ratio")


plotting_module.plot_ratio(
    TAS.modules["TS_Cyto"], invert=True, manual_ax_modification=ax_modification
)


def ax_modification(ax, ax2):
    ax.set_title("Ratio of IL-1β:IL-18 Concentration over Time")
    ax.set_ylabel("Speck Formation")
    ymin, ymax = ax.get_ylim()
    ax.set_ylim(ymin, ymax * 1.25)
    ymin, ymax = ax2.get_ylim()
    ax2.set_ylim(ymin, ymax * 1.25)
    ax2.set_ylabel("IL-1β:IL-18 Concentration Ratio")


plotting_module.plot_count_against_ratio(
    TAS, ratio_name="IL18:IL1b", invert=True, manual_ax_modification=ax_modification
)

# SAME, BUT FOR EACH TREATMENT
for treatment in cyto.data["Treatment"].unique():

    def ax_modification(ax, ax2):
        ax.set_title(f"Ratio of IL-1β:IL-18 Concentration over Time - {treatment}")
        ax.set_ylabel("Speck Formation")
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.15)
        ymin, ymax = ax2.get_ylim()
        ax2.set_ylim(ymin, ymax * 1.15)
        ax2.set_ylabel("IL-1β:IL-18 Concentration Ratio")

    plotting_module.plot_count_against_ratio(
        TAS,
        ratio_name="IL18:IL1b",
        treatments=[treatment],
        invert=True,
        manual_ax_modification=ax_modification,
    )

## Line Plots of Speck data for only each MCC950 treatment
for treatment in ["MCC950_ATP", "MCC950_MSU", "MCC950_Nigericin"]:
    treatment_name = treatment.replace("_", " + ")

    def ax_modification(ax):
        ax.set_title(str("Inflammasome Speck Count over Time - " + treatment_name))
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.1)
        # set y axis label
        ax.set_ylabel("ASC-Speck Formation")
        handles, labels = ax.get_legend_handles_labels()
        new_labels = [label.replace("_", " + ") for label in labels]
        ax.legend(handles, new_labels)

    plotting_module.plot_lineplot(
        **TAS.prepare_lineplot(
            TAS.modules["TS_Speck"],
            treatments=[treatment],
            measurement_type="Measurement",
        ),
        manual_ax_modification=ax_modification,
    )

## Line Plots of Speck data for only each non-MCC950 treatment
for treatment in ["ATP", "MSU", "Nigericin"]:
    treatment_name = treatment.replace("_", " + ")

    def ax_modification(ax):
        ax.set_title(
            str("Inflammasome Speck Count over T#EDCvfr4ime - " + treatment_name)
        )
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.1)
        # set y axis label
        ax.set_ylabel("ASC-Speck Formation")
        handles, labels = ax.get_legend_handles_labels()
        new_labels = [label.replace("_", " + ") for label in labels]
        ax.legend(handles, new_labels)

    plotting_module.plot_lineplot(
        **TAS.prepare_lineplot(
            TAS.modules["TS_Speck"],
            treatments=[treatment],
            measurement_type="Measurement",
        ),
        manual_ax_modification=ax_modification,
    )


# for treatment in module.data["Treatment"].unique():
#     dict = TAS.prepare_lineplot(module, [treatment])
#     plotting_module.plot_lineplot(**dict)


# for treatment in module.comp.treatments:
#     mod_treat = "_".join([module.comp.modifier, treatment])
#     TAS.plot_speck_lineplots(
#         module,
#         [treatment, mod_treat],
#         output_path=Path("plots", "speck_plots", "modifier_comparison"),
#         show_p=True,
#     )

# ##Line Plots of Cytokine Data
# module = TAS.modules["TS_Cyto"]
# for treatment in module.data["Treatment"].unique():
#     TAS.plot_speck_lineplots(
#         module, [treatment], output_path=Path("plots", "cytokine_plots"), show_p=True
#     )

# ##Multimodel Plots
# TAS.plot_multimodel(
#     TAS.modules["TS_Speck"],
#     TAS.modules["TS_Cyto"],
#     output_directory=Path("plots", "multimodel_plots"),
# )

# # Arbitrary Maximums
# for module in TAS.modules:
#     TAS.find_arbitrary_max(TAS.modules[module])
#     nested_dict = {}
#     for key, value in TAS.modules[module].arbitrary_maximums.__dict__.items():
#         value_dict = value.__dict__.copy()
#         if "window" in value_dict:
#             del value_dict["window"]
#         nested_dict[key] = value_dict
#     with open(Path("reports", "arbitrary_max_" + module + ".json"), "w") as outfile:
#         json.dump(nested_dict, outfile)

# ## Is there a difference between the modifier and the base treatment at the arbitrary maximum time?
# module = TAS.modules["TS_Speck"]
# arb_max_diff = {}
# for treatment in module.comp.treatments:
#     mod_treat = "_".join([module.comp.modifier, treatment])
#     arbitrary_max_time_treat = module.arbitrary_maximums.__dict__[treatment].time
#     arbitrary_max_time_mod_treat = module.arbitrary_maximums.__dict__[mod_treat].time

#     # The measurements for the treatment, no mod, at the arbitrary maximum for the treatment without the modifier (A1)
#     A1 = module.data[
#         (module.data["Treatment"] == treatment)
#         & (module.data["Time (hrs)"] == arbitrary_max_time_treat)
#     ]["Measurement"]
#     # The measurements for the treatment, no mod, at the time of the arbitrary maximum for the treatment with the modifier (A2)
#     A2 = module.data[
#         (module.data["Treatment"] == treatment)
#         & (module.data["Time (hrs)"] == arbitrary_max_time_mod_treat)
#     ]["Measurement"]
#     # The measurements for the treatment, with mod, at the time of the arbitrary maximum for the treatment without the modifier (B1)
#     B1 = module.data[
#         (module.data["Treatment"] == mod_treat)
#         & (module.data["Time (hrs)"] == arbitrary_max_time_treat)
#     ]["Measurement"]
#     # The measurements for the treatment, with mod, at the arbitrary maximum for the treatment with the modifier (B2)
#     B2 = module.data[
#         (module.data["Treatment"] == mod_treat)
#         & (module.data["Time (hrs)"] == arbitrary_max_time_mod_treat)
#     ]["Measurement"]
#     p, ci = bootstrap_t_test(A1, B1)
#     arb_max_diff[treatment] = {
#         "Time1": arbitrary_max_time_treat,
#         "Time2": arbitrary_max_time_mod_treat,
#         "A1-B1": {
#             "p": p,
#             "ci": ci,
#             "mean(treatment)": round(A1.mean(), 3),
#             "mean(modified_treatment)": round(B1.mean(), 3),
#         },
#     }
#     if arbitrary_max_time_treat == arbitrary_max_time_mod_treat:
#         continue
#     p, ci = bootstrap_t_test(A2, B2)
#     arb_max_diff[treatment]["A2-B2"] = {
#         "p": p,
#         "ci": ci,
#         "mean(treatment)": round(A2.mean(), 3),
#         "mean(modified_treatment)": round(B2.mean(), 3),
#     }
#     p, ci = bootstrap_t_test(A1, B2)
#     arb_max_diff[treatment]["A1-B2"] = {
#         "p": p,
#         "ci": ci,
#         "mean(treatment)": round(A1.mean(), 3),
#         "mean(modified_treatment)": round(B2.mean(), 3),
#     }


# # Test if there is a difference between the times of the max speck and max cyto for each treatment
# report = TAS.compare_max_time_distances(TAS.modules["TS_Speck"], TAS.modules["TS_Cyto"])
# pprint.pprint(report)
# # Test if the cytokine FC max is different between treatments
# report = TAS.compare_max_normalized_measurement_anova(TAS.modules["TS_Cyto"].data)
# pprint.pprint(report)
# for analyte in report.keys():
#     print(analyte)
#     df = pd.DataFrame(report[analyte]["post_hoc"])
#     print(df)
