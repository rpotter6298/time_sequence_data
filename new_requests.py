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
cyto.structure_data("Concentration")

speck = time_sequence_module(
    Count_Aggregator(
        "Speck", base_treatments, modifier=modifier, short_name_dict=short_name_dict
    )
)
speck.structure_data("SpeckFormation")
speck.data = speck.limit_data(speck.data, time_limits=time_limits)
speck.data = speck.data[
    ~(
        (speck.data["Treatment"].str.contains("MCC950"))
        & (speck.data["Experimental_Replicate"] == "1")
    )
]

TAS = analysis_module([cyto, speck])
TAS.compute_ratio("TS_Cyto")


def ax_modification(ax):
    ax.set_title(str("Inflammasome Speck Count over Time"))


## SPECK CURVES IN THE SAME PLOT WITH RAW COUNTS
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments="all",
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
    ratio_name="IL18:IL1b", invert=True, manual_ax_modification=ax_modification
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
        ratio_name="IL18:IL1b",
        treatments=[treatment],
        invert=True,
        manual_ax_modification=ax_modification,
    )


## Line Plots of Speck data for only MCC950 treatments
def ax_modification(ax):
    ax.set_title(str("Inflammasome Speck Count over Time - MCC950 Treatments"))
    handles, labels = ax.get_legend_handles_labels()
    new_labels = [label.replace("_", " + ") for label in labels]
    ax.legend(handles, new_labels)


plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments=["MCC950_ATP", "MCC950_MSU", "MCC950_Nigericin"],
        measurement_type="Measurement",
    ),
    manual_ax_modification=ax_modification,
)
## Line Plots of Speck data for only each MCC950 treatment
for treatment in ["MCC950_ATP", "MCC950_MSU", "MCC950_Nigericin"]:
    treatment_name = treatment.replace("_", " + ")

    def ax_modification(ax):
        ax.set_title(str("Inflammasome Speck Count over Time - " + treatment_name))
        ymin, ymax = ax.get_ylim()
        ax.set_ylim(ymin, ymax * 1.1)
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

def ax_modification(ax):
    ax.set_title("Ratio of IL-1β:IL-18 Concentration over Time")
    ax.set_ylabel("IL-1β:IL-18 Concentration Ratio")
plotting_module.plot_ratio(cyto, invert=True, manual_ax_modification=ax_modification)