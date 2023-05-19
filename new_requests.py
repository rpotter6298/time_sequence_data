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

## SPECK CURVES IN THE SAME PLOT WITH RAW COUNTS
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"], treatments="all", measurement_type="Measurement"
    )
)
## SEPARATE ATP SPECK CURVE WITH RAW COUNTS
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"], treatments=["ATP"], measurement_type="Measurement"
    )
)
## CYTOKINES PLOTS, ALL TREATMENTS, SINGLE ANALYTE
for analyte in cyto.data["Analyte"].unique():

    def ax_modification(ax):
        ax.set_title(str(analyte + " Concentration over Time"))
        ax.set_ylabel("Concentration (pg/mL)")

    plotting_module.plot_lineplot(
        **TAS.prepare_lineplot(
            TAS.modules["TS_Cyto"],
            treatments="all",
            measurement_type="Measurement",
            analyte_selection=[analyte],
        ),
        manual_ax_modification=ax_modification
    )

## SPECK COUNT, CYTO RATIO, EACH TREATMENT, and ALL TREATMENTS


def ax_modification(ax):
    ax.set_title("Ratio of IL1b:IL18 Concentration over Time")
    ax.set_ylabel("IL1b:IL18 Concentration Ratio")


plotting_module.plot_ratio(
    TAS.modules["TS_Cyto"], invert=True, manual_ax_modification=ax_modification
)
plotting_module.plot_count_against_ratio(ratio_name="IL18:IL1b", invert=True)
plotting_module.plot_count_against_ratio(
    ratio_name="IL18:IL1b", treatments=["ATP"], invert=True
)
plotting_module.plot_count_against_ratio(
    ratio_name="IL18:IL1b", treatments=["MSU"], invert=True
)
plotting_module.plot_count_against_ratio(
    ratio_name="IL18:IL1b", treatments=["Nigericin"], invert=True
)

## Line Plots of Speck data for only MCC950 treatments
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments=["MCC950_ATP", "MCC950_MSU", "MCC950_Nigericin"],
        measurement_type="Measurement",
    )
)
## Line Plots of Speck data for only each MCC950 treatment
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments=["MCC950_ATP"],
        measurement_type="Measurement",
    )
)
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments=["MCC950_MSU"],
        measurement_type="Measurement",
    )
)
plotting_module.plot_lineplot(
    **TAS.prepare_lineplot(
        TAS.modules["TS_Speck"],
        treatments=["MCC950_Nigericin"],
        measurement_type="Measurement",
    )
)
