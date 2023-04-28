from modules import *
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


def main():
    ##Example Usage:
    file_path = Path("data", "count_all_data.xlsx")
    base_treatments = ["ATP", "MSU", "Nigericin"]
    modifier = "MCC950"
    short_name_dict = {"nig": "Nigericin"}
    aggregated = Count_Aggregator(
        file_path, base_treatments, modifier, short_name_dict=short_name_dict
    )
    aggregated_mod = aggregated
    filtered_data = aggregated_mod.data[
        ~(
            (aggregated_mod.data["Treatment"].str.contains("MCC950"))
            & (aggregated_mod.data["ExperimentalReplicate"] == str(1))
        )
    ]
    aggregated_mod.data = filtered_data

    # file_path = "speck_data.xlsx"
    # summary = Summary_Adapter(file_path, base_treatments, modifier, short_name_dict)
    ts_module = time_sequence_module(aggregated_mod, time_limits={"ALL": 21})
    working_model = ts_analysis_module(ts_module)
    return working_model


base_treatments = ["ATP", "MSU", "Nigericin"]
modifier = "MCC950"
short_name_dict = {"nig": "Nigericin"}

working_model = main()
data_module = Cytokine_Adapter(base_treatments, short_name_dict=short_name_dict)
data_module.data
data_module.data["Treatment"].unique()
data_module.data.sort_values(by=["Treatment", "ExperimentalReplicate"])

untreated = working_model.comp.comp.treatments
treated = [
    working_model.comp.comp.modifier + "_" + t
    for t in working_model.comp.comp.treatments
]
combined = untreated + treated

plt, ax = working_model.line_plot_selected(combined)
plt.savefig("plots/line_plots/combined.png", dpi=300)
for name, tlist in {"untreated": untreated, "treated": treated}.items():
    plt, ax = working_model.line_plot_selected(tlist, show_p=True)
    plt.savefig(f"plots/line_plots/{name}.png", dpi=300)
for combo in [[untreated[i], treated[i]] for i in range(len(untreated))]:
    plt, ax = working_model.line_plot_selected(combo, show_p=True)
    plt.savefig(f"plots/line_plots/{combo[0]}_{combo[1]}.png", dpi=300)


for treatment in combined:
    print(treatment)
    plt, info = working_model.compare_standardized_replicate_curves(treatment)
    plt.savefig(f"plots/curve_comparisons/{treatment}.png", dpi=300)
    plt.close()


plt, ax = working_model.line_plot_selected(["ATP", "MCC950_ATP"], show_p=True)
plt.savefig("ATP.png", dpi=300)
plt, ax = working_model.line_plot_selected(["MSU", "MCC950_MSU"], show_p=True)
plt.savefig("MSU.png", dpi=300)
plt, ax = working_model.line_plot_selected(
    ["Nigericin", "MCC950_Nigericin"], show_p=True
)
plt.savefig("Nigericin.png", dpi=300)

##CHECK IF THE MODIFIER CHANGES THE BEHAVIOR AND BY HOW MUCH
impact = working_model.check_modifier_impact()
impact.T.to_excel("MCC950_impact.xlsx", index=True, header=False)


working_model.compare_mean_replicates_table("ATP")
working_model.compare_mean_replicates_table("MSU")
working_model.compare_mean_replicates_table("Nigericin")
# plt, ax = working_model.line_plot_selected()
# plt.show()
##COMPARE DIFF BETWEEN CURVES
working_model.quantify_curve_difference("MSU", "ATP")

##CHECK IF THERE IS ANY DIFFERENCE IN GROWTH ACROSS ALL TIMES
working_model.compare_mean_replicates_table("ATP")


##CHECK GROWTH RATE OF EACH PAIR
working_model.compare_growth_rate()
##BUNCH OF STUFF
working_model.create_summary_table()
##PLOT TIME COURSE ANALYSIS
working_model.plot_time_course_analysis("ATP")
working_model.plot_time_course_analysis("ATP", "MCC950_ATP")
working_model.plot_time_course_analysis("Nigericin", "MCC950_Nigericin")
working_model.plot_time_course_analysis("MSU", "MCC950_MSU")
working_model.plot_time_course_analysis("ATP", "Nigericin", points_only=True)
##CHECKS ALL REPLICATES OF A TREATMENT
working_model.compare_replicate_curves("MSU")

working_model.line_plot()
plt, ax = working_model.line_plot()
plt.show()


treatment = "MCC950_ATP"
treatment = "MCC950_MSU"
treatment = "MCC950_Nigericin"
treatment = "ATP"
treatment = "MSU"
treatment = "Nigericin"
shift_x = True
shift_x = False
compare_standardized_replicate_curves(working_model, "ATP", shift_x=True)
