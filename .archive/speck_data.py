from modules import *
from pathlib import Path


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


model = main()
impact = model.check_modifier_impact()
impact.T.to_excel("MCC950_impact.xlsx", index=True, header=False)

summary = model.create_summary_table()
summary.to_excel("Speck_Summary.xlsx", index=False, header=True)

plt, ax = model.line_plot_selected(["ATP", "MCC950_ATP"], show_p=True)
plt.savefig("ATP.png", dpi=300)
plt, ax = model.line_plot_selected(["MSU", "MCC950_MSU"], show_p=True)
plt.savefig("MSU.png", dpi=300)
plt, ax = model.line_plot_selected(["Nigericin", "MCC950_Nigericin"], show_p=True)
plt.savefig("Nigericin.png", dpi=300)

model.
model.compare_time_sequence_count_to_zero_time("ATP")

[(model.data["Treatment"] == "MSU") & (model.data["ExperimentalReplicate"] == "3")]
model.compare_replicate_curves("MSU")
model.check_modifier_impact().T
model.compare_mean_replicates_table("ATP")

plt, ax = model.line_plot_selected(["ATP", "MSU"], show_p=True)
plt.show()
model.max_growth_rate
plt, ax = model.line_plot()
plt.show()
