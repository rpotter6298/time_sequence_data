from modules import *
from pathlib import Path
import numpy as np


##Example Usage:
file_path = Path("data", "count_all_data.xlsx")
base_treatments = ["ATP", "MSU", "Nigericin"]
modifier = "MCC950"
short_name_dict = {"nig": "Nigericin"}
aggregated = Count_Aggregator(
    file_path, base_treatments, modifier, short_name_dict=short_name_dict
)
# file_path = "speck_data.xlsx"
# summary = Summary_Adapter(file_path, base_treatments, modifier, short_name_dict)
ts_module = time_sequence_module(aggregated, time_limits={"ALL": 21})
working_model = ts_analysis_module(ts_module)

# plt, ax = working_model.line_plot_selected()
# plt.show()
##COMPARE DIFF BETWEEN CURVES
working_model.quantify_curve_difference("MSU", "MCC950_MSU")

##CHECK IF THERE IS ANY DIFFERENCE IN GROWTH ACROSS ALL TIMES
working_model.compare_mean_replicates_table("MCC950_ATP")
##CHECK IF THE MODIFIER CHANGES THE BEHAVIOR AND BY HOW MUCH
working_model.check_modifier_impact()
##CHECK GROWTH RATE OF EACH PAIR
working_model.compare_growth_rate()
##BUNCH OF STUFF
working_model.create_summary_table()
##PLOT TIME COURSE ANALYSIS
working_model.plot_time_course_analysis("ATP")
working_model.plot_time_course_analysis("ATP", "MCC950_ATP")
working_model.plot_time_course_analysis("Nigericin", "MCC950_Nigericin")
working_model.plot_time_course_analysis("MSU", "MCC950_MSU")
working_model.plot_time_course_analysis("ATP", "MSU", points_only=True)
##CHECKS ALL REPLICATES OF A TREATMENT
working_model.compare_replicate_curves("MSU")

working_model.line_plot()
plt, ax = working_model.line_plot()
plt.show()


working_model.compare_standardized_replicate_curves("MCC950_ATP", shift_x=True)
working_model.compare_standardized_replicate_curves("MCC950_MSU", shift_x=True)
working_model.compare_standardized_replicate_curves("MCC950_Nigericin", shift_x=True)
working_model.compare_standardized_replicate_curves("MCC950_ATP", shift_x=False)
working_model.compare_standardized_replicate_curves("MCC950_MSU", shift_x=False)
working_model.compare_standardized_replicate_curves("MCC950_Nigericin", shift_x=False)

working_model.compare_standardized_replicate_curves("ATP")
working_model.compare_standardized_replicate_curves("MSU")
working_model.compare_standardized_replicate_curves("Nigericin")

treatment = "MCC950_ATP"
treatment = "MCC950_MSU"
treatment = "MCC950_Nigericin"
treatment = "ATP"
treatment = "MSU"
treatment = "Nigericin"
shift_x = True
shift_x = False
compare_standardized_replicate_curves(working_model, "ATP", shift_x=True)
