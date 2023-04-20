from modules import *


##Example Usage:
file_path = "count_all_data.xlsx"
base_treatments = ["ATP", "MSU", "Nigericin"]
modifier = "MCC950"
short_name_dict = {"nig": "Nigericin"}
aggregated = Count_Aggregator(
    file_path, base_treatments, modifier, short_name_dict=short_name_dict
)
file_path = "speck_data.xlsx"
summary = Summary_Adapter(file_path, base_treatments, modifier, short_name_dict)

working_model = ts_analysis_module(time_sequence_module(aggregated))
