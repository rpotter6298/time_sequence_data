from cytokine_modules import *
from pathlib import Path

base_treatments = ["ATP", "MSU", "Nigericin"]
short_name_dict = {"nig": "Nigericin"}
data_module = Cytokine_Adapter(base_treatments, short_name_dict=short_name_dict)
ts = cyt_time_sequence_module(data_module=data_module, variable="Concentration")
tsa = cyt_ts_analysis_module(ts)
tsa.aggregate_ctsc0t()

tsa.split_line_plot()
