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
cyto.data = cyto.limit_data(cyto.data, time_limits=time_limits)

speck = time_sequence_module(
    Count_Aggregator(
        "Speck", base_treatments, modifier=modifier, short_name_dict=short_name_dict
    )
)
speck.structure_data("SpeckFormation")
speck.data = speck.limit_data(speck.data, time_limits=time_limits)

TAS = analysis_module([cyto, speck])

## Create Summaries
for module in TAS.modules:
    TAS.create_summary_table(TAS.modules[module]).to_excel(
        Path("reports", "".join(["summary_", module, ".xlsx"])), index=False
    )
for module in [module for module in TAS.modules if TAS.modules[module].comp.modifier]:
    TAS.check_modifier_impact(TAS.modules[module]).T.to_excel(
        Path("reports", "".join(["modifier_impact_", module, ".xlsx"])), header=False
    )


##Line Plots of Speck Data With and Without Modifiers
module = TAS.modules["TS_Speck"]
for treatment in module.data["Treatment"].unique():
    TAS.plot_speck_lineplots(
        module, [treatment], output_path=Path("plots", "speck_plots"), show_p=True
    )
for treatment in module.comp.treatments:
    mod_treat = "_".join([module.comp.modifier, treatment])
    TAS.plot_speck_lineplots(
        module,
        [treatment, mod_treat],
        output_path=Path("plots", "speck_plots", "modifier_comparison"),
        show_p=True,
    )

##Line Plots of Cytokine Data
module = TAS.modules["TS_Cyto"]
for treatment in module.data["Treatment"].unique():
    TAS.plot_speck_lineplots(
        module, [treatment], output_path=Path("plots", "cytokine_plots"), show_p=True
    )

##Multimodel Plots
TAS.plot_multimodel(
    TAS.modules["TS_Speck"],
    TAS.modules["TS_Cyto"],
    output_directory=Path("plots", "multimodel_plots"),
)

# Arbitrary Maximums
for module in TAS.modules:
    TAS.find_arbitrary_max(TAS.modules[module])
    nested_dict = {}
    for key, value in TAS.modules[module].arbitrary_maximums.__dict__.items():
        value_dict = value.__dict__.copy()
        if "window" in value_dict:
            del value_dict["window"]
        nested_dict[key] = value_dict
    with open(Path("reports", "arbitrary_max_" + module + ".json"), "w") as outfile:
        json.dump(nested_dict, outfile)

## Is there a difference between the modifier and the base treatment at the arbitrary maximum time?
module = TAS.modules["TS_Speck"]
arb_max_diff = {}
for treatment in module.comp.treatments:
    mod_treat = "_".join([module.comp.modifier, treatment])
    arbitrary_max_time_treat = module.arbitrary_maximums.__dict__[treatment].time
    arbitrary_max_time_mod_treat = module.arbitrary_maximums.__dict__[mod_treat].time

    # The measurements for the treatment, no mod, at the arbitrary maximum for the treatment without the modifier (A1)
    A1 = module.data[
        (module.data["Treatment"] == treatment)
        & (module.data["Time (hrs)"] == arbitrary_max_time_treat)
    ]["Measurement"]
    # The measurements for the treatment, no mod, at the time of the arbitrary maximum for the treatment with the modifier (A2)
    A2 = module.data[
        (module.data["Treatment"] == treatment)
        & (module.data["Time (hrs)"] == arbitrary_max_time_mod_treat)
    ]["Measurement"]
    # The measurements for the treatment, with mod, at the time of the arbitrary maximum for the treatment without the modifier (B1)
    B1 = module.data[
        (module.data["Treatment"] == mod_treat)
        & (module.data["Time (hrs)"] == arbitrary_max_time_treat)
    ]["Measurement"]
    # The measurements for the treatment, with mod, at the arbitrary maximum for the treatment with the modifier (B2)
    B2 = module.data[
        (module.data["Treatment"] == mod_treat)
        & (module.data["Time (hrs)"] == arbitrary_max_time_mod_treat)
    ]["Measurement"]
    p, ci = bootstrap_t_test(A1, B1)
    arb_max_diff[treatment] = {
        "Time1": arbitrary_max_time_treat,
        "Time2": arbitrary_max_time_mod_treat,
        "A1-B1": {
            "p": p,
            "ci": ci,
            "mean(treatment)": round(A1.mean(), 3),
            "mean(modified_treatment)": round(B1.mean(), 3),
        },
    }
    if arbitrary_max_time_treat == arbitrary_max_time_mod_treat:
        continue
    p, ci = bootstrap_t_test(A2, B2)
    arb_max_diff[treatment]["A2-B2"] = {
        "p": p,
        "ci": ci,
        "mean(treatment)": round(A2.mean(), 3),
        "mean(modified_treatment)": round(B2.mean(), 3),
    }
    p, ci = bootstrap_t_test(A1, B2)
    arb_max_diff[treatment]["A1-B2"] = {
        "p": p,
        "ci": ci,
        "mean(treatment)": round(A1.mean(), 3),
        "mean(modified_treatment)": round(B2.mean(), 3),
    }


# Test if there is a difference between the times of the max speck and max cyto for each treatment
report = TAS.compare_max_time_distances(TAS.modules["TS_Speck"], TAS.modules["TS_Cyto"])
pprint.pprint(report)
# Test if the cytokine FC max is different between treatments
report = TAS.compare_max_normalized_measurement_anova(TAS.modules["TS_Cyto"].data)
pprint.pprint(report)
for analyte in report.keys():
    print(analyte)
    df = pd.DataFrame(report[analyte]["post_hoc"])
    print(df)
