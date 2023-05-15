import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

cytos = tsa
specks = model

specks.data
cytos.data[(cytos.data["Analyte"] == "IL18") & (cytos.data["Treatment"] == "MSU")][
    "Concentration"
]
cytos.data[(cytos.data["Analyte"] == "IL1b") & (cytos.data["Treatment"] == "MSU")][
    "Concentration"
]


for treatment in ["ATP", "MSU", "Nigericin"]:
    cyto_df = cytos.data[cytos.data["Treatment"] == treatment]
    speck_df = specks.data[specks.data["Treatment"] == treatment]

    fig, ax = plt.subplots(figsize=(10, 6))
    sns.lineplot(
        data=cyto_df,
        x="Time (hrs)",
        y="Normalized_Measurement",
        hue="Analyte",
        errorbar="se",
    )
    sns.lineplot(
        data=speck_df,
        x="Time (hrs)",
        y="NormalizedSpeckFormation",
        color="Green",
        errorbar="se",
    )

    plt.xlabel("Time (hrs)")
    plt.ylabel("Normalized Value")
    plt.title(
        str(
            "Normalized Measurements for IL1b, IL18, and Normalized Speck Formation - "
            + treatment
        )
    )

    handles, labels = ax.get_legend_handles_labels()

    # Add white and black lines for IL18 and IL1b
    from matplotlib.lines import Line2D

    handles.append(Line2D([0], [0], color="green", linewidth=2))
    labels.append("Speck Formation")

    ax.legend(handles=handles, labels=labels)
    fig.show()
