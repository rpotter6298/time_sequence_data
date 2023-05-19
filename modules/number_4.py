def plot_count_against_ratio(ratio_name, treatments=None, invert=False):
    speck_info = TAS.modules["TS_Speck"]
    ratio_info = TAS.modules["TS_Cyto"].ratio_data[ratio_name]
    ratio_name_x = ratio_name.split(":")[0]
    ratio_name_y = ratio_name.split(":")[1]
    if treatments is None:
        treatments = speck_info.comp.treatments
    if len(treatments) == 1:
        treatment = treatments[0]
    elif treatments == speck_info.comp.treatments:
        treatment = "All Treatments"
    speck_data = speck_info.data[speck_info.data["Treatment"].isin(treatments)]
    ratio_data = ratio_info[ratio_info["Treatment"].isin(treatments)]
    fig, ax = plt.subplots(figsize=(10, 6))
    ax.set_ylabel("Normalized Speck Formation")
    ax2 = ax.twinx()
    lineplot_ax1 = sns.lineplot(
        data=speck_data,
        x="Time (hrs)",
        y="Measurement",
        hue="Treatment",
        ax=ax,
        errorbar="se",
    )
    if invert == True:
        y_name = f"Ratio_{ratio_name_y}_{ratio_name_x}"
        ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_y}:{ratio_name_x}")
    else:
        y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
        ax2.set_ylabel(f"Ratio of Cytokine Measurements, {ratio_name_x}:{ratio_name_y}")
    pointplot = sns.pointplot(
        data=ratio_data,
        x="Time (hrs)",
        y=y_name,
        hue="Treatment",
        ax=ax2,
        palette="pastel",
    )
    # ax.set_ylim(1, None)  # Set the lower limit of ax's y-axis to 1
    # ax2.set_ylim(1, None)
    ax.get_legend().remove()  # Set the lower limit of ax2's y-axis to 1
    # for line in pointplot.lines:
    #     line.set_linestyle("dotted")
    plt.xlim(0, 21)
    plt.xlabel("Time (hrs)")

    plt.title(
        str(
            "Normalized Measurements for IL1b, IL18, and Normalized Speck Formation - "
            + treatment
        )
    )
    plt.show()
