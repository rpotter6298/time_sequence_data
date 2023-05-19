import pandas as pd
import numpy as np
from typing import Type, List, Dict, Callable
import matplotlib.pyplot as plt
import seaborn as sns
from modules.stats_functions import bootstrap_t_test


class plotting_module:
    @staticmethod
    def plot_lineplot(
        dataframe: pd.DataFrame,
        y: str,
        x: str = "Time (hrs)",
        hue: str = "Treatment",
        hue_order: List = None,
        errorbar: str = "se",
        filepath: str = None,
        show_p=False,
        manual_ax_modification=None,
    ):
        fig, ax = plt.subplots(figsize=(10, 6))
        lineplot = sns.lineplot(
            data=dataframe,
            x=x,
            y=y,
            hue=hue,
            hue_order=hue_order,
            errorbar=errorbar,
            ax=ax,
        )
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        ax.set_title(f"{y} Over Time")
        ax.legend()
        if manual_ax_modification is not None:
            manual_ax_modification(
                ax
            )  # Call the callback function with the ax parameter
        if show_p == True:
            ax = plotting_module.show_p(data=dataframe, y=y, separator=hue, ax=ax)
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
            return fig, ax

    @staticmethod
    def plot_pointplot(
        dataframe: pd.DataFrame,
        y: str,
        x: str = "Time (hrs)",
        hue: str = "Treatment",
        hue_order: List = None,
        errorbar: str = "se",
        filepath: str = None,
        show_p=False,
        manual_ax_modification=None,
    ):
        fig, ax = plt.subplots(figsize=(10, 6))
        pointplot = sns.pointplot(
            data=dataframe,
            x=x,
            y=y,
            hue=hue,
            hue_order=hue_order,
            errorbar=errorbar,
            ax=ax,
        )
        ax.set_xlabel(x)
        ax.set_ylabel(y)
        # for line in pointplot.lines:
        #     line.set_linestyle("dashed")
        ax.set_title(f"{y} Over Time")
        ax.legend()
        if manual_ax_modification is not None:
            manual_ax_modification(
                ax
            )  # Call the callback function with the ax parameter
        if show_p == True:
            ax = plotting_module.show_p(data=dataframe, y=y, separator=hue, ax=ax)
        if filepath is not None:
            plt.savefig(filepath, dpi=300, bbox_inches="tight")
            plt.close()
        else:
            plt.show()
            return fig, ax

    @staticmethod
    def plot_multimodel(
        modelA,
        modelB,
        measurement_type="Normalized_Measurement",
        output_directory=None,
    ):
        ##Confirm that ModelA and ModelB have the same treatments
        if modelA.comp.treatments != modelB.comp.treatments:
            raise ValueError("Models do not have the same treatments.")
        for treatment in modelA.comp.treatments:
            dfA = modelA.data[modelA.data["Treatment"] == treatment]
            dfB = modelB.data[modelB.data["Treatment"] == treatment]
            fig, ax = plt.subplots(figsize=(10, 6))
            ax.set_ylabel("Normalized Speck Formation")
            ax2 = ax.twinx()
            ax2.set_ylabel("Normalized Cytokine Measurement")
            for model in [dfA, dfB]:
                ax_to_use = ax if model is dfA else ax2
                if len(model["Analyte"].unique()) > 1:
                    sns.lineplot(
                        data=model,
                        x="Time (hrs)",
                        y=measurement_type,
                        hue="Analyte",
                        ax=ax_to_use,
                        errorbar="se",
                    )
                else:
                    sns.lineplot(
                        data=model,
                        x="Time (hrs)",
                        y=measurement_type,
                        color="green",
                        ax=ax_to_use,
                        errorbar="se",
                    )
            ax.set_ylim(1, None)  # Set the lower limit of ax's y-axis to 1
            ax2.set_ylim(1, None)  # Set the lower limit of ax2's y-axis to 1
            plt.xlabel("Time (hrs)")
            # plt.ylabel("Normalized Value")
            plt.title(
                str(
                    "Normalized Measurements for IL1b, IL18, and Normalized Speck Formation - "
                    + treatment
                )
            )
            if output_directory is not None:
                filepath = output_directory / f"{treatment}_multimodel.png"
                plt.savefig(filepath, dpi=300, bbox_inches="tight")
                plt.close()
            else:
                plt.show()

    @staticmethod
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
            ratio_name_x, ratio_name_y = ratio_name_y, ratio_name_x

        y_name = f"Ratio_{ratio_name_x}_{ratio_name_y}"
        
        pointplot = sns.pointplot(
            data=ratio_data,
            x="Time (hrs)",
            y=y_name,
            hue="Treatment",
            ax=ax2,
            palette="pastel",
        )
        ax2.set_ylabel(
            f"Ratio: {ratio_name_x}:{ratio_name_y}"
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
                f"Ratio of {ratio_name_x}:{ratio_name_y} and Speck Formation Counts - "
                + treatment
            )
        )
        plt.show()

    @staticmethod
    def plot_ratio_old(
        module, measurement_type="Measurement", normalize_start=False, invert=False
    ):
        ##If the module's data does not have two analytes, return an error
        if len(module.data["Analyte"].unique()) < 2:
            raise ValueError("Model does not have two analytes.")
        split_data = plotting_module.split_analytes(module)
        df1, df2 = split_data.values()
        if invert == True:
            df1, df2 = df2, df1
        analyte1, analyte2 = split_data.keys()
        merged_data = pd.merge(
            df1, df2, on=["Treatment", "Time (hrs)", "Experimental_Replicate"]
        )
        merged_data[f"Ratio_{analyte1}_{analyte2}"] = (
            merged_data[f"{measurement_type}_x"] / merged_data[f"{measurement_type}_y"]
        )
        # If normalize_start is True, normalize the ratio to the first time point for each treatment and experimental replicate
        if normalize_start:
            merged_data[f"Ratio_{analyte1}_{analyte2}"] = merged_data.groupby(
                ["Treatment", "Experimental_Replicate"]
            )[f"Ratio_{analyte1}_{analyte2}"].transform(lambda x: x / x.iloc[0])

        plotting_module.plot_lineplot(
            merged_data, y=f"Ratio_{analyte1}_{analyte2}", errorbar="se"
        )
        return merged_data

    @staticmethod
    def plot_ratio(
        module,
        invert=False,
        analyte_pair=None,
        manual_ax_modification=None,
    ):
        ##If module does not have .ratio_data, return an error
        if module.ratio_data is None:
            raise ValueError("Module does not have ratio_data.")
        ##If analyte_pair is None, use the first key in the ratio_data dictionary
        if analyte_pair is None:
            analyte_pair = list(module.ratio_data.keys())[0]
        merged_data = module.ratio_data[analyte_pair]
        if invert == True:
            analyte1 = analyte_pair.split(":")[1]
            analyte2 = analyte_pair.split(":")[0]
        else:
            analyte1 = analyte_pair.split(":")[0]
            analyte2 = analyte_pair.split(":")[1]
        fig, ax = plotting_module.plot_pointplot(
            merged_data, y=f"Ratio_{analyte1}_{analyte2}", errorbar="se", manual_ax_modification=manual_ax_modification
        )
        #return fig, ax

    @staticmethod
    def split_analytes(module):
        """
        Split the data into separate dataframes for each analyte.

        Args:
            model (time_sequence): The time_sequence object to split.

        Returns:
            dict: A dictionary of dataframes for each analyte.
        """
        analytes = module.data["Analyte"].unique()
        split_data = {}
        for analyte in analytes:
            split_data[analyte] = module.data[module.data["Analyte"] == analyte]
        return split_data

    @staticmethod
    def show_p(data, y, separator, ax):
        # Get the current ylim values
        ymin, ymax = ax.get_ylim()
        # Set the new ylim with an increased upper limit
        ax.set_ylim(ymin, ymax * 1.25)
        # Get the treatments in the data
        treatments = data[separator].unique()
        p_values_text = ""
        for treatment in treatments:
            treatment_data = data[data[separator] == treatment]
            max_time = treatment_data.loc[treatment_data[y].idxmax(), "Time (hrs)"]

            # max_time_per_replicate = (
            #     treatment_data.groupby("Experimental_Replicate")[y]
            #     .idxmax()
            #     .reset_index()
            # )
            # # Merge the original treatment_data with max_time_per_replicate to get the values for each replicate at their max times
            # max_values_per_replicate = treatment_data.loc[max_time_per_replicate[y]]
            # # Filter the max_values_per_replicate DataFrame for the current treatment
            # treatment_max_values = max_values_per_replicate[
            #     max_values_per_replicate[separator] == treatment
            # ]

            time_zero_data = np.array(
                treatment_data[treatment_data["Time (hrs)"] == 0][y]
            )
            max_time_data = np.array(
                treatment_data[treatment_data["Time (hrs)"] == max_time][y]
            )
            # print(time_zero_data, max_time_data)
            p, ci = bootstrap_t_test(max_time_data, time_zero_data)
            # print(p)
            p_values_text += f"{treatment} (Max = {max_time} hrs):\n p = {p:.4f}\n"

        ax.text(
            0.02,
            0.98,
            p_values_text,
            transform=ax.transAxes,
            fontsize=10,
            verticalalignment="top",
            horizontalalignment="left",
        )
        return ax
