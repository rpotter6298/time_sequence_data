import pandas as pd


class time_sequence_module:
    """
    A class for processing and analyzing time sequence data from a data_class object.
    """

    def __init__(self, data_class, time_limits: dict = None):
        """
        Initializes the time_sequence_module object.

        :param data_class: A data_class object containing the raw data to be analyzed.
        :param time_limits: An optional dictionary specifying time limits for the data. Defaults to None.
        """
        self.comp = data_class
        self.raw_data = data_class.data
        self.time_limits = time_limits
        self.data = self.normalize()

    def normalize(self, robust=True):
        """
        Normalizes the data by dividing the SpeckFormation by the initial values at time 0.

        :param robust: Optional boolean to calculate growth rate. Defaults to True.
        :return: DataFrame with normalized data.
        """
        tech_merge = input("Merge Technical Replicates? (y/n): ")
        if tech_merge == "y":
            data = self.merge_technical_replicates()
        elif tech_merge == "n":
            data = self.preserve_technical_replicates()
        else:
            print("Invalid input.")
            return
        initial_values = data[data["Time (hrs)"] == 0].set_index("Treatment_Replicate")[
            "SpeckFormation"
        ]
        data["NormalizedSpeckFormation"] = data.apply(
            lambda row: row["SpeckFormation"]
            / initial_values[row["Treatment_Replicate"]],
            axis=1,
        )
        if self.time_limits is not None:
            data = self.limit_data(data, self.time_limits)
        if robust == True:
            data = self.calc_growth_rate(data)
        return data

    def limit_data(self, data, time_limits):
        """
        Limits the data based on the provided time limits.

        :param data: DataFrame containing the data to be limited.
        :param time_limits: Dictionary specifying the time limits for the data.
        :return: DataFrame with limited data.
        """
        limited_data = data
        if time_limits is None:
            print("No time limits specified. Returning all data.")
            return limited_data
        else:
            for treatment in time_limits.keys():
                if treatment == "ALL":
                    limited_data = limited_data[
                        limited_data["Time (hrs)"] <= time_limits[treatment]
                    ]
                else:
                    limited_data = limited_data[
                        (limited_data["Treatment"].str.contains(treatment, case=False))
                        & (limited_data["Time (hrs)"] <= time_limits[treatment])
                        | (
                            ~limited_data["Treatment"].str.contains(
                                treatment, case=False
                            )
                        )
                    ]
            return limited_data

    def merge_technical_replicates(self):
        """
        Merges technical replicates by averaging their SpeckFormation values.

        :return: DataFrame with merged technical replicates.
        """
        data = self.raw_data.copy()
        data["Treatment_Replicate"] = (
            data["Treatment"] + "_" + data["ExperimentalReplicate"]
        )
        merged_data = (
            data.groupby(["Treatment_Replicate", "Time (hrs)"])["SpeckFormation"]
            .mean()
            .reset_index()
        )
        # Merge the grouped data back with the original data, dropping the "SpeckFormation" column
        data = data.drop("SpeckFormation", axis=1).merge(
            merged_data, on=["Treatment_Replicate", "Time (hrs)"]
        )
        data["TechnicalReplicate"] = pd.to_numeric(
            data["TechnicalReplicate"], errors="coerce"
        )
        data = data[data["TechnicalReplicate"] <= 1]
        # Set the remaining "TechnicalReplicate" values to 0
        data["TechnicalReplicate"] = 0
        # Assign the merged data back to self.raw_data
        return data

    def preserve_technical_replicates(self):
        """
        Preserves the technical replicates in the data.

        :return: DataFrame with preserved technical replicates.
        """
        data = self.raw_data.copy()
        data["Treatment_Replicate"] = (
            data["Treatment"]
            + "_"
            + data["ExperimentalReplicate"].astype(str)
            + "_"
            + data["TechnicalReplicate"].astype(str)
        )
        return data

    def calc_growth_rate(self, data):
        """
        Calculates the growth rate by taking the difference in NormalizedSpeckFormation between consecutive time points.

        :param data: DataFrame containing the normalized data.
        :return: DataFrame with growth rate data.
        """
        growth_rate = data.copy()
        growth_rate["GrowthRate"] = growth_rate.groupby(
            ["Treatment", "Treatment_Replicate"]
        )["NormalizedSpeckFormation"].diff()
        data["GrowthRate"] = growth_rate["GrowthRate"]
        return data
