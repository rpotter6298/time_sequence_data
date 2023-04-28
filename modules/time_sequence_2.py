import pandas as pd


class time_sequence_module:
    """
    A class for processing and analyzing time sequence data from a data_class object.
    """

    def __init__(self, data_module, variable, time_limits: dict = None):
        """
        Initializes the time_sequence_module object.

        :param data_module: A data_class object containing the raw data to be analyzed.
        :param time_limits: An optional dictionary specifying time limits for the data. Defaults to None.
        """
        self.comp = data_module
        self.raw_data = data_module.data
        self.time_limits = time_limits
        self.data = self.normalize(self.raw_data, variable)

    def normalize(self, data, variable, robust=True):
        """
        Normalizes the data by dividing the SpeckFormation by the initial values at time 0.
        :param data: The dataframe containing the data to be normalized.
        :param variable: The variable to be normalized. Should be a column name in the data.
        :param robust: Optional boolean to calculate growth rate. Defaults to True.
        :return: DataFrame with normalized data.
        """

        # if self.raw_data.columns.contains("TechnicalReplicate"):
        #     tech_merge = input("Merge Technical Replicates? (y/n): ")
        #     if tech_merge == "y":
        #         data = self.merge_technical_replicates()
        #     elif tech_merge == "n":
        #         data = self.preserve_technical_replicates()
        #     else:
        #         print("Invalid input.")
        #         return
        def norm(data, variable):
            initial_values = data[data["Time (hrs)"].astype(int) == 0].set_index(
                "Treatment_Replicate"
            )[variable]
            data["NormalizedMeasurement"] = data.apply(
                lambda row: row[variable] / initial_values[row["Treatment_Replicate"]],
                axis=1,
            )
            return data

        if data["Analyte"].nunique() > 1:
            fullset = pd.DataFrame()
            for analyte in data["Analyte"].unique():
                subdata = norm(data[data["Analyte"] == analyte], variable=variable)
                fullset = pd.concat([fullset, subdata])
            data = fullset
        else:
            data = norm(data, variable=variable)

        if self.time_limits is not None:
            data = self.limit_data(data, self.time_limits)
        if robust == True:
            data = self.calc_change_rate(data)
        return data

    def calc_change_rate(self, data):
        """
        Calculates the growth rate by taking the difference in NormalizedSpeckFormation between consecutive time points.

        :param data: DataFrame containing the normalized data.
        :return: DataFrame with growth rate data.
        """
        growth_rate = data.copy()
        growth_rate["RateofChange"] = growth_rate.groupby(
            ["Treatment", "Treatment_Replicate"]
        )["NormalizedMeasurement"].diff()
        data["RateofChange"] = growth_rate["RateofChange"]
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
        else:
            for treatment in time_limits.keys():
                limited_data = limited_data[
                    (limited_data["Treatment"] != treatment)
                    | (
                        limited_data["Time (hrs)"].astype(float)
                        <= time_limits[treatment]
                    )
                ]
        return limited_data
