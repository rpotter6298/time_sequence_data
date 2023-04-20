import pandas as pd


class Summary_Adapter:
    def __init__(self, file_path, treatments, modifier=None, short_name_dict={}):
        """
        Initializes the summary_adapter class.

        Parameters:
        file_path (str): Path to the Excel file containing the data.
        treatments (list): A list of treatment names.
        modifier (str, optional): The treatment modifier. Default is None.
        short_name_dict (dict, optional): A dictionary of short names for treatments. Default is an empty dictionary.
        """

        self.xls = pd.read_excel(file_path, header=0)
        self.treatments = treatments
        self.modifier = modifier
        self.short_name_dict = short_name_dict
        self.data = self.restructure()

    def restructure(self):
        """
        Restructures the data from the Excel file into a DataFrame compatible with time_sequence_module.

        Returns:
        data (DataFrame): The restructured data from the Excel file.
        """
        data = self.xls.set_index("Time").stack().reset_index()
        data.columns = ["Time (hrs)", "Treatment_Replicate", "SpeckFormation"]
        data["Treatment"] = data["Treatment_Replicate"].str.rsplit("_", n=1).str[0]
        data["ExperimentalReplicate"] = (
            data["Treatment_Replicate"].str.rsplit("_", n=1).str[1]
        )
        data["TechnicalReplicate"] = 0
        data.drop("Treatment_Replicate", axis=1, inplace=True)
        self.validate_treatments(data)
        return data

    def validate_treatments(self, data):
        """
        Validates the treatment names in the DataFrame against class treatments.

        Parameters:
        data (DataFrame): The DataFrame containing the treatment names.

        Raises:
        ValueError: If an unrecognized treatment name is found.
        """
        for treatment in data["Treatment"].unique():
            if (
                treatment not in self.treatments
                and treatment not in [f"{self.modifier}_{t}" for t in self.treatments]
                and treatment not in self.short_name_dict.values()
            ):
                raise ValueError(f"Unrecognized treatment: {treatment}")


# # Example usage:
# # file_path = "speck_data.xlsx"
# base_treatments = ["ATP", "MSU", "Nigericin"]
# modifier = "MCC950"
# short_name_dict = {"nig": "Nigericin"}

# adapter = summary_adapter(file_path, base_treatments, modifier, short_name_dict)
# adapter.data
