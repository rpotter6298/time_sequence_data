import pandas as pd
from pathlib import Path


class Cytokine_Adapter:
    def __init__(
        self,
        treatments,
        directory=Path("data", "cytokine_data", "plates"),
        modifier=None,
        short_name_dict={},
    ):
        self.treatments = treatments
        self.modifier = modifier
        self.short_name_dict = short_name_dict
        self.directory = directory
        self.data = self.validate_treatments(self.structure_data())

    def structure_data(self):
        compiled_plates = pd.DataFrame()
        for csv in self.directory.glob("*.csv"):
            print(csv)
            data = pd.read_csv(csv, skiprows=1)
            data["Analyte"] = csv.name.split(" ")[-1].split(".")[0]
            data["Treatment"] = data["Sample"].str.split(" ", expand=True)[0]
            data["Time (hrs)"] = (
                data["Sample"]
                .str.split(" ", expand=True)[1]
                .str.split(".", expand=True)[0]
            )
            data["ExperimentalReplicate"] = (
                data["Sample"]
                .str.split(" ", expand=True)[1]
                .str.split(".", expand=True)[1]
            )
            data["Concentration"] = data["Calc. Concentration"]
            data["Treatment_Replicate"] = (
                (data["Treatment"]) + "_" + (data["ExperimentalReplicate"])
            )
            data.drop(
                columns=[
                    "Sample",
                    "Mean",
                    "CV",
                    "Calc. Concentration",
                    "Calc. Conc. Mean",
                    "Calc. Conc. CV",
                    "Detection Range",
                ],
                inplace=True,
            )
            compiled_plates = pd.concat([compiled_plates, data])
        return compiled_plates

    def validate_treatments(self, data):
            """
            Validates the treatment names in the DataFrame against class treatments.

            Parameters:
            data (DataFrame): The DataFrame containing the treatment names.

            Raises:
            ValueError: If an unrecognized treatment name is found.
            """
            valid_treatments = self.treatments + [f"{self.modifier}_{t}" for t in self.treatments] + list(self.short_name_dict.keys())
            data = data[data["Treatment"].str.lower().isin([t.lower() for t in valid_treatments])]
            for key, value in self.short_name_dict.items():
                data.loc[data["Treatment"].str.lower() == key.lower(), "Treatment"] = value
            return(data)
    