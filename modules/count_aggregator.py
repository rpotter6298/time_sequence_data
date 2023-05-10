import pandas as pd
import re
from pathlib import Path


class Count_Aggregator:
    def __init__(
        self,
        name,
        treatments,
        modifier=None,
        short_name_dict={},
        file_path=Path("data", "speck_data", "count_all_data.xlsx"),
    ):
        """
        Initializes the count_aggregator class.

        Parameters:
        file_path (str): Path to the Excel file containing the data.
        treatments (list): A list of treatment names.
        modifier (str, optional): The treatment modifier. Default is None.
        short_name_dict (dict, optional): A dictionary of short names for treatments. Default is an empty dictionary.
        """
        self.name = name
        self.xls = pd.read_excel(file_path, sheet_name=None, header=1)
        self.treatments = treatments
        self.modifier = modifier
        self.short_name_dict = short_name_dict
        self.dflist = self.build_df()
        self.data = self.aggregate_data()

    def match_treatment(self, sheet_name):
        """
        Matches the treatment based on the sheet name.

        Parameters:
        sheet_name (str): Name of the sheet in the Excel file.

        Returns:
        matched_treatment (str): The matched treatment name.
        """
        sheet_name_lower = sheet_name.lower().replace("_", "").replace(" ", "")
        for short_name, full_name in self.short_name_dict.items():
            if short_name.lower() in sheet_name_lower:
                sheet_name_lower = sheet_name_lower.replace(
                    short_name.lower(), full_name.lower()
                )
        matched_treatment = ""
        for treatment in self.treatments:
            if treatment.lower() in sheet_name_lower:
                matched_treatment = treatment
                break
        if self.modifier.lower() in sheet_name_lower:
            matched_treatment = f"{self.modifier}_{matched_treatment}"
        return matched_treatment

    def list_short_names(self, query: str):
        """
        Lists short names for a given treatment.

        Parameters:
        query (str): The treatment name to find short names for.

        Returns:
        key_match (list): A list of short names matching the treatment.
        """
        key_match = []
        for key, value in self.short_name_dict.items():
            if value == query:
                key_match.append(key)
        return key_match

    def extract_replicates(self, sheet_name):
        """
        Extracts experimental and technical replicates from the sheet name.

        Parameters:
        sheet_name (str): Name of the sheet in the Excel file.

        Returns:
        experimental_replicate (str or None): The experimental replicate number or None if not found.
        technical_replicate (str or None): The technical replicate number or None if not found.
        """
        sheet_name_lower = sheet_name.lower().replace("_", "").replace(" ", "")
        for treatment in self.treatments:
            if treatment.lower() in sheet_name_lower:
                sheet_name_lower = sheet_name_lower.replace(treatment.lower(), "")
                break
            elif treatment in self.short_name_dict.values():
                qlist = self.list_short_names(treatment)
                for q in qlist:
                    if q.lower() in sheet_name_lower:
                        sheet_name_lower = sheet_name_lower.replace(q.lower(), "")
                        break
        if self.modifier.lower() in sheet_name_lower:
            sheet_name_lower = sheet_name_lower.replace(self.modifier.lower(), "")

        decimal_match = re.search(r"(?<!\d)(\d+)\.?(\d+)?", sheet_name_lower)
        if decimal_match:
            experimental_replicate = decimal_match.group(1)
            technical_replicate = (
                decimal_match.group(2) if decimal_match.group(2) else 1
            )
        else:
            experimental_replicate, technical_replicate = None, None

        return experimental_replicate, technical_replicate

    def build_df(self):
        """
        Builds a DataFrame from the Excel file sheets.

        Returns:
        all_sheets (list): A list of DataFrames containing the data from each sheet.
        """
        all_sheets = []
        for sheet_name, df in self.xls.items():
            # print(sheet_name)
            # Add a column with the sheet name
            df["SheetName"] = sheet_name

            # Create a column for treatment
            matched_treatment = self.match_treatment(sheet_name)
            df["Treatment"] = matched_treatment

            # Create a column for experimental replicate
            experimental_replicate, technical_replicate = self.extract_replicates(
                sheet_name
            )
            df["Experimental_Replicate"] = experimental_replicate
            df["Technical_Replicate"] = technical_replicate
            all_sheets.append(df)
        return all_sheets

    def aggregate_data(self):
        """
        Aggregates the data from all sheets into a single DataFrame compatible with time_sequence_module.

        Returns:
        aggregate_data (DataFrame): The aggregated data from all sheets.
        """
        aggregate_data = pd.DataFrame()
        for df in self.dflist:
            addition = df.filter(
                [
                    "Time (hrs)",
                    "Treatment",
                    "Experimental_Replicate",
                    "Technical_Replicate",
                ],
                axis=1,
            )
            pos_columns = df.filter(like="Pos", axis=1)
            addition["SpeckFormation"] = pos_columns.sum(axis=1)
            aggregate_data = pd.concat([aggregate_data, addition])
        return aggregate_data


# ##Example Usage:
# file_path = "count_all_data.xlsx"
# base_treatments = ["ATP", "MSU", "Nigericin"]
# modifier = "MCC950"
# short_name_dict = {"nig": "Nigericin"}
# structure = count_aggregator(
#     file_path, base_treatments, modifier, short_name_dict=short_name_dict
# )
# data = structure.data
