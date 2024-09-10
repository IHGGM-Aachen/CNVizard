"""
File which contains export functions of the CNVizard
@author: Jeremias Krause , Carlos Classen, Matthias Begemann. Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd
from io import BytesIO


class CNVExporter:
    """
    Class used to export the filtered dataframe as excel tables.
    """

    def __init__(self):
        """
        Constructor of the class CNVExporter.
        """

    def save_tables_as_excel(
        self,
        total_df: pd.DataFrame,
        bintest_df: pd.DataFrame,
        hom_del_df: pd.DataFrame,
        total_candi_df: pd.DataFrame,
        bintest_candi_df: pd.DataFrame,
        consecutive_del_df: pd.DataFrame,
        consecutive_dup_df: pd.DataFrame,
    ) -> bytes:
        """
        Function which writes the preset DataFrames to an excel file.

        Each preset represents an excel sheet.

        Args:
            total_df (pd.DataFrame): .cnr DataFrame
            bintest_df (pd.DataFrame): bintest DataFrame
            hom_del_df (pd.DataFrame): .cnr DataFrame filtered for homozygous deletions
            total_candi_df (pd.DataFrame): .cnr DataFrame filtered for candigenes
            bintest_candi_df (pd.DataFrame): bintest DataFrame filtered for candigenes
            consecutive_del_df (pd.DataFrame): .cnr DataFrame filtered for consecutive deletions
            consecutive_dup_df (pd.DataFrame): .cnr DataFrame filtered for consecutive duplications

        Returns:
            bytes: Processed data to be passed to the streamlit download button object.
        """
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine="xlsxwriter")
        list_of_saved_results = [
            "total",
            "bintest",
            "hom_del",
            "total_candi",
            "bintest_candi",
            "consecutive_del",
            "consecutive_dup",
        ]
        list_of_selections = [
            total_df,
            bintest_df,
            hom_del_df,
            total_candi_df,
            bintest_candi_df,
            consecutive_del_df,
            consecutive_dup_df,
        ]

        for name, df in zip(list_of_saved_results, list_of_selections):
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            num_columns = len(df.columns.tolist())
            workbook = writer.book
            worksheet = writer.sheets[name]
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "fg_color": "#D7E4BC",
                    "border": 1,
                }
            )
            format_yellow = workbook.add_format({"bg_color": "#FFFF00"})
            format_red = workbook.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, num_columns, 11)
            worksheet.conditional_format(
                "G1:G1048576",
                {
                    "type": "3_color_scale",
                    "min_type": "num",
                    "min_value": 0,
                    "mid_type": "num",
                    "mid_value": 0.5,
                    "max_type": "num",
                    "max_value": 1,
                },
            )
            worksheet.conditional_format(
                "I1:I1048576",
                {
                    "type": "cell",
                    "criteria": "less than or equal to",
                    "value": -0.65,
                    "format": format_yellow,
                },
            )
            column_length = str(len(df))
            area_to_color = "F1:" + "F" + column_length
            worksheet.conditional_format(
                area_to_color,
                {
                    "type": "cell",
                    "criteria": "equal to",
                    "value": 0,
                    "format": format_red,
                },
            )
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            bg_format1 = workbook.add_format({"bg_color": "#EEEEEE"})
            bg_format2 = workbook.add_format({"bg_color": "#FFFFFF"})
            for rn in range(len(df)):
                worksheet.set_row(
                    rn, cell_format=(bg_format1 if rn % 2 == 0 else bg_format2)
                )

        writer.close()
        processed_data = output.getvalue()
        return processed_data

    def save_filtered_table_as_excel(
        self, filtered_df: pd.DataFrame, filtered_name: str
    ) -> bytes:
        """
        Function which writes the filtered DataFrame to an excel file.

        Args:
            filtered_df (pd.DataFrame): CNV DataFrame which was filtered with the filter criteria provided by the user.
            filtered_name (str): Name which will be used to name the excel sheet.

        Returns:
            bytes: Processed data to be passed to the streamlit download button object.
        """
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine="xlsxwriter")
        list_of_saved_results = [filtered_name]
        list_of_selections = [filtered_df]

        for name, df in zip(list_of_saved_results, list_of_selections):
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            num_columns = len(df.columns.tolist())
            workbook = writer.book
            worksheet = writer.sheets[name]
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "fg_color": "#D7E4BC",
                    "border": 1,
                }
            )
            format_yellow = workbook.add_format({"bg_color": "#FFFF00"})
            format_red = workbook.add_format({"bg_color": "#ff0000"})
            worksheet.set_column(1, num_columns, 11)
            worksheet.conditional_format(
                "G1:G1048576",
                {
                    "type": "3_color_scale",
                    "min_type": "num",
                    "min_value": 0,
                    "mid_type": "num",
                    "mid_value": 0.5,
                    "max_type": "num",
                    "max_value": 1,
                },
            )
            worksheet.conditional_format(
                "I1:I1048576",
                {
                    "type": "cell",
                    "criteria": "less than or equal to",
                    "value": -0.65,
                    "format": format_yellow,
                },
            )
            column_length = str(len(df))
            area_to_color = "F1:" + "F" + column_length
            worksheet.conditional_format(
                area_to_color,
                {
                    "type": "cell",
                    "criteria": "equal to",
                    "value": 0,
                    "format": format_red,
                },
            )
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            bg_format1 = workbook.add_format({"bg_color": "#EEEEEE"})
            bg_format2 = workbook.add_format({"bg_color": "#FFFFFF"})
            for rn in range(len(df)):
                worksheet.set_row(
                    rn, cell_format=(bg_format1 if rn % 2 == 0 else bg_format2)
                )

        writer.close()
        processed_data = output.getvalue()
        return processed_data

    def save_tables_as_excel_tsv(self, filtered_tsv: pd.DataFrame) -> bytes:
        """
        Function which writes the preset DataFrames to an excel file. Each preset represents an excel sheet.

        Args:
            filtered_tsv (pd.DataFrame): .tsv file annotated by AnnotSV, subsequently filtered and formatted.

        Returns:
            bytes: Processed data to be passed to the streamlit download button object.
        """
        output = BytesIO()
        writer = pd.ExcelWriter(output, engine="xlsxwriter")
        list_of_saved_results = ["filtered"]
        list_of_selections = [filtered_tsv]

        for name, df in zip(list_of_saved_results, list_of_selections):
            df.to_excel(writer, sheet_name=name, header=False, index=False, startrow=1)
            num_columns = len(df.columns.tolist())
            workbook = writer.book
            worksheet = writer.sheets[name]
            header_format = workbook.add_format(
                {
                    "bold": True,
                    "text_wrap": True,
                    "valign": "top",
                    "fg_color": "#D7E4BC",
                    "border": 1,
                }
            )
            worksheet.set_column(1, num_columns, 11)
            for col_num, value in enumerate(df.columns.values):
                worksheet.write(0, col_num, value, header_format)
            bg_format1 = workbook.add_format({"bg_color": "#EEEEEE"})
            bg_format2 = workbook.add_format({"bg_color": "#FFFFFF"})
            for rn in range(len(df)):
                worksheet.set_row(
                    rn, cell_format=(bg_format1 if rn % 2 == 0 else bg_format2)
                )

        writer.close()
        processed_data = output.getvalue()
        return processed_data
