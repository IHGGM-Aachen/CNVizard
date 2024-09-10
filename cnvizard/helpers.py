"""
File which contains little helpers of the CNVizard
@author: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd


def filter_tsv(
    tsv: pd.DataFrame,
    chromosome_list_cnv: list,
    cnv_type: list,
    acmg_class: list,
    entered_cnv_chrom: list,
    entered_cnv_type: list,
    entered_acmg_class: list,
):
    """
    Function which applies previously defined filters on the .tsv DataFrame.

    Args:
        tsv (pd.DataFrame): .tsv DataFrame
        chromosome_list_cnv (list): List which contains all the possible chromosomes (used to negate an empty filter)
        cnv_type (list): List which contains all the possible cnv_types (used to negate an empty filter)
        acmg_class (list): List which contains all the possible acmg_classes (used to negate an empty filter)
        entered_cnv_chrom (list): List which contains the chromosomes selected by the user
        entered_cnv_type (list): List which contains the cnv types selected by the user
        entered_acmg_class (list): List which contains the acmg classes selected by the user

    Returns:
        pd.DataFrame: Filtered .tsv DataFrame
    """
    # Perform casts to ensure the dataframe columns contain the correct types
    # This cast is performed to ensure compatibility with different AnnotSV
    # versions
    tsv["ACMG_class"] = tsv["ACMG_class"].fillna("NA")
    tsv["ACMG_class"] = tsv["ACMG_class"].astype(str)
    tsv_with_equal = tsv[(tsv["ACMG_class"].str.contains("="))]
    tsv_with_equal["ACMG_class"] = tsv_with_equal["ACMG_class"].str.split("=").str[1]
    tsv_without_equal = tsv[~(tsv["ACMG_class"].str.contains("="))]
    tsv = pd.concat([tsv_without_equal,tsv_with_equal]).reset_index(drop=True)
    tsv["SV_chrom"] = tsv["SV_chrom"].astype(str)
    tsv["SV_type"] = tsv["SV_type"].astype(str)
    # If entered_cnv_chrom is empty the empty filter is negated by assining
    # chromosome_list_cnv
    if entered_cnv_chrom is None or entered_cnv_chrom == []:
        entered_cnv_chrom = chromosome_list_cnv
    # If entered_cnv_type is empty the empty filter is negated by assining
    # cnv_type
    if entered_cnv_type is None or entered_cnv_type == []:
        entered_cnv_type = cnv_type
    # If entered_acmg_class is empty the empty filter is negated by assining
    # acmg_class
    if entered_acmg_class is None or entered_acmg_class == []:
        entered_acmg_class = acmg_class
    # Apply filters
    filtered_tsv = tsv[tsv["SV_chrom"].isin(entered_cnv_chrom)]
    filtered_tsv = filtered_tsv[filtered_tsv["SV_type"].isin(entered_cnv_type)]
    filtered_tsv = filtered_tsv[filtered_tsv["ACMG_class"].isin(entered_acmg_class)]

    return filtered_tsv
