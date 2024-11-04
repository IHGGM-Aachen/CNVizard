"""
File which contains little helpers of the CNVizard
@author: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd
from cnvizard import (
    vcfMerger, 
    CNVExporter,
)
import streamlit as st
import numpy as np

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
    
def merge_annotsv_with_dbvar_study(path_to_annotsv_tsv:str,path_to_dbvar_study_file:str,path_to_output_xlsx:str):
    """
    Function to merge annotSV .tsv file with a dbvar study_file

    Args:
        path_to_annotsv_tsv (str): path to .tsv AnnotSV file
        path_to_dbvar_study_file (str): path to dvar study file vcf.gz
        path_to_output_xlsx (str): path to output xlsx

    Returns:
        pd.DataFrame: Filtered .tsv DataFrame
    """
    current_merger = vcfMerger()
    dbvar_df = current_merger.read_vcf(path_to_dbvar_study_file)
    annotsv_df = pd.read_csv(path_to_annotsv_tsv,sep="\t")
    dbvar_df["END"] = dbvar_df["INFO"].str.split(";").str[3].str.split("END=").str[1]
    dbvar_df = dbvar_df[dbvar_df["END"].notna()]
    annotsv_df["SV_start"] = annotsv_df["SV_start"].astype(int)
    annotsv_df["SV_end"] = annotsv_df["SV_end"].astype(int)
    dbvar_df["POS"] = dbvar_df["POS"].astype(int)
    dbvar_df["END"] = dbvar_df["END"].astype(int)
    annotsv_df["SV_chrom"] = annotsv_df["SV_chrom"].astype(str)
    dbvar_df["SV_chrom"] = dbvar_df["SV_chrom"].astype(str)

    list_of_overlap_tests = []
    for chrom_sign in ["1","2","3","4","5","6","7","8","9","10","11","12","13","14","15","16","17","18","19","20","21","X","Y"]:
        try:
            annotsv_df = annotsv_df[annotsv_df["SV_chrom"]==chrom_sign].reset_index(drop=True)
            dbvar_df = dbvar_df[dbvar_df["#CHROM"]==chrom_sign].reset_index(drop=True)

            annotSV_start = annotsv_df.SV_start.values
            dbvar_end = dbvar_df.END.values
            dbvar_start = dbvar_df.POS.values

            i, j = np.where((annotSV_start[:, None] >= dbvar_start) & (annotSV_start[:, None] <= dbvar_end))

            start_inbetween = pd.concat([annotsv_df.loc[i, :].reset_index(drop=True),dbvar_df.loc[j, :].reset_index(drop=True)], axis=1)

            annotSV_end = annotsv_df.SV_end.values
            dbvar_end = dbvar_df.END.values
            dbvar_start = dbvar_df.POS.values

            i, j = np.where((annotSV_end[:, None] >= dbvar_start) & (annotSV_end[:, None] <= dbvar_end))

            stop_inbetween = pd.concat([annotsv_df.loc[i, :].reset_index(drop=True),dbvar_df.loc[j, :].reset_index(drop=True)], axis=1)

            total_overlap = pd.concat([start_inbetween,stop_inbetween]).reset_index(drop=True)

            total_overlap = total_overlap.drop_duplicates(subset = ['AnnotSV_ID', 'ID'], keep = 'last').reset_index(drop = True)

            total_overlap["distance_start_start"] = total_overlap["SV_start"]-total_overlap["POS"]
            total_overlap["distance_start_end"] = total_overlap["SV_start"]-total_overlap["END"]
            total_overlap["distance_end_end"] = total_overlap["SV_end"]-total_overlap["END"]
            total_overlap["distance_end_start"] = total_overlap["SV_end"]-total_overlap["POS"]

            list_of_overlap_tests.append(total_overlap)
        except Exception:
            st.write("merging failed, please provide properly formatted .vcf and .tsv files (chrom without chr)")

        pd.concat(list_of_overlap_tests)
        CNVExporter()
        CNVExporter.save_tables_as_excel_tsv(path_to_output_xlsx)
