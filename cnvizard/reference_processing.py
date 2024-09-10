"""
Utility functions for reference processing in CNVizard
Authors: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
Company: University Hospital Aachen
Email: jerkrause@ukaachen.de
"""

import pandas as pd
import os
import numpy as np


def prepare_cnv_table(df: pd.DataFrame, df2: pd.DataFrame) -> pd.DataFrame:
    """
    Format and reorder the .cnr DataFrame.

    Parameters:
        df (pd.DataFrame): .cnr DataFrame
        df2 (pd.DataFrame): OMIM DataFrame

    Returns:
        pd.DataFrame: Reordered and formatted DataFrame
    """
    df = df[~df["gene"].str.contains("Antitarget")]
    df["squaredvalue"] = 2 ** df["log2"]
    df[["gene", "exon"]] = df["gene"].str.split("_", expand=True)
    df["exon"] = df["exon"].astype(int)
    cols = df.columns.tolist()
    cols = cols[:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
    df = df[cols]
    df.insert(loc=7, column="call", value="")
    df.loc[df["log2"] <= -1.1, "call"] = 0
    df.loc[(df["log2"] <= -0.4) & (df["log2"] > -1.1), "call"] = 1
    df.loc[(df["log2"] <= 0.3) & (df["log2"] > -0.4), "call"] = 2
    df.loc[df["log2"] > 0.3, "call"] = 3
    df = pd.merge(df, df2, on="gene", how="left")
    df.fillna("-", inplace=True)
    df["comments"] = "."
    return df

def prepare_ref_cnv_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Format and reorder the .cnr DataFrame.

    Parameters:
        df (pd.DataFrame): .cnr DataFrame

    Returns:
        pd.DataFrame: Reordered and formatted DataFrame
    """
    df = df[~df["gene"].str.contains("Antitarget")]
    df["squaredvalue"] = 2 ** df["log2"]
    df[["gene", "exon"]] = df["gene"].str.split("_", expand=True)
    df["exon"] = df["exon"].astype(int)
    cols = df.columns.tolist()
    cols = cols[:4] + cols[-1:] + cols[5:6] + cols[6:7] + cols[4:5] + cols[7:-1]
    df = df[cols]
    df.insert(loc=7, column="call", value="")
    df.loc[df["log2"] <= -1.1, "call"] = 0
    df.loc[(df["log2"] <= -0.4) & (df["log2"] > -1.1), "call"] = 1
    df.loc[(df["log2"] <= 0.3) & (df["log2"] > -0.4), "call"] = 2
    df.loc[df["log2"] > 0.3, "call"] = 3
    return df


def explode_cnv_table(df: pd.DataFrame) -> pd.DataFrame:
    """
    Split redundant exons and explode the DataFrame.

    Parameters:
        df (pd.DataFrame): .cnr DataFrame

    Returns:
        pd.DataFrame: Formatted .cnr DataFrame
    """
    df["gene"] = df["gene"].str.split(",")
    df = df.explode("gene")
    return df


def _get_call_counts(call: list) -> list:
    """
    Calculate the amount of deletions, duplications, wild types from a call list.

    Parameters:
        call (list): List inside a pandas DataFrame column which stores all calls (per exon) inside the reference DataFrame.

    Returns:
        list: List inside a pandas DataFrame column which stores the sums of all call-types (per exon) inside the reference DataFrame.
    """
    dels_het, dels_hom, dups, wt = 0, 0, 0, 0
    for x in call:
        if x == 1:
            dels_het += 1
        elif x == 0:
            dels_hom += 1
        elif x == 2:
            wt += 1
        else:
            dups += 1
    all_calls = dels_het + dels_hom + wt + dups
    return [dels_het, dels_hom, dups, wt, all_calls]


def _get_frequency(call_counts: list, position1: int, position2: int) -> float:
    """
    Calculate the frequency of different call types inside the reference DataFrame.

    Parameters:
        call_counts (list): List which contains the previously calculated sums for each call-types.
        position1 (int): This position selects the call-type for which the frequency is to be calculated.
        position2 (int): This position selects the total sum of individuals inside the reference DataFrame.

    Returns:
        float: Frequency of a CNV-type per exon inside the reference DataFrame.
    """
    number_of_cnvs = call_counts[position1]
    number_of_individuals = call_counts[position2]
    frequency = number_of_cnvs / number_of_individuals
    return frequency


def _get_frequency_bintest(
    call_counts: list, position1: int, position2: int, call_counts_ref: list
) -> float:
    """
    Calculate the frequency of different call types inside the reference DataFrame for bintest.

    Parameters:
        call_counts (list): List which contains the previously calculated sums for each call-types.
        position1 (int): This position selects the call-type for which the frequency is to be calculated.
        position2 (int): This position selects the total sum of individuals inside the reference DataFrame.
        call_counts_ref (list): List which contains the previously calculated sums for each call-types (from the "normal" reference).

    Returns:
        float: Frequency of a CNV-type per exon inside the reference DataFrame.
    """
    number_of_cnvs = call_counts[position1]
    try:
        number_of_individuals = call_counts_ref[position2]
    except Exception:
        number_of_individuals = call_counts[position2]
    frequency = number_of_cnvs / number_of_individuals
    return frequency


def merge_reference_files(
    path_to_input: str, path_to_output: str, path_to_bintest: str
):
    """
    Merge previously created individual reference files into a single reference file for plotting with CNVizard.

    Parameters:
        path_to_input (str): Path to the input directory containing the individual reference files.
        path_to_output (str): Path to the output directory where the merged reference file will be saved.
        path_to_bintest (str): Path to the directory containing the bintest reference files.
    """
    # Load individual reference files
    reference_files = [
        os.path.join(path_to_input, f)
        for f in os.listdir(path_to_input)
        if f.endswith(".parquet")
    ]
    reference_dfs = [pd.read_parquet(file) for file in reference_files]

    # Concatenate all individual reference DataFrames
    reference_df = pd.concat(reference_dfs, ignore_index=True)

    # Drop unnecessary columns
    reference_df.drop(
        columns=[
            "chromosome",
            "start",
            "end",
            "weight",
            "squaredvalue",
            #"OMIMG",
            #"Disease",
            #"OMIMP",
            #"Inheritance",
            #"comments",
        ],
        inplace=True,
    )

    # Group by gene and exon, aggregate columns
    reference_df = (
        reference_df.groupby(["gene", "exon"])
        .agg({"depth": list, "log2": list, "call": list})
        .reset_index()
    )

    # Apply call count and frequency functions
    reference_df["call_counts"] = reference_df["call"].apply(_get_call_counts)
    reference_df["het_del_frequency"] = reference_df["call_counts"].apply(
        lambda x: _get_frequency(x, 0, 4)
    )
    reference_df["hom_del_frequency"] = reference_df["call_counts"].apply(
        lambda x: _get_frequency(x, 1, 4)
    )
    reference_df["dup_frequency"] = reference_df["call_counts"].apply(
        lambda x: _get_frequency(x, 2, 4)
    )

    # Calculate additional statistics
    reference_df["max_log2"] = reference_df["log2"].apply(max)
    reference_df["min_log2"] = reference_df["log2"].apply(min)
    reference_df["mean_depth"] = reference_df["depth"].apply(np.mean)
    reference_df["mean_log2"] = reference_df["log2"].apply(np.mean)
    reference_df["median_depth"] = reference_df["depth"].apply(
        lambda x: np.quantile(x, 0.5)
    )
    reference_df["median_log2"] = reference_df["log2"].apply(
        lambda x: np.quantile(x, 0.5)
    )
    reference_df["q1_depth"] = reference_df["depth"].apply(
        lambda x: np.quantile(x, 0.25)
    )
    reference_df["q1_log2"] = reference_df["log2"].apply(lambda x: np.quantile(x, 0.25))
    reference_df["q3_depth"] = reference_df["depth"].apply(
        lambda x: np.quantile(x, 0.75)
    )
    reference_df["q3_log2"] = reference_df["log2"].apply(lambda x: np.quantile(x, 0.75))
    reference_df["std_depth"] = reference_df["depth"].apply(np.std)
    reference_df["std_log2"] = reference_df["log2"].apply(np.std)
    reference_df["box_size"] = (reference_df["q3_log2"] - reference_df["q1_log2"]) * 1.5
    reference_df["actual_minimum_log2"] = (
        reference_df["q1_log2"] - reference_df["box_size"]
    )
    reference_df["actual_maximum_log2"] = (
        reference_df["q3_log2"] + reference_df["box_size"]
    )
    reference_df["box_size_depth"] = (
        reference_df["q3_depth"] - reference_df["q1_depth"]
    ) * 1.5
    reference_df["actual_minimum_depth"] = (
        reference_df["q1_depth"] - reference_df["box_size_depth"]
    )
    reference_df["actual_maximum_depth"] = (
        reference_df["q3_depth"] + reference_df["box_size_depth"]
    )

    # Create a new DataFrame with call counts from normal references
    ref_counts_df = reference_df[["gene", "exon", "call_counts"]]

    # Load bintest reference files
    bintest_files = [
        os.path.join(path_to_bintest, f)
        for f in os.listdir(path_to_bintest)
        if f.endswith(".parquet")
    ]
    bintest_dfs = [pd.read_parquet(file) for file in bintest_files]

    # Concatenate all bintest reference DataFrames
    bintest_df = pd.concat(bintest_dfs, ignore_index=True)
    bintest_df.drop(
        columns=[
            "chromosome",
            "start",
            "end",
            "weight",
            "squaredvalue",
            #"OMIMG",
            #"Disease",
            #"OMIMP",
            #"Inheritance",
            #"comments",
        ],
        inplace=True,
    )

    # Group by gene and exon, aggregate columns
    bintest_df = (
        bintest_df.groupby(["gene", "exon"])
        .agg({"depth": list, "log2": list, "call": list})
        .reset_index()
    )
    bintest_df["call_counts"] = bintest_df["call"].apply(_get_call_counts)

    # Merge with bintest reference DataFrame
    bintest_df = pd.merge(bintest_df, ref_counts_df, on=["gene", "exon"], how="left")
    bintest_df["het_del_frequency"] = bintest_df.apply(
        lambda x: _get_frequency_bintest(x["call_counts_x"], 0, 4, x["call_counts_y"]),
        axis=1,
    )
    bintest_df["hom_del_frequency"] = bintest_df.apply(
        lambda x: _get_frequency_bintest(x["call_counts_x"], 1, 4, x["call_counts_y"]),
        axis=1,
    )
    bintest_df["dup_frequency"] = bintest_df.apply(
        lambda x: _get_frequency_bintest(x["call_counts_x"], 2, 4, x["call_counts_y"]),
        axis=1,
    )
    bintest_df.drop(
        columns=["depth", "log2", "call", "call_counts_x", "call_counts_y"],
        inplace=True,
    )

    # Write bintest reference to parquet file
    bintest_df.to_parquet(os.path.join(path_to_output, "cnv_reference_bintest.parquet"))

    # Drop unnecessary columns from reference DataFrame
    reference_df.drop(
        columns=[
            "depth",
            "log2",
            "call",
            "call_counts",
            "max_log2",
            "min_log2",
            "box_size",
            "box_size_depth",
        ],
        inplace=True,
    )

    # Write merged and formatted reference to parquet file
    reference_df.to_parquet(os.path.join(path_to_output, "cnv_reference.parquet"))


def create_reference_files(
       path_to_input,path_to_output, reference_type
):
    """
    Create individual reference files for CNV visualization.

    Parameters:
        path_to_input (str): Path to the input directory.
        path_to_output (str): Path to the output directory.
        reference_type (str): Type of reference ('normal' or 'bintest').
    """
    # Define run name
    run_name = os.path.basename(os.path.normpath(path_to_input))

    # Load OMIM reference
    #omim_all = os.path.join(omim_path, "omim.txt")
    #omim_df = pd.read_csv(omim_all, delimiter="\t")

    list_of_dfs = []
    if reference_type == "normal":
        for cnr_file in os.listdir(path_to_input):
            if cnr_file.endswith(".cnr"):
               # st.write(path_to_input+cnr_file)
                current_total_df = pd.read_csv(path_to_input+cnr_file, delimiter="\t")
                list_of_dfs.append(current_total_df)
    else:
        for cnr_file in os.listdir(path_to_input):
            if cnr_file.endswith("_bintest.tsv"):
                current_total_df = pd.read_csv(path_to_input+cnr_file, delimiter="\t")
                list_of_dfs.append(current_total_df)

    # Concatenate all individual dataframes to a reference dataframe
    reference = pd.concat(list_of_dfs)
    exploded_reference = explode_cnv_table(reference)
    ordered_reference = prepare_ref_cnv_table(exploded_reference)

    # Define export name
    export_name_parquet = os.path.join(
        path_to_output,
        f'CNV_reference_{
            "bintest_" if reference_type != "normal" else ""}{run_name}.parquet',
    )

    # Export reference to parquet file
    ordered_reference.to_parquet(export_name_parquet, index=False)

def convert_genomics_england_panel_to_txt(
        path_to_input:str,path_to_output:str
    ):
    """
    Convert a tab delimited file from the genomics england panel app to a .txt file 
    that is compatible with the CNVizard.

    Parameters:
        path_to_input (str): Path to the input directory.
        path_to_output (str): Path to the output directory.
    """
    # read panel app list
    panel_app_list = pd.read_csv(path_to_input, sep="\t")
    # extract gene list from panel app dataframe
    gene_list = panel_app_list["Gene Symbol"]
    gene_list = gene_list.drop_duplicates()
    # write extracted gene list into newly created .txt file
    gene_list.to_csv(path_to_output, index=False, header=False)
