"""
File which is used to run the CNVizard app
@author: Jeremias Krause , Carlos Classen, Matthias Begemann. Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import os
import dotenv
import argparse
import streamlit as st
import pandas as pd
import seaborn as sns
import subprocess
import tempfile
import importlib
from cnvizard import (
    make_pretty,
    CNVExporter,
    CNVPlotter,
    filter_tsv,
    CNVVisualizer,
    merge_reference_files,
    create_reference_files,
)
from pathlib import Path
from cyvcf2 import VCF
from tempfile import NamedTemporaryFile
from reference_processing import convert_genomics_england_panel_to_txt
sns.set_theme()


def get_resources_base_path():
    """
    Get the base path to the 'resources' directory within the package using importlib.resources.

    Returns:
        Path: Path to the 'resources' directory.
    """
    try:
        # Use importlib.resources.files to access the 'resources' directory in the 'cnvizard' package
        resources_path = importlib.resources.files("cnvizard").joinpath("resources")

        # Use as_file to access the resources as an actual file path
        with importlib.resources.as_file(resources_path) as path:
            if not path.exists():
                raise FileNotFoundError(f"Resources directory not found: {path}")
            return path

    except Exception as e:
        raise FileNotFoundError(f"Error locating resources directory: {e}")



def validate_env_file(env_file):
    """
    Validates the contents of an environment file.

    Args:
        env_file (str): Path to the environment file.

    Returns:
        tuple: A tuple containing two lists:
            - missing_keys (list): List of missing keys.
            - invalid_paths (list): List of invalid paths.
    """
    required_keys = [
        "CANDIDATE_LIST_DIR",
        "REFERENCE_FILES_DIR",
        "ANNOTS_SV_FORMAT_PATH",
    ]
    optional_keys = ["APPSETTING_IGV_OUTLINK","OMIM_ANNOTATION_PATH"]
    missing_keys = []
    invalid_paths = []

    with open(env_file, "r") as file:
        lines = file.readlines()
        env_dict = dict(line.strip().split("=") for line in lines if "=" in line)

    for key in required_keys:
        if key not in env_dict:
            missing_keys.append(key)
        elif not Path(env_dict[key]).exists():
            invalid_paths.append(env_dict[key])

    for key in optional_keys:
        if key in env_dict and not Path(env_dict[key]).exists():
            invalid_paths.append(env_dict[key])

    return missing_keys, invalid_paths


def write_env_file(omim_path, candidate_list_dir, reference_files_dir, annotsv_format_path):
    """
    Writes a temporary .env file with the provided paths.
    
    Args:
        omim_path (str): Path to the OMIM annotation file.
        candidate_list_dir (str): Path to the candidate list directory.
        reference_files_dir (str): Path to the reference files directory.
        annotsv_format_path (str): Path to the Annotsv format file.
    
    Returns:
        str: Path to the created .env file.
    """
    env_content = f"""OMIM_ANNOTATION_PATH={omim_path}
CANDIDATE_LIST_DIR={candidate_list_dir}
REFERENCE_FILES_DIR={reference_files_dir}
ANNOTS_SV_FORMAT_PATH={annotsv_format_path}
"""
    # Create a temporary file
    temp_dir = tempfile.mkdtemp()
    env_file_path = os.path.join(temp_dir, ".env")

    # Write the content to the .env file
    with open(env_file_path, "w") as env_file:
        env_file.write(env_content)
    
    return env_file_path

def load_and_select_env(unique_id="1"):
    """
    Streamlit function to select between using bundled resources or uploading an environment file.
    
    This function allows users to choose between using bundled resources, uploading an environment file,
    or uploading an additional OMIM annotation file.

    Returns:
        str: Path to the loaded or created environment file, or None if there was an error.
    """
    st.markdown("### Select Environment Setup")
    
    # Get the base path to the resources directory
    resources_base_path = get_resources_base_path()

    # Ensure unique key for each radio widget
    env_choice = st.radio(
        "Select environment configuration:",
        ("Use Default Bundled Resources", "Upload Custom Environment File"),
        key=f"env_choice_{unique_id}"  # Add a unique ID to the key
    )

    # Option 1: Use Default Bundled Resources but require OMIM file upload
    if env_choice == "Use Default Bundled Resources":
        st.markdown("#### Upload OMIM Annotation File")
        omim_file = st.file_uploader("Upload OMIM annotation file", type=["txt"], key=f"omim_file_{unique_id}")
        if omim_file:
            st.success(f"Uploaded OMIM file: {omim_file.name}")
            # When OMIM file is uploaded, show a confirmation button
            try:
                # Save the uploaded OMIM file locally
                omim_file_path = omim_file.name
                with open(omim_file_path, "wb") as file:
                    file.write(omim_file.getbuffer())
            except Exception as e:
                st.error(f"Error loading OMIM file: {e}")
                return None  # Return None in case of error
            
        if st.button("Use default Resources", key=f"load_resources_{unique_id}"):
            if omim_file:
                env_file_path = write_env_file(
                    omim_file_path,
                    str(resources_base_path / "candidate_lists"),
                    str(resources_base_path / "references"),
                    str(resources_base_path / "annotsv_format.txt")
                )
            else:
                env_file_path = write_env_file(
                    "",
                    str(resources_base_path / "candidate_lists"),
                    str(resources_base_path / "references"),
                    str(resources_base_path / "annotsv_format.txt")
                )

            st.success("Bundled resources loaded successfully.")
            return env_file_path 

    # Option 2: Upload custom .env file and use the loader and creator as it is
    elif env_choice == "Upload Custom Environment File":
        uploaded_env_file = st.file_uploader("Upload environment file")

        if uploaded_env_file:
            try:
                env_file = uploaded_env_file.name
                with open(env_file, "wb") as file:
                    file.write(uploaded_env_file.getbuffer())

                missing_keys, invalid_paths = validate_env_file(env_file)
                if missing_keys or invalid_paths:
                    error_message = "Error in environment file:\n"
                    if missing_keys:
                        error_message += f"Missing keys: {', '.join(missing_keys)}\n"
                    if invalid_paths:
                        error_message += f"Invalid paths: {', '.join(invalid_paths)}"
                    raise ValueError(error_message)

                dotenv.load_dotenv(env_file)
                st.success(f"Loaded environment file: {env_file}")
                return env_file  # Return the loaded env file path
            except Exception as e:
                st.error(f"Error loading environment file: {e}")
                return None

    return None


def create_new_reference_file():
    """
    Streamlit function to create a new reference file.

    This functions uses the create_reference_files function .
    """
    st.title("CNVizard Reference creator")
    st.markdown(
        """
    ### Create a new Reference file by passing an input path, an output path and a reference type (normal or bintest)
    """
    )
    path_to_individual_cnr_files = st.text_input(
            "path to cnr files:" ,"/path/to/cnrs/"
        )
    path_to_output_for_created_reference = st.text_input(
            "output path for created reference:" ,"/path/to/output/"
        )
    reference_type = st.multiselect(
            "reference type", ["normal","bintest"],max_selections=1
        )
    if st.button("Create new reference files"):
        create_reference_files(path_to_individual_cnr_files,path_to_output_for_created_reference,reference_type)

def merge_new_reference_file():
    """
    Streamlit function to merge a new reference file.

    This functions uses the merge_reference_files function.
    """
    st.title("CNVizard Reference merger")
    st.markdown(
        """
    ### Merge and format previously created Reference files
    """
    )
    path_to_cnr_references = st.text_input(
            "path to cnr references:" ,"/path/to/cnrreferences/"
        )
    path_to_output_for_merged_reference = st.text_input(
            "output path for merged and reformatted references:" ,"/path/to/output/"
        )
    path_to_bintest_references = st.text_input(
            "path to bintest references:" ,"/path/to/bintestreferences/"
        )
    if st.button("Merge created reference files"):
        merge_reference_files(path_to_cnr_references,path_to_output_for_merged_reference,path_to_bintest_references)

def convert_genomics_england_panel():
    """
    Streamlit function to make genomics england panels compatible with CNVizard.

    This functions uses the convert_genomics_england_panel_to_txt function.
    """
    st.title("Panel Converter")
    st.markdown(
        """
    ### Make genomics england panels compatible with CNVizard
    """
    )
    path_to_txt_panel = st.text_input(
            "path to .tsv panel file:" ,"/path/to/cnrreferences/genomics_panel.tsv"
        )
    path_to_output_for_panel = st.text_input(
            "output path for panel file:" ,"/path/to/output/new_gene_panel.txt"
        )
    if st.button("convert panel list"):
        convert_genomics_england_panel_to_txt(path_to_txt_panel,path_to_output_for_panel)


def prepare_filter_for_consecutive_cnvs(self, df: pd.DataFrame) -> pd.DataFrame:
    """
    Function which is used to filter for consecutively deleted/duplicated exons.

    Args:
        df (pd.DataFrame): .cnr DataFrame

    Returns:
        pd.DataFrame: .cnr DataFrame filtered for consecutively deleted/duplicated exons.
    """
    df = df.copy()  # Create a copy to avoid SettingWithCopyWarning
    df["difference_previous"] = df.groupby("gene")["exon"].diff()
    df["difference_previous"] = df["difference_previous"].fillna(method="bfill")
    df["difference_next"] = df.groupby("gene")["exon"].diff(periods=-1)
    df["difference_next"] = df["difference_next"].fillna(method="ffill")
    return df


def main(env_file_path):
    # Streamlit App Title and icon
    st.set_page_config(
    page_title="CNVizard - Copy Number Variant Visualization Tool", 
    page_icon="./CNVizard.png",
    layout="wide"
    )
    st.markdown(
        "This is a Streamlit web app providing analysis tools for genetic copy number variants. Please first make an environment selection by pressing the button 'Use default Resources' or upload a custom environment file. Using an OMIM file is optional in both cases."
    )

    if env_file_path is None:
        env_file_path = load_and_select_env(unique_id="main")
    
    # Load environment variables
    dotenv.load_dotenv(env_file_path)
    igv_string = os.getenv("APPSETTING_IGV_OUTLINK")

    # Load paths from environment variables
    omim_annotation_path = Path(
        os.getenv("OMIM_ANNOTATION_PATH", "")
    )
    resources_base_path = get_resources_base_path() # For default candidate lists
    candidate_list_dir = Path(
        os.getenv("CANDIDATE_LIST_DIR", str(resources_base_path / "candidate_lists"))
    )
    reference_files_dir = Path(
        os.getenv("REFERENCE_FILES_DIR", "")
    )
    annotsv_format_path = Path(
        os.getenv("ANNOTS_SV_FORMAT_PATH", "")
    )

    # Filter options
    chrom_list = [f"chr{i}" for i in range(1, 23)] + ["chrX", "chrY"]
    call_list = [0, 1, 2, 3]

    # File Uploads Section
    st.subheader("Upload")
    st.markdown(
        "Please upload a reference file, an individual .cnr file, and an individual .bintest file, provided by CNVkit."
    )

    omim_annotation_file = None
    gene_list = []
    omim_annotation_file = omim_annotation_path

    #if omim_annotation_path.exists():
    #    omim_annotation_file = omim_annotation_path
    #    gene_list = pd.read_csv(omim_annotation_file, sep="\t")["gene"].tolist()

    cols = st.columns(3)
    entered_cnr = cols[0].file_uploader(".cnr", type=["txt", "cnr"])
    entered_bintest = cols[1].file_uploader("bintest", type=["txt", "tsv"])
    ngs_type = cols[2].radio("Choose NGS Type", ["WES", "WGS"], key="ngs_type")

    reference_path = reference_files_dir / (
        "genome_cnv_reference_large.parquet"
        if ngs_type == "WGS"
        else "exome_cnv_reference_large.parquet"
    )
    reference_bintest_path = reference_files_dir / (
        "genome_cnv_reference_bintest_large.parquet"
        if ngs_type == "WGS"
        else "exome_cnv_reference_bintest_large.parquet"
    )

    sample_name = ""
    reference_df = None
    cnr_df = None
    bintest_df = None
    if entered_cnr:
        sample_name = entered_cnr.name.split(".")[0]
        try:
            cnr_df = pd.read_csv(entered_cnr, delimiter="\t")
            gene_list = list(set(cnr_df["gene"].str.split("_").str[0].tolist()))
        except Exception as e:
            st.error(f"Error reading .cnr file: {e}")
            cnr_df = None
        if igv_string:
            igv_string = igv_string.replace("samplename", sample_name)
    if entered_bintest:
        try:
            bintest_df = pd.read_csv(entered_bintest, delimiter="\t")
        except Exception as e:
            st.error(f"Error reading bintest file: {e}")
            bintest_df = None
    if reference_path.exists():
        try:
            reference_df = pd.read_parquet(reference_path)
        except Exception as e:
            st.error(f"Error reading reference file: {e}")
            reference_df = None
    try:
        reference_bintest_df = (
            pd.read_parquet(reference_bintest_path)
            if reference_bintest_path.exists()
            else pd.DataFrame()
        )
    except Exception as e:
        st.error(f"Error reading reference bintest file: {e}")
        reference_bintest_df = pd.DataFrame()

    st.subheader("Configurations")
    st.markdown(
        "Provide a sample name and define the number of consecutive exons to display (default value = 2)."
    )

    cols2 = st.columns(2)
    entered_del_size = cols2[0].text_input("deletion_size", value="2")
    entered_dup_size = cols2[1].text_input("duplication_size", value="2")

    # Sidebar configuration
    st.sidebar.title("About")
    st.sidebar.markdown(
        "This Streamlit web app enables you to visualize copy-number-variant data using dataframes and plots."
    )

    st.sidebar.subheader("Candidate Gene Selection")
    candidate_options = (
        sorted(os.listdir(candidate_list_dir)) if candidate_list_dir.exists() else []
    )
    selected_candidate = st.sidebar.radio(
        "Select the desired candidate gene list", candidate_options, key="candidate"
    )
    st.sidebar.write(selected_candidate)
    st.sidebar.markdown("Create a new env file, reference file or a new panel list")
    # Expander to edit env_file
    with st.sidebar.expander("change_env_file"):
        load_and_select_env(unique_id="sidebar")
    
    #Expander for reference creation
    with st.sidebar.expander("create reference files"):
        create_new_reference_file()
    
    #Expander for reference merging
    with st.sidebar.expander("merge reference files"):
        merge_new_reference_file()
    
    #Expander for genomics england panel conversion
    with st.sidebar.expander("convert panel files"):
        convert_genomics_england_panel()

    st.subheader("Dataframe Visualization")
    with st.expander("Filter"):
        cols3 = st.columns(3)
        chrom_selection = cols3[0].multiselect("Chromosome", chrom_list)
        start_selection = cols3[1].text_input("start", value="")
        end_selection = cols3[2].text_input("end", value="")

        cols4 = st.columns(3)
        depth_selection = cols4[0].text_input("min_depth", value="")
        weight_selection = cols4[1].text_input("min_weight", value="")
        call_selection = cols4[2].multiselect("call", call_list)

        cols5 = st.columns(3)
        log2_selection = cols5[0].text_input("min_log2", value="")
        gene_selection = cols5[1].multiselect("gene", gene_list)
        het_del_selection = cols5[2].text_input("max_het_del_freq", value="")

        cols7 = st.columns(3)
        hom_del_selection = cols7[0].text_input("max_hom_del_freq", value="")
        dup_selection = cols7[1].text_input("max_dup_freq", value="")

    list_of_possible_dataframes = [
        "total",
        "bintest",
        "hom_del",
        "total_candidate",
        "bintest_candidate",
        "consecutive_del",
        "consecutive_dup",
    ]
    df_to_be_displayed = st.selectbox(
        "Select dataframe to display", list_of_possible_dataframes
    )

    if reference_df is not None and cnr_df is not None and bintest_df is not None:
        cnv_visualizer_instance = CNVVisualizer(reference_df, cnr_df, bintest_df)
        omim_df, candidate_df, cnr_db, bintest_db = cnv_visualizer_instance.format_df(
            omim_annotation_file,
            (
                os.path.join(candidate_list_dir, selected_candidate)
                if candidate_list_dir.exists()
                else None
            ),
        )

        call_df = reference_df[
            ["gene", "exon", "het_del_frequency", "hom_del_frequency", "dup_frequency"]
        ]
        bintest_inhouse_df = (
            reference_bintest_df if not reference_bintest_df.empty else pd.DataFrame()
        )

        cnr_db = pd.merge(cnr_db, call_df, on=["gene", "exon"], how="left")
        if not bintest_inhouse_df.empty:
            bintest_db = pd.merge(
                bintest_db, bintest_inhouse_df, on=["gene", "exon"], how="left"
            )
            bintest_db = bintest_db.infer_objects().fillna(0)

        cnr_db_filtered = cnv_visualizer_instance.apply_filters(
            cnr_db,
            start_selection,
            end_selection,
            depth_selection,
            weight_selection,
            chrom_selection,
            call_selection,
            log2_selection,
            gene_selection,
            chrom_list,
            call_list,
            gene_list,
            het_del_selection,
            hom_del_selection,
            dup_selection,
        )

        cnr_db_filtered.rename(columns={"chromosome": "chr"}, inplace=True)
        cnr_db_filtered = cnr_db_filtered.round(2)
        bintest_db.rename(columns={"chromosome": "chr"}, inplace=True)
        bintest_db = bintest_db.round(2)

        if igv_string:
            cnr_db_filtered["IGV_outlink"] = (
                igv_string
                + cnr_db_filtered["chr"]
                + ":"
                + cnr_db_filtered["start"].astype(str)
            )
            bintest_db["IGV_outlink"] = (
                igv_string + bintest_db["chr"] + ":" + bintest_db["start"].astype(str)
            )

        # Dataframe display logic
        display_mapping = {
            "total": cnr_db_filtered,
            "bintest": bintest_db,
            "hom_del": cnv_visualizer_instance.filter_for_deletions_hom(
                cnr_db_filtered
            ),
            "total_candidate": cnv_visualizer_instance.filter_for_candi_cnvs(
                cnr_db_filtered, candidate_df
            ),
            "bintest_candidate": cnv_visualizer_instance.filter_for_candi_cnvs(
                bintest_db, candidate_df
            ),
            "consecutive_del": cnv_visualizer_instance.filter_for_consecutive_cnvs(
                cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.filter_for_deletions(cnr_db_filtered)
                ),
                "del",
                entered_del_size,
                entered_dup_size,
            ),
            "consecutive_dup": cnv_visualizer_instance.filter_for_consecutive_cnvs(
                cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.filter_for_duplications(cnr_db_filtered)
                ),
                "dup",
                entered_del_size,
                entered_dup_size,
            ),
        }

        download_filter = display_mapping[df_to_be_displayed]
        st.dataframe(
            display_mapping[df_to_be_displayed].style.pipe(make_pretty)
            if df_to_be_displayed != "total"
            else display_mapping[df_to_be_displayed]
        )

        # Download buttons
        download_columns = st.columns(2)
        download_preparator_all = download_columns[0].button(
            "Prepare for download (all)"
        )
        download_preparator_filtered = download_columns[1].button(
            "Prepare for download (filtered)"
        )

        download_message_columns = st.columns(2)
        download_button_columns = st.columns(2)

        if download_preparator_all:
            tables_to_export = [
                cnr_db,
                bintest_db,
                cnv_visualizer_instance.filter_for_deletions_hom(cnr_db),
                cnv_visualizer_instance.filter_for_candi_cnvs(cnr_db, candidate_df),
                cnv_visualizer_instance.filter_for_candi_cnvs(bintest_db, candidate_df),
                cnv_visualizer_instance.filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                        cnv_visualizer_instance.filter_for_deletions(cnr_db)
                    ),
                    "del",
                    entered_del_size,
                    entered_dup_size,
                ),
                cnv_visualizer_instance.filter_for_consecutive_cnvs(
                    cnv_visualizer_instance.prepare_filter_for_consecutive_cnvs(
                        cnv_visualizer_instance.filter_for_duplications(cnr_db)
                    ),
                    "dup",
                    entered_del_size,
                    entered_dup_size,
                ),
            ]
            table_exporter = CNVExporter()
            export_data = table_exporter.save_tables_as_excel(*tables_to_export)
            download_button_columns[0].download_button(
                label="Download", data=export_data, file_name=f"{sample_name}_df.xlsx"
            )
        else:
            download_message_columns[0].write("Click button to prepare download")

        if download_preparator_filtered:
            table_exporter = CNVExporter()
            export_data_filtered = table_exporter.save_filtered_table_as_excel(
                download_filter, df_to_be_displayed
            )
            download_button_columns[1].download_button(
                label="Download filtered",
                data=export_data_filtered,
                file_name=f"{sample_name}_filtered_df.xlsx",
            )
        else:
            download_message_columns[1].write("Click button to prepare download")

    st.subheader("Plot log2 and depth for selected genes")
    entered_gene = st.multiselect("gene", gene_list, max_selections=1)
    entered_gene = entered_gene[0].upper() if entered_gene else None

    if (
        reference_df is not None
        and cnr_df is not None
        and bintest_df is not None
        and entered_gene
    ):
        gene_plotter = CNVPlotter()
        gene_plotter.plot_log2_for_gene_precomputed(
            entered_gene, cnr_db, reference_df, sample_name
        )
        gene_plotter.plot_depth_for_gene_precomputed(
            entered_gene, cnr_db, reference_df, sample_name
        )

    st.subheader("Load additional .cnr files from the index patient's parents")
    with st.expander("Filter for Trio"):
        cols_for_trio_filter = st.columns(3)
        call_selection_index = cols_for_trio_filter[0].multiselect(
            "call index", call_list
        )
        call_selection_father = cols_for_trio_filter[1].multiselect(
            "call father", call_list
        )
        call_selection_mother = cols_for_trio_filter[2].multiselect(
            "call mother", call_list
        )

    cols_trio = st.columns(3)
    father_cnr = cols_trio[0].file_uploader("Father .cnr", type=["txt", "cnr"])
    mother_cnr = cols_trio[1].file_uploader("Mother .cnr", type=["txt", "cnr"])

    if (
        reference_df is not None
        and cnr_df is not None
        and bintest_df is not None
        and father_cnr
        and mother_cnr
    ):
        try:
            father_cnr_df = cnv_visualizer_instance.prepare_parent_cnv(
                pd.read_csv(father_cnr, delimiter="\t")
            )
            mother_cnr_df = cnv_visualizer_instance.prepare_parent_cnv(
                pd.read_csv(mother_cnr, delimiter="\t")
            )
        except Exception as e:
            st.error(f"Error reading parent .cnr files: {e}")
            st.stop()

        father_cnr_df = father_cnr_df.drop(["chromosome", "start", "end"], axis=1)
        father_cnr_df = father_cnr_df.rename(
            columns={
                "gene": "gene_f",
                "exon": "exon_f",
                "depth": "depth_f",
                "weight": "weight_f",
                "call": "call_f",
                "log2": "log2_f",
                "squaredvalue": "squaredvalue_f",
            }
        )
        mother_cnr_df = mother_cnr_df.drop(["chromosome", "start", "end"], axis=1)
        mother_cnr_df = mother_cnr_df.rename(
            columns={
                "gene": "gene_m",
                "exon": "exon_m",
                "depth": "depth_m",
                "weight": "weight_m",
                "call": "call_m",
                "log2": "log2_m",
                "squaredvalue": "squaredvalue_m",
            }
        )

        trio_cnr_df = cnr_db.merge(
            father_cnr_df,
            left_on=["gene", "exon"],
            right_on=["gene_f", "exon_f"],
            how="left",
        ).merge(
            mother_cnr_df,
            left_on=["gene", "exon"],
            right_on=["gene_m", "exon_m"],
            how="left",
        )

        trio_cnr_df.rename(columns={"chromosome": "chr"}, inplace=True)
        trio_cnr_df = trio_cnr_df.round(2)
        trio_cnr_df_filtered = cnv_visualizer_instance.apply_trio_filters(
            trio_cnr_df,
            call_selection_index,
            call_selection_father,
            call_selection_mother,
            call_list,
        )
        st.dataframe(trio_cnr_df_filtered)

    st.subheader("Plot genome-wide or chromosome-wide scatter plot")
    cols_cns_upload = st.columns(3)
    entered_cns = cols_cns_upload[0].file_uploader(".cns file", type=["txt", "cns"])
    entered_cnr_new = cols_cns_upload[1].file_uploader(".cnr file", type=["txt", "cnr"])
    entered_vcf_new = cols_cns_upload[2].file_uploader(".vcf file", type=["vcf.gz"])

    scatter_options = chrom_list + ["All"]
    with st.expander("Select chromosome"):
        selected_scatter = st.radio("Select chromosome or all", scatter_options, key="scatter")

    if entered_cns and entered_cnr_new and not entered_vcf_new:
        test_cnr = pd.read_csv(entered_cnr_new,sep="\t")
        test_cns = pd.read_csv(entered_cns,sep="\t")
        scatter_plotter = CNVPlotter()
        scatter_plotter.plot_scatter(test_cnr,test_cns,selected_scatter)
    elif entered_cns and entered_cnr_new and entered_vcf_new:
        st.write(entered_vcf_new.name)
        #Import vcf.gz and write it to temp dir to import it with cyvcf2
        #This step requires tabix
        try:
            with NamedTemporaryFile("wb",suffix=".vcf.gz",delete=False) as temp_file:
                temp_file.write(entered_vcf_new.getvalue())
            st.write(temp_file.name)
            subprocess.run(["tabix", temp_file.name])
            vcf = VCF(temp_file.name)
        except Exception as e:
            st.write(e)
            st.write("VCF processing failed. Is Tabix installed ?")

        st.write(entered_cns)

        scatter_plotter = CNVPlotter()

        test_cnr = pd.read_csv(entered_cnr_new,sep="\t")
        test_cns = pd.read_csv(entered_cns,sep="\t") 

        scatter_plotter.plot_scatter(test_cnr,test_cns,selected_scatter)
        scatter_plotter.plot_baf(test_cnr,test_cns,vcf,selected_scatter)

    chromosome_list_cnv = [str(i) for i in range(1, 23)] + ["X", "Y"]
    cnv_type = ["DEL", "DUP","TRA","INS","INV"]
    acmg_class = ["1", "2", "3", "4", "5","NA"]
    st.subheader("Load annotated .tsv file created by AnnotSV")
    entered_tsv_file = st.file_uploader(".tsv file", type=["tsv"])

    with st.expander("Filter for .tsv file"):
        cols6 = st.columns(3)
        entered_cnv_chrom = cols6[0].multiselect("Chromosome", chromosome_list_cnv)
        entered_cnv_type = cols6[1].multiselect("CNV_Type", cnv_type)
        entered_acmg_class = cols6[2].multiselect("ACMG_Class", acmg_class)

    if entered_tsv_file:
        try:
            tsv_df = pd.read_csv(entered_tsv_file, delimiter="\t")
        except Exception as e:
            st.error(f"Error reading .tsv file: {e}")
            st.stop()

        with open(annotsv_format_path, "r") as column_file:
            columns_to_keep = [line.strip() for line in column_file]

        tsv_df = tsv_df[columns_to_keep]

        filtered_tsv = filter_tsv(
            tsv_df,
            chromosome_list_cnv,
            cnv_type,
            acmg_class,
            entered_cnv_chrom,
            entered_cnv_type,
            entered_acmg_class,
        )
        filtered_tsv["SV_chrom"] = "chr" + filtered_tsv["SV_chrom"]
        filtered_tsv["SV_chrom"] = pd.Categorical(filtered_tsv["SV_chrom"], chrom_list)
        filtered_tsv = filtered_tsv.sort_values("SV_chrom")

        st.write("Filtered AnnotSV DataFrame:")
        st.write(filtered_tsv)

        if st.button("Prepare for download of annotated tsv data"):
            table_exporter = CNVExporter()
            to_be_exported = table_exporter.save_tables_as_excel_tsv(filtered_tsv)
            st.download_button(
                label="Download",
                data=to_be_exported,
                file_name=f"{sample_name}_annotated_df.xlsx",
            )
        else:
            st.write("Click button to prepare download")




if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("env", nargs="?", default=None)
    args = parser.parse_args()
    main(args.env)
