"""
Init file for the CNVizard package.
"""
# Import and expose functions and classes from submodules
from .reference_processing import (
    prepare_cnv_table,
    explode_cnv_table,
    merge_reference_files,
    create_reference_files,
)
from .styler import make_pretty
from .exporter import CNVExporter
from .plotter import CNVPlotter
from .visualizer import CNVVisualizer
from .vcf_Merger import vcfMerger
from .helpers import ( 
    filter_tsv,
    merge_annotsv_with_dbvar_study,
)
from .omim import (
    phenotypes_to_gene,
    gene_to_phenotypes,
    read_genemap2,
    read_mim2gene,
    read_phenotype_to_gene,
    read_gene_to_phenotype,
    read_mimTitles
)

# Expose functions and classes
__all__ = [
    "prepare_cnv_table",
    "explode_cnv_table",
    "merge_reference_files",
    "create_reference_files",
    "make_pretty",
    "CNVExporter",
    "CNVPlotter",
    "filter_tsv",
    "merge_annotsv_with_dbvar_study",
    "CNVVisualizer",
    "vcfMerger",
    "phenotypes_to_gene",
    "gene_to_phenotypes",
    "read_genemap2",
    "read_mim2gene",
    "read_phenotype_to_gene",
    "read_gene_to_phenotype",
    "read_mimTitles"
]
