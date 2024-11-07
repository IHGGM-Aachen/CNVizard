# Changelog

All notable changes to this project will be documented in this file.

## [0.3] - 04/11/2024

### Initial release of CNVizard

- Streamlit-based interactive interface for visualization and analysis of germline CNVs.
- Interactive plots using Plotly for detailed CNV data exploration.
- Support for uploading and processing genetic data files, including `.cnr` and `.bintest` files.
- Filtering options for chromosome, CNV type (DUP, DEL), and ACMG classification.
- Functions for preparing and merging reference files.
- Export capabilities for processed data in Excel and Parquet formats.
- Flexible environment setup with optional `.env` file for configuration.
- IGV outlink configuration support for genome browser integration.
- Basic cross-analysis functionality utilizing dbVar.

## [0.3.1] - 07/11/2024

### Hotfix for vcf_Merger

- Fixed a `TypeError` in within `vcfMerger` class