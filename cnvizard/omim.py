"""
File which contains helper functions to ensure compatibility with OMIM files and CNVizard
@author: Jeremias Krause, Carlos Classen, Daniela Dey
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd
import re 

def phenotypes_to_gene(omim_files):
    """
    Function which is used to create a phenotype file, which lists individual phenotypes.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files
    """
    df = read_genemap2(omim_files)
    df['split_p'] = df['Phenotypes'].str.split('; ')
    data = []
    for index, row in df.iterrows():
        if isinstance(row['split_p'], list):
            for item in row['split_p']:
                split_phenotypes = re.split(' \(\d\)', item)

                try:
                    int(split_phenotypes[0][-6:])
                    omim_p = split_phenotypes[0][-6:]
                except ValueError:
                    omim_p = ''

                if omim_p:
                    phenotypes = split_phenotypes[0][:-8]
                else:
                    phenotypes = split_phenotypes[0]

                inheritance = split_phenotypes[1].lstrip(', ')
                data.append([phenotypes, omim_p, row['Approved Symbol'], row['MIM Number'], inheritance])

    new_df = pd.DataFrame(data, columns=['Phenotype', 'OMIMP', 'gene', 'OMIMG', 'Inheritance'])
    new_df = new_df.replace({'Inheritance': r'Autosomal dominant?'}, {'Inheritance': 'AD'}, regex=True)
    new_df = new_df.replace({'Inheritance': r'Autosomal recessive?'}, {'Inheritance': 'AR'}, regex=True)
    new_df = new_df.replace({'Inheritance': r'Pseudoautosomal dominant?'}, {'Inheritance': 'Pseudo AD'}, regex=True)
    new_df = new_df.replace({'Inheritance': r'Pseudoautosomal recessive?'}, {'Inheritance': 'Pseudo AR'}, regex=True)
    new_df = new_df.replace({'Inheritance': r'X-linked dominant?'}, {'Inheritance': 'XLD'}, regex=True)
    new_df = new_df.replace({'Inheritance': r'X-linked recessive?'}, {'Inheritance': 'XLR'}, regex=True)
    new_df.to_csv(omim_files['phenotypes_file'], sep='\t', index=False)

def gene_to_phenotypes(omim_files):
    """
    Function which is used to create a gene file, which lists individual genes.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files
    """
    phenotypes = read_phenotype_to_gene(omim_files)
    genes = list(set(phenotypes['gene'].to_list()))
    genes = [x for x in genes if str(x) != 'nan']
    phenotypes['OMIMP'] = phenotypes['OMIMP'].astype('Int64')

    new_list = []
    for i in genes:
        phenotypes_list = []
        omimp_list = []
        inheritance_list = []
        sub_df = phenotypes[phenotypes['gene'] == i]
        for index, row in sub_df.iterrows():
            phenotypes_list.append(row['Phenotype'])
            if str(row['OMIMP']) == 'nan':
                pass
            else:
                omimp_list.append(str(row['OMIMP']))

            if str(row['Inheritance']) == 'nan':
                pass
            else:
                inheritance_list.append(str(row['Inheritance']))
            omimg = row['OMIMG']
        inheritance_list = list(set(inheritance_list))
        new_list.append([i, omimg, '; '.join(phenotypes_list), '; '.join(omimp_list), '; '.join(inheritance_list)])
    new_df = pd.DataFrame(new_list, columns=['gene', 'OMIMG', 'Disease', 'OMIMP', 'Inheritance'])
    new_df.to_csv(omim_files['genes_file'], sep='\t', index=False)

def read_genemap2(omim_files):
    """
    Function which is used to read the genemap omim file.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files

    Returns:
        pd.DataFrame: Dataframe containing the genemap omim file information
    """
    df = pd.read_csv(omim_files['genemap2_file'], comment='#', sep='\t',
                     names=['Chromosome', 'Genomic Position Start', 'Genomic Position End', 'Cyto Location',
                            'Computed Cyto Location', 'MIM Number', 'Gene Symbols', 'Gene Name', 'Approved Symbol',
                            'Entrez Gene ID', 'Ensembl Gene ID', 'Comments', 'Phenotypes', 'Mouse Gene Symbol/ID'])
    return df


def read_mim2gene(omim_files):
    """
    Function which is used to read the mim2gene omim file.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files

    Returns:
        pd.DataFrame: Dataframe containing the mim2gene omim file information
    """
    df = pd.read_csv(omim_files['mim2gene_file'], comment='#',
                     names=['MIM Number', 'type', 'entrez', 'gene', 'ensembl'], sep='\t')
    return df


def read_phenotype_to_gene(omim_files):
    """
    Function which is used to read the previously created phenotype omim file.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files

    Returns:
        pd.DataFrame: Dataframe containing the phenotype omim file information
    """
    df = pd.read_csv(omim_files['phenotypes_file'], sep='\t')
    return df


def read_gene_to_phenotype(omim_files):
    """
    Function which is used to read the previously created gene omim file.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files

    Returns:
        pd.DataFrame: Dataframe containing the gene omim file information
    """
    df = pd.read_csv(omim_files['genes_file'], sep='\t')
    return df


def read_mimTitles(omim_files):
    """
    Function which is used to read the mimTitles omim file.

    Args:
        omim_files (dictionary): dictionary containing paths to the omim files

    Returns:
        pd.DataFrame: Dataframe containing the mimTitles omim file information
    """
    df = pd.read_csv(omim_files['mimTitles_file'], comment='#',
                     names=['Prefix', 'mim_number', 'Title', 'alt_title', 'included_titles'], sep='\t')
    return df

