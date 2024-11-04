"""
File which contains the vcf cnr merger class used for the CNVizard
@author: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""
import pandas as pd 
import csv
import gzip

class vcfMerger:
    """
    Class used to merge .cnr with vcf files.
    """

    def __init__(self):
        """
        Constructor of the class CNVExporter.
        """

    def read_vcf_header(self, vcf_file_path: str) -> list:
        """
        Function which extracts the header from a zipped vcf file and returns it .

        Args:
            vcf_file_path : path (str) to .vcf.gz file.

        Returns:
            lines: list of tab seperated lines.
        """
        with gzip.open (vcf_file_path, 'rt') as vcf_file: 
            lines = []
            for line in vcf_file: 
                if line.startswith('#'):
                    lines.append(line)
                else:
                    break
        return lines
    
    def read_vcf(self, vcf_file_path: str):
        """
        Function which reads a zipped vcf file as a pandas dataframe.

        Args:
            vcf_file_path : path (str) to .vcf.gz file.

        Returns:
            lines: list of tab seperated lines.
        """
        with gzip.open (vcf_file_path, 'rt') as vcf_file: 
            lines = []
            for line in vcf_file: 
                if line.startswith('#'):
                    lines.append(line)
                else:
                    break
        last_line = lines[-1]
        header = last_line.strip().split('\t')
        vcf_df = pd.read_csv(vcf_file_path, sep ="\t", comment='#', names = header,compression="gzip")
        return vcf_df
    
    def explode_df(df:pd.DataFrame) -> pd.DataFrame:
        """
        Function which splits the gene column of a pandas Dataframe on a comme and subsequently explodes the column
        Input Arguments:
            df (pandas DataFrame) : DataFrame to be exploded
        Output Arguments: 
            df (pandas DataFrame) : exploded DataFrame
        """
        df['gene'] = df['gene'].str.split(',')
        df = df.explode('gene')
        df['gene'] = df['gene'].str.split('_')
        df['exon'] = df['gene'].str[1]
        df['gene'] = df['gene'].str[0]
        
        df.loc[df['log2'] <= -1.1, 'call'] = int(0)
        df.loc[(df['log2'] <= -0.4) & (df['log2'] > -1.1), 'call'] = int(1)
        df.loc[(df['log2'] <= 0.3) & (df['log2'] > -0.4), 'call'] = int(2)
        df.loc[df['log2'] > 0.3, 'call'] = int(3)
        df.loc[df['log2'] > 1.8, 'call'] = int(4)
        
        df.loc[:, 'squaredvalue'] = 2**df['log2']
        
        return df 
    
    def prepare_filter_for_consecutive_cnvs(self,df:pd.DataFrame) -> pd.DataFrame :
        """
        Function which is used to filter for consecutively deleted/duplicated exons
        Input Arguments:
            df (pandas DataFrame) : .cnr DataFrame
        Ouput Arguments: 
            df (pandas DataFrame) : .cnr DataFrame filtered for consecutively deleted/duplicated exons
        """
        #Group CNVs by gene and calculate the difference between the exons in each group (consecutive exons will have 1/-1)
        # Because exons might be listed as 1,2,3 or 3,2,1 the calculation is performed both ways
        df['exon'] = df['exon'].astype(int)
        df['difference_previous'] = df.groupby('gene')['exon'].diff()
        df['difference_previous'] = df.groupby('gene')['difference_previous'].fillna(method='backfill')
        df['difference_next'] = df.groupby('gene')['exon'].diff(periods=-1)
        df['difference_next'] = df.groupby('gene')['difference_next'].fillna(method='ffill')
        return df
    
    def filter_for_consecutive_cnvs(self,df:pd.DataFrame,del_or_dup:str,del_size:int,dup_size:int) -> pd.DataFrame :
        """
        Function which extends the prepare_filter_for_consecutive_cnvs function
        It takes the output from the aforementioned function and selects those cnvs, which have a difference to previous/next of 1/-1
        Input Arguments:
            df (pandas DataFrame) : .cnr DataFrame annotated with consecutive cnvs 
            del_or_up (str) : String that determines wether to filter for deletions or duplications 
            del_size (int) : Int that determines the necessary amount of consecutive deletions
            dup_size (int) : Int that determines the necessary amount of consecutive duplications
        Output Arguments:
            df_cons (pandas DataFrame) : .cnr DataFrame with applied filter for consecutive deletions/duplications

        """
        #Set standard value : 
        if del_size is None or del_size=="":
            del_size = 2
        if dup_size is None or dup_size=="":
            dup_size = 2
        #Filter for consecutive deletions
        del_size = int(del_size)
        dup_size = int(dup_size)
        if del_or_dup == 'del':
            df_cons = df[(((df['call'] == 0) | (df['call'] ==1)) & ((df['difference_previous'] == -1) | (df['difference_previous'] ==1))) | (((df['call'] == 0) | (df['call'] ==1)) & ((df['difference_next'] == -1) | (df['difference_next'] ==1)))]
            #Calculate the size of consecutive deletions per gene 
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
            #Select only those genes which contain the defined amount of deleted exons or which are smaller than the defined number
            df_cons = df_cons[(df_cons['counts']>=del_size)]
        #Filter for consecutive duplications
        else:
            df_cons = df[(((df['call'] == 3) | (df['call'] == 4)) & ((df['difference_previous'] == -1) | (df['difference_previous'] ==1))) | (((df['call'] == 3)|(df['call'] == 4)) & ((df['difference_next'] == -1) | (df['difference_next'] ==1)))]
            #Calculate the size of consecutive duplications per gene 
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
            #Select only those genes which contain the defined amount of duplicated exons or which are smaller than the defined number
            df_cons = df_cons[(df_cons['counts']>=dup_size)]
        return df_cons
    

    def filter_for_nonconsecutive_cnvs(self,df:pd.DataFrame,del_or_dup:str,del_size:int,dup_size:int) -> pd.DataFrame :
        """
        Function which extends the prepare_filter_for_consecutive_cnvs function
        It takes the output from the aforementioned function and selects those cnvs, which have a difference to previous/next of 1/-1
        Input Arguments:
            df (pandas DataFrame) : .cnr DataFrame annotated with consecutive cnvs 
            del_or_up (str) : String that determines wether to filter for deletions or duplications 
            del_size (int) : Int that determines the necessary amount of consecutive deletions
            dup_size (int) : Int that determines the necessary amount of consecutive duplications
        Output Arguments:
            df_cons (pandas DataFrame) : .cnr DataFrame with applied filter for consecutive deletions/duplications

        """
        #Set standard value : 
        if del_size is None or del_size=="":
            del_size = 2
        if dup_size is None or dup_size=="":
            dup_size = 2
        del_size = int(del_size)
        dup_size = int(dup_size)
        if del_or_dup == 'del':
            df_cons = df[(((df['call'] == 0) | (df['call'] ==1)) & ((df['difference_previous'] != -1) & (df['difference_previous'] !=1))) & (((df['call'] == 0) | (df['call'] ==1)) & ((df['difference_next'] != -1) & (df['difference_next'] !=1)))]
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
        else:
            df_cons = df[(((df['call'] == 3) | (df['call'] == 4)) & ((df['difference_previous'] != -1) & (df['difference_previous'] !=1))) & (((df['call'] == 3)|(df['call'] == 4)) & ((df['difference_next'] != -1) & (df['difference_next'] !=1)))]
            affected_size = df_cons.groupby('gene')['gene'].size().reset_index(name='counts')
            df_cons = pd.merge(df_cons, affected_size, on='gene', how='left')
        return df_cons
    
    def merge_cnr_with_vcf(self, path_to_bintest_tsv:str, path_to_vcf:str, path_to_output:str):
        """
        Funciton with merges a bintest_tsv or .cnr file with a vcf file
        Input Arguments:
            path_to_bintest_tsv (str) : path to bintest or .cnr file
            path_to_vcf (str) : path to .vcf.gz file
            path_to_output (str) : path to output.vcf
        """

        sample_number = path_to_bintest_tsv.split("/")[-1]
        sample_number = sample_number.split("_")[0] # if bintest file get str before _
        sample_number = sample_number.split(".")[0] # if cnr file get str before .

        #Import input_vcf_df 
        input_vcf_df = self.read_vcf(path_to_vcf)

        #Import Input bintest file 
        input_bintest_df = pd.read_csv(path_to_bintest_tsv, header=0,delimiter="\t")

        #Format bintest 
        input_bintest_df = self.explode_df(input_bintest_df)

        #Seperate bintest in del and dup 
        input_bintest_df_dels = input_bintest_df[(input_bintest_df["call"]==1)|(input_bintest_df["call"]==0)]
        input_bintest_df_dups = input_bintest_df[(input_bintest_df["call"]==3)|(input_bintest_df["call"]==4)]

        #Seperate bintest even further in consecutive and non consecutive 
        consecutive_dels_prep = self.prepare_filter_for_consecutive_cnvs(input_bintest_df_dels)
        consecutive_dups_prep = self.prepare_filter_for_consecutive_cnvs(input_bintest_df_dups)
        consecutive_dels = self.filter_for_consecutive_cnvs(consecutive_dels_prep,"del",2,2)
        consecutive_dups = self.filter_for_consecutive_cnvs(consecutive_dups_prep,"dup",2,2)
        nonconsecutive_dels = self.filter_for_nonconsecutive_cnvs(consecutive_dels_prep,"del",2,2)
        nonconsecutive_dups = self.filter_for_nonconsecutive_cnvs(consecutive_dups_prep,"dup",2,2)

        #Group consecutive_dels and consecutive_dups in individual groups, if a gene contains multiple consecutive cnvs e.g. TTN Exon 1,2,3,4 are deleted and Exon 71,72,73 -> [1,2,3,4] = Group 1 & [71,72,73] = Group 2 
        consecutive_dels["group"] = (consecutive_dels['gene'].ne(consecutive_dels['gene'].shift())|(consecutive_dels['gene'].eq(consecutive_dels['gene'].shift())&consecutive_dels['difference_previous'].ne(consecutive_dels['difference_previous'].shift()))&(consecutive_dels['gene'].eq(consecutive_dels['gene'].shift())&consecutive_dels['difference_next'].ne(consecutive_dels['difference_next'].shift()))).cumsum()
        consecutive_dups["group"] = (consecutive_dups['gene'].ne(consecutive_dups['gene'].shift())|(consecutive_dups['gene'].eq(consecutive_dups['gene'].shift())&consecutive_dups['difference_previous'].ne(consecutive_dups['difference_previous'].shift()))&(consecutive_dups['gene'].eq(consecutive_dups['gene'].shift())&consecutive_dups['difference_next'].ne(consecutive_dups['difference_next'].shift()))).cumsum()

        #Calculate Start and end of consecutive dels / dups
        gene_start_df = consecutive_dels.groupby(['gene',"group"])['start']
        gene_end_df = consecutive_dels.groupby(['gene',"group"])['end']
        consecutive_dels = consecutive_dels.assign(min=gene_start_df.transform(min), max=gene_end_df.transform(max))

        gene_start_df = consecutive_dups.groupby(['gene',"group"])['start']
        gene_end_df = consecutive_dups.groupby(['gene',"group"])['end']
        consecutive_dups = consecutive_dups.assign(min=gene_start_df.transform(min), max=gene_end_df.transform(max))

        #Calculate mean log2 of consecutive dels / dups
        mean_log2_df = consecutive_dels.groupby(['gene',"group"])['log2'].mean().reset_index()
        mean_log2_df = mean_log2_df.rename(columns={"log2":"log2_m"})
        mean_square_df = consecutive_dels.groupby(['gene',"group"])['squaredvalue'].mean().reset_index()
        mean_square_df = mean_square_df.rename(columns={"squaredvalue":"squaredvalue_m"})
        consecutive_dels = pd.merge(consecutive_dels,mean_log2_df,how="left",left_on=["gene","group"],right_on=["gene","group"])
        consecutive_dels = pd.merge(consecutive_dels,mean_square_df,how="left",left_on=["gene","group"],right_on=["gene","group"])

        mean_log2_df = consecutive_dups.groupby(['gene',"group"])['log2'].mean().reset_index()
        mean_log2_df = mean_log2_df.rename(columns={"log2":"log2_m"})
        mean_square_df = consecutive_dups.groupby(['gene',"group"])['squaredvalue'].mean().reset_index()
        mean_square_df = mean_square_df.rename(columns={"squaredvalue":"squaredvalue_m"})
        consecutive_dups = pd.merge(consecutive_dups,mean_log2_df,how="left",left_on=["gene","group"],right_on=["gene","group"])
        consecutive_dups = pd.merge(consecutive_dups,mean_square_df,how="left",left_on=["gene","group"],right_on=["gene","group"])

        #Create a list containing the dataframes that will be concatenated 
        to_be_concatenated = []

        #Prepare consecutive_dels
        if consecutive_dels.empty:
            print("this sample does not contain consecutive deletions")
        else:
            consecutive_dels.loc[((consecutive_dels['call'] == 0)|(consecutive_dels['call'] == 1)), 'CNV_TYPE'] = "DEL"
            consecutive_dels.loc[((consecutive_dels['call'] == 3)|(consecutive_dels['call'] == 4)), 'CNV_TYPE'] = "DUP"
            consecutive_dels.loc[((consecutive_dels['call'] == 0)), 'GT'] = "1/1"
            consecutive_dels.loc[((consecutive_dels['call'] == 1)), 'GT'] = "0/1"
            consecutive_dels.loc[((consecutive_dels['call'] == 3)), 'GT'] = "0/1"
            consecutive_dels.loc[((consecutive_dels['call'] == 4)), 'GT'] = "0/1"

            consecutive_dels["SVLEN"] = consecutive_dels["max"] - consecutive_dels["min"]
            consecutive_dels["#CHROM"] = consecutive_dels["chromosome"]
            consecutive_dels["POS"] = consecutive_dels["min"]
            consecutive_dels["ID"] = "."
            consecutive_dels["REF"] = "N"
            consecutive_dels.loc[((consecutive_dels['call'] == 0)|(consecutive_dels['call'] == 1)), 'ALT'] = "<DEL>"
            consecutive_dels["QUAL"] = "."
            consecutive_dels["FILTER"] = "."
            consecutive_dels["INFO"] = "IMPRECISE;SVTYPE="+consecutive_dels["CNV_TYPE"] + ";END=" + consecutive_dels["max"].astype(str) + ";SVLEN=" + consecutive_dels["SVLEN"].astype(str) + ";FOLD_CHANGE=" + consecutive_dels["squaredvalue_m"].astype(str) + ";FOLD_CHANGE_LOG=" + consecutive_dels["log2_m"].astype(str) + ";PROBES=" + consecutive_dels["counts"].astype(str)
            consecutive_dels.loc[((consecutive_dels['call'] == 0)|(consecutive_dels['call'] == 1)), 'FORMAT'] = "GT:GQ"
            consecutive_dels.loc[((consecutive_dels['call'] == 3)|(consecutive_dels['call'] == 4)), 'FORMAT'] = "GT:GQ:CN:CNQ"
            consecutive_dels[sample_number] = consecutive_dels["GT"] + ":" + consecutive_dels["counts"].astype(str)
            consecutive_dels = consecutive_dels.drop(["chromosome","start","end","gene","log2","depth","weight","p_bintest","exon",
                                                    "call","squaredvalue","difference_previous","difference_next","counts","min","max",
                                                    "CNV_TYPE","SVLEN","group","log2_m","squaredvalue_m","GT"],axis=1)
            consecutive_dels = consecutive_dels.drop_duplicates()
            to_be_concatenated.append(consecutive_dels)

        #Prepare nonconsecutive dels 
        if nonconsecutive_dels.empty:
            print("this sample does not contain nonconsecutive deletions")
        else:
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 0)|(nonconsecutive_dels['call'] == 1)), 'CNV_TYPE'] = "DEL"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 3)|(nonconsecutive_dels['call'] == 4)), 'CNV_TYPE'] = "DUP"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 0)), 'GT'] = "1/1"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 1)), 'GT'] = "0/1"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 3)), 'GT'] = "0/1"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 4)), 'GT'] = "0/1"

            nonconsecutive_dels["SVLEN"] = nonconsecutive_dels["start"] - nonconsecutive_dels["end"]
            nonconsecutive_dels["#CHROM"] = nonconsecutive_dels["chromosome"]
            nonconsecutive_dels["POS"] = nonconsecutive_dels["start"]
            nonconsecutive_dels["ID"] = "."
            nonconsecutive_dels["REF"] = "N"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 0)|(nonconsecutive_dels['call'] == 1)), 'ALT'] = "<DEL>"
            nonconsecutive_dels["QUAL"] = "."
            nonconsecutive_dels["FILTER"] = "."
            nonconsecutive_dels["INFO"] = "IMPRECISE;SVTYPE="+nonconsecutive_dels["CNV_TYPE"] + ";END=" + nonconsecutive_dels["end"].astype(str) + ";SVLEN=" + nonconsecutive_dels["SVLEN"].astype(str) + ";FOLD_CHANGE=" + nonconsecutive_dels["squaredvalue"].astype(str) + ";FOLD_CHANGE_LOG=" + nonconsecutive_dels["log2"].astype(str) + ";PROBES=" + nonconsecutive_dels["counts"].astype(str)
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 0)|(nonconsecutive_dels['call'] == 1)), 'FORMAT'] = "GT:GQ"
            nonconsecutive_dels.loc[((nonconsecutive_dels['call'] == 3)|(nonconsecutive_dels['call'] == 4)), 'FORMAT'] = "GT:GQ:CN:CNQ"
            nonconsecutive_dels[sample_number] = nonconsecutive_dels["GT"] + ":" + nonconsecutive_dels["counts"].astype(str)
            nonconsecutive_dels = nonconsecutive_dels.drop(["chromosome","start","end","gene","log2","depth","weight","p_bintest","exon",
                                                    "call","squaredvalue","difference_previous","difference_next","counts",
                                                    "CNV_TYPE","SVLEN","log2","squaredvalue","GT"],axis=1)
            nonconsecutive_dels = nonconsecutive_dels.drop_duplicates()
            to_be_concatenated.append(nonconsecutive_dels)

        #Prepare nonconsecutive dups 
        if nonconsecutive_dups.empty:
            print("this sample does not contain nonconsecutive duplications")
        else:
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 0)|(nonconsecutive_dups['call'] == 1)), 'CNV_TYPE'] = "DEL"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 3)|(nonconsecutive_dups['call'] == 4)), 'CNV_TYPE'] = "DUP"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 0)), 'GT'] = "1/1"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 1)), 'GT'] = "0/1"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 3)), 'GT'] = "0/1"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 4)), 'GT'] = "0/1"

            nonconsecutive_dups["SVLEN"] = nonconsecutive_dups["start"] - nonconsecutive_dups["end"]

            nonconsecutive_dups["#CHROM"] = nonconsecutive_dups["chromosome"]
            nonconsecutive_dups["POS"] = nonconsecutive_dups["start"]
            nonconsecutive_dups["ID"] = "."
            nonconsecutive_dups["REF"] = "N"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 3)|(nonconsecutive_dups['call'] == 4)), 'ALT'] = "<DUP>"
            nonconsecutive_dups["QUAL"] = "."
            nonconsecutive_dups["FILTER"] = "."
            nonconsecutive_dups["INFO"] = "IMPRECISE;SVTYPE="+nonconsecutive_dups["CNV_TYPE"] + ";END=" + nonconsecutive_dups["end"].astype(str) + ";SVLEN=" + nonconsecutive_dups["SVLEN"].astype(str) + ";FOLD_CHANGE=" + nonconsecutive_dups["squaredvalue"].astype(str) + ";FOLD_CHANGE_LOG=" + nonconsecutive_dups["log2"].astype(str) + ";PROBES=" + nonconsecutive_dups["counts"].astype(str)
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 0)|(nonconsecutive_dups['call'] == 1)), 'FORMAT'] = "GT:GQ"
            nonconsecutive_dups.loc[((nonconsecutive_dups['call'] == 3)|(nonconsecutive_dups['call'] == 4)), 'FORMAT'] = "GT:GQ:CN:CNQ"
            nonconsecutive_dups[sample_number] = nonconsecutive_dups["GT"] + ":" + nonconsecutive_dups["counts"].astype(str)
            nonconsecutive_dups = nonconsecutive_dups.drop(["chromosome","start","end","gene","log2","depth","weight","p_bintest","exon",
                                                    "call","squaredvalue","difference_previous","difference_next","counts",
                                                    "CNV_TYPE","SVLEN","log2","squaredvalue","GT"],axis=1)
            nonconsecutive_dups = nonconsecutive_dups.drop_duplicates()
            to_be_concatenated.append(nonconsecutive_dups)

        #Prepare consecutive dups 
        if consecutive_dups.empty:
            print("this sampe does not contain consecutive duplication")
        else:
            consecutive_dups.loc[((consecutive_dups['call'] == 0)|(consecutive_dups['call'] == 1)), 'CNV_TYPE'] = "DEL"
            consecutive_dups.loc[((consecutive_dups['call'] == 3)|(consecutive_dups['call'] == 4)), 'CNV_TYPE'] = "DUP"
            consecutive_dups.loc[((consecutive_dups['call'] == 0)), 'GT'] = "1/1"
            consecutive_dups.loc[((consecutive_dups['call'] == 1)), 'GT'] = "0/1"
            consecutive_dups.loc[((consecutive_dups['call'] == 3)), 'GT'] = "0/1"
            consecutive_dups.loc[((consecutive_dups['call'] == 4)), 'GT'] = "0/1"

            consecutive_dups["SVLEN"] = consecutive_dups["max"] - consecutive_dups["min"]
            consecutive_dups["#CHROM"] = consecutive_dups["chromosome"]
            consecutive_dups["POS"] = consecutive_dups["min"]
            consecutive_dups["ID"] = "."
            consecutive_dups["REF"] = "N"
            consecutive_dups.loc[((consecutive_dups['call'] == 3)|(consecutive_dups['call'] == 4)), 'ALT'] = "<DUP>"
            consecutive_dups["QUAL"] = "."
            consecutive_dups["FILTER"] = "."
            consecutive_dups["INFO"] = "IMPRECISE;SVTYPE="+consecutive_dups["CNV_TYPE"] + ";END=" + consecutive_dups["max"].astype(str) + ";SVLEN=" + consecutive_dups["SVLEN"].astype(str) + ";FOLD_CHANGE=" + consecutive_dups["squaredvalue_m"].astype(str) + ";FOLD_CHANGE_LOG=" + consecutive_dups["log2_m"].astype(str) + ";PROBES=" + consecutive_dups["counts"].astype(str)
            consecutive_dups.loc[((consecutive_dups['call'] == 0)|(consecutive_dups['call'] == 1)), 'FORMAT'] = "GT:GQ"
            consecutive_dups.loc[((consecutive_dups['call'] == 3)|(consecutive_dups['call'] == 4)), 'FORMAT'] = "GT:GQ:CN:CNQ"
            consecutive_dups[sample_number] = consecutive_dups["GT"] + ":" + consecutive_dups["counts"].astype(str)
            consecutive_dups = consecutive_dups.drop(["chromosome","start","end","gene","log2","depth","weight","p_bintest","exon",
                                                    "call","squaredvalue","difference_previous","difference_next","counts","min","max",
                                                    "CNV_TYPE","SVLEN","group","log2_m","squaredvalue_m","GT"],axis=1)
            consecutive_dups = consecutive_dups.drop_duplicates()
            to_be_concatenated.append(consecutive_dups)

        total_vcf_df = pd.concat(to_be_concatenated)
        #added
        total_vcf_df = total_vcf_df.drop(["probes"],axis=1)
        total_vcf_df = pd.concat([total_vcf_df,input_vcf_df])
        sort_list_for_categorical = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22","chrX","chrY"]
        total_vcf_df["#CHROM"] = pd.Categorical(total_vcf_df["#CHROM"], sort_list_for_categorical)
        total_vcf_df = total_vcf_df.sort_values(by=['#CHROM', 'POS'])
        input_vcf_df_header = self.read_vcf_header(path_to_vcf)

        open(path_to_output + sample_number + "_merged.vcf", "a")
        #for line in input_vcf_df_header[:-1]:
        #    f.write(line)
        total_vcf_df.to_csv(path_to_output + sample_number + "_merged.vcf", sep ="\t", quoting = csv.QUOTE_NONE, header=True, index=False)


        reversed_header = list(reversed(input_vcf_df_header[:-1]))

        with open (path_to_output + sample_number + "_merged.vcf", "r+") as fp:
            all_lines = fp.readlines()
            for line in reversed_header:
                all_lines.insert(0, line)
            fp.seek(0)
            fp.writelines(all_lines)
            
