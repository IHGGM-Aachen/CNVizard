"""
File which contains plotting functions of the CNVizard
@author: Jeremias Krause, Carlos Classen, Matthias Begemann, Florian Kraft
@company: UKA Aachen (RWTH)
@mail: jerkrause@ukaachen.de
"""

import pandas as pd
import streamlit as st
import plotly.graph_objects as go
import seaborn as sns
import cyvcf2
from matplotlib import pyplot as plt
sns.set_theme()

class CNVPlotter:
    """
    Class containing functions for plotting CNV data.
    """

    def __init__(self):
        """
        Constructor of the Class CNVPlotter.
        """

    def df_to_list(self, df: pd.DataFrame):
        """
        Function which transforms the log2 and depth columns to lists, used for plotting.

        Args:
            df (pd.DataFrame): .cnr DataFrame

        Returns:
            list: List containing the log2 values
            list: List containing the depth values
        """
        col_of_log2 = df.groupby("exon")["log2"].apply(list)
        col_of_depth = df.groupby("exon")["depth"].apply(list)
        list_of_log2 = col_of_log2.tolist()
        list_of_depth = col_of_depth.tolist()
        return list_of_log2, list_of_depth
    
    def index_ref_processor(self, df: pd.DataFrame, selected_gene: str):
        """
        Function used to extract the log2 and depth values from the index .cnr DataFrame for a specific gene.

        Args:
            df (pd.DataFrame): Index .cnr DataFrame
            selected_gene (str): Gene for which to filter

        Returns:
            list: List of log2 values from index for selected gene
            list: List of exons from index for selected gene
            list: List of depth values from index for selected gene
        """
        index_gene_df = df[df["gene"] == selected_gene]
        index_gene_df = index_gene_df.sort_values(by=["exon"], ascending=True)
        list_of_log2_index, list_of_depth_index = self.df_to_list(index_gene_df)
        exon_list = index_gene_df["exon"].unique().tolist()
        return list_of_log2_index, exon_list, list_of_depth_index

    def plot_log2_for_gene_precomputed(
        self,
        gene: str,
        df_total: pd.DataFrame,
        reference_df: pd.DataFrame,
        sample_name: str,
    ):
        """
        Function used to create boxplot for log2 values, extracted from the reference DataFrame.
        Subsequently, the individual log2 values from the index cnr are plotted "on top" of the boxplot to show
        how the index log2 values compare to the reference log2 values.

        Args:
            gene (str): Selected gene
            df_total (pd.DataFrame): Index .cnr DataFrame
            reference_df (pd.DataFrame): DataFrame containing the precomputed statistics, created by reference_builder
            sample_name (str): Index sample name, automatically extracted after uploading the index cnr file
        """
        fig = go.Figure()
        selected_gene = reference_df[reference_df["gene"] == gene]
        listed_q1 = selected_gene["q1_log2"].tolist()
        listed_median = selected_gene["median_log2"].tolist()
        listed_q3 = selected_gene["q3_log2"].tolist()
        listed_min = selected_gene["actual_minimum_log2"].tolist()
        listed_max = selected_gene["actual_maximum_log2"].tolist()
        listed_mean = selected_gene["mean_log2"].tolist()
        selected_gene["std_log2"].tolist()
        listed_exons = selected_gene["exon"].tolist()

        fig.add_trace(
            go.Box(
                q1=listed_q1,
                median=listed_median,
                q3=listed_q3,
                lowerfence=listed_min,
                upperfence=listed_max,
                mean=listed_mean,
                x=listed_exons,
                showlegend=False,
                fillcolor="white",
                line={"color": "black"},
            )
        )
        fig.update_layout(title=f"log2-plot-{gene}-{sample_name}")
        fig.update_xaxes(title_text="Exons", tickmode="linear")
        fig.update_yaxes(title_text="log-2 value", range=[-2, 2])
        list_of_log2_index, exon_list, list_of_depth_index = self.index_ref_processor(
            df_total, gene
        )

        for log2_value, exon_name in zip(list_of_log2_index, listed_exons):
            for value in log2_value:
                y_list = [exon_name]
                x_list = [value]
                fig.add_trace(
                    go.Scatter(
                        y=x_list,
                        x=y_list,
                        mode="markers",
                        marker=dict(color="red"),
                        showlegend=False,
                    )
                )

        fig.add_hline(
            y=0.3,
            line_width=1,
            line_dash="dash",
            line_color="blue",
            showlegend=True,
            name="CN>2",
        )
        fig.add_hline(
            y=-0.4,
            line_width=1,
            line_dash="dash",
            line_color="red",
            showlegend=True,
            name="CN<2",
        )
        fig.add_hline(
            y=-1.1,
            line_width=1,
            line_dash="dash",
            line_color="darkred",
            showlegend=True,
            name="CN<1",
        )
        fig.update_layout(
            legend=dict(
                orientation="h", yanchor="middle", y=1.02, xanchor="left", x=0.05
            )
        )

        if len(exon_list) > 30:
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.plotly_chart(fig)

    def plot_depth_for_gene_precomputed(
        self,
        gene: str,
        df_total: pd.DataFrame,
        reference_df: pd.DataFrame,
        sample_name: str,
    ):
        """
        Function used to create boxplot for depth values, extracted from the reference DataFrame.
        Subsequently, the individual depth values from the index cnr are plotted "on top" of the boxplot to show
        how the index depth values compare to the reference depth values.

        Args:
            gene (str): Selected gene
            df_total (pd.DataFrame): Index .cnr DataFrame
            reference_df (pd.DataFrame): DataFrame containing the precomputed statistics, created by reference_builder
            sample_name (str): Index sample name, automatically extracted after uploading the index cnr file
        """
        fig = go.Figure()
        selected_gene = reference_df[reference_df["gene"] == gene]
        listed_q1 = selected_gene["q1_depth"].tolist()
        listed_median = selected_gene["median_depth"].tolist()
        listed_q3 = selected_gene["q3_depth"].tolist()
        listed_min = selected_gene["actual_minimum_depth"].tolist()
        listed_max = selected_gene["actual_maximum_depth"].tolist()
        listed_mean = selected_gene["mean_depth"].tolist()
        selected_gene["std_depth"].tolist()
        listed_exons = selected_gene["exon"].tolist()

        fig.add_trace(
            go.Box(
                q1=listed_q1,
                median=listed_median,
                q3=listed_q3,
                lowerfence=listed_min,
                upperfence=listed_max,
                mean=listed_mean,
                x=listed_exons,
                showlegend=False,
                fillcolor="white",
                line={"color": "black"},
            )
        )
        fig.update_layout(title=f"depth-plot-{gene}-{sample_name}")
        fig.update_xaxes(title_text="Exons", tickmode="linear")
        fig.update_yaxes(title_text="Depth value")
        list_of_log2_index, exon_list, list_of_depth_index = self.index_ref_processor(
            df_total, gene
        )

        for depth_value, exon_name in zip(list_of_depth_index, listed_exons):
            for value in depth_value:
                y_list = [exon_name]
                x_list = [value]
                fig.add_trace(
                    go.Scatter(
                        y=x_list,
                        x=y_list,
                        mode="markers",
                        marker=dict(color="red"),
                        showlegend=False,
                    )
                )

        if len(exon_list) > 30:
            st.plotly_chart(fig, use_container_width=True)
        else:
            st.plotly_chart(fig)
    
    def set_up_merged_frame(
            self,
            chrom_to_selected: str,
            sample_vcf: pd.DataFrame,
            cnr_df: pd.DataFrame,
            cns_df : pd.DataFrame,
                            
    ):
        """
        Function used to filter VCF files for bins defined in the cnr file.
        Used to plot the BAF plot.

        Args:
            chrom_to_selected (String): Chromosome selected by the user 
            sample_vcf (pd.DataFrame): DataFrame created from the VCf
            cnr_df (pd.DataFrame): DataFrame created from the cnr file
        """
        sample_vcf_chr = sample_vcf[sample_vcf["chromosome"]==chrom_to_selected]
        sample_vcf_chr = sample_vcf_chr.dropna()
        filter_frame = pd.DataFrame()
        cnr_chrom_filtered = cnr_df[cnr_df["chromosome"]==chrom_to_selected]
        chrom_cnr = list(cnr_chrom_filtered.chromosome)
        start_cnr = list(cnr_chrom_filtered.start)
        end_cnr = list(cnr_chrom_filtered.end)
        filter_frame["chromosome"] = chrom_cnr
        filter_frame["start"] = start_cnr
        filter_frame["end"] = end_cnr
        filter_frame1 = filter_frame[filter_frame["chromosome"]==chrom_to_selected]
        merged = pd.merge_asof(left=sample_vcf_chr,right=filter_frame1.rename(columns={"start": "coord","chromosome": "chrom"}),
                left_on="start",right_on="coord",direction="backward")
        merged[merged["end"].notna()]
        merged2 = pd.merge_asof(left=merged,right=filter_frame1.rename(columns={"end": "coord2","chromosome": "chrom2"}),
                            left_on="start",right_on="coord2",direction="forward")
        merged2 = merged2[merged2["end"].notna()]
        merged2 = merged2[merged2["coord2"].notna()]
        merged2_f = merged2[merged2["coord"]==merged2["start_y"]]

        #Filter for cns
        cns_df_selected = cns_df[cns_df["chromosome"]==chrom_to_selected] 
        cns_df_selected = cns_df_selected[cns_df_selected["cn1"]!=cns_df_selected["cn2"]]
        return merged2_f,cns_df_selected


    def plot_scatter(
        self,
        cnr_df: pd.DataFrame,
        cns_df: pd.DataFrame,
        selected_scatter : str,
    ):
        """
        Function used to plot a scatterplot based on the log2 values of the cnr file and 
        the segments from the cns file.

        Args:
            cnr_df (pd.DataFrame): DataFrame created from the .cnr file
            cns_df (pd.DataFrame): DataFrame created from the .cns file 
        """
        if selected_scatter == "All":
            #Filter the copy number segment file for only 
            cns_call_filtered = cns_df[cns_df["cn"]!=2]

            #Split cnr file into chromosomes
            cnr_1 = cnr_df[cnr_df["chromosome"]=="chr1"]
            cnr_2 = cnr_df[cnr_df["chromosome"]=="chr2"]
            cnr_3 = cnr_df[cnr_df["chromosome"]=="chr3"]
            cnr_4 = cnr_df[cnr_df["chromosome"]=="chr4"]
            cnr_5 = cnr_df[cnr_df["chromosome"]=="chr5"]
            cnr_6 = cnr_df[cnr_df["chromosome"]=="chr6"]
            cnr_7 = cnr_df[cnr_df["chromosome"]=="chr7"]
            cnr_8 = cnr_df[cnr_df["chromosome"]=="chr8"]
            cnr_9 = cnr_df[cnr_df["chromosome"]=="chr9"]
            cnr_10 = cnr_df[cnr_df["chromosome"]=="chr10"]
            cnr_11 = cnr_df[cnr_df["chromosome"]=="chr11"]
            cnr_12 = cnr_df[cnr_df["chromosome"]=="chr12"]
            cnr_13 = cnr_df[cnr_df["chromosome"]=="chr13"]
            cnr_14 = cnr_df[cnr_df["chromosome"]=="chr14"]
            cnr_15 = cnr_df[cnr_df["chromosome"]=="chr15"]
            cnr_16 = cnr_df[cnr_df["chromosome"]=="chr16"]
            cnr_17 = cnr_df[cnr_df["chromosome"]=="chr17"]
            cnr_18 = cnr_df[cnr_df["chromosome"]=="chr18"]
            cnr_19 = cnr_df[cnr_df["chromosome"]=="chr19"]
            cnr_20 = cnr_df[cnr_df["chromosome"]=="chr20"]
            cnr_21 = cnr_df[cnr_df["chromosome"]=="chr21"]
            cnr_22 = cnr_df[cnr_df["chromosome"]=="chr22"]
            cnr_X = cnr_df[cnr_df["chromosome"]=="chrX"]
            cnr_Y = cnr_df[cnr_df["chromosome"]=="chrY"]

            #Split cns file into chromosomes
            cns_1_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr1"]
            cns_2_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr2"]
            cns_3_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr3"]
            cns_4_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr4"]
            cns_5_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr5"]
            cns_6_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr6"]
            cns_7_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr7"]
            cns_8_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr8"]
            cns_9_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr9"]
            cns_10_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr10"]
            cns_11_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr11"]
            cns_12_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr12"]
            cns_13_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr13"]
            cns_14_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr14"]
            cns_15_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr15"]
            cns_16_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr16"]
            cns_17_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr17"]
            cns_18_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr18"]
            cns_19_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr19"]
            cns_20_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr20"]
            cns_21_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr21"]
            cns_22_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr22"]
            cns_X_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chrX"]
            cns_Y_call = cns_call_filtered[cns_call_filtered["chromosome"]=="chr2Y"]

            redpalette = ["#ff0000","#ff0000"]

            fig, axes = plt.subplots(1, 24, sharey=True, figsize=(15,5))
            fig.suptitle('Genome wide Scatterplot')
            axes[0].set_title('chr1')
            axes[1].set_title('chr2')
            axes[2].set_title('chr3')
            axes[3].set_title('chr4')
            axes[4].set_title('chr5')
            axes[5].set_title('chr6')
            axes[6].set_title('chr7')
            axes[7].set_title('chr8')
            axes[8].set_title('chr9')
            axes[9].set_title('chr10')
            axes[10].set_title('chr11')
            axes[11].set_title('chr12')
            axes[12].set_title('chr13')
            axes[13].set_title('chr14')
            axes[14].set_title('chr15')
            axes[15].set_title('chr16')
            axes[16].set_title('chr17')
            axes[17].set_title('chr18')
            axes[18].set_title('chr19')
            axes[19].set_title('chr20')
            axes[20].set_title('chr21')
            axes[21].set_title('chr22')
            axes[22].set_title('chrX')
            axes[23].set_title('chrY')
            g= sns.scatterplot(ax=axes[0],x="start",y="log2",data=cnr_1,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[0],x="start",y="log2",data=cns_1_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[0],x="end",y="log2",data=cns_1_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

    
            g2= sns.scatterplot(ax=axes[1],x="start",y="log2",data=cnr_2,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[1],x="start",y="log2",data=cns_2_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[1],x="end",y="log2",data=cns_2_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g3= sns.scatterplot(ax=axes[2],x="start",y="log2",data=cnr_3,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[2],x="start",y="log2",data=cns_3_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[2],x="end",y="log2",data=cns_3_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g4= sns.scatterplot(ax=axes[3],x="start",y="log2",data=cnr_4,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[3],x="start",y="log2",data=cns_4_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[3],x="end",y="log2",data=cns_4_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g5= sns.scatterplot(ax=axes[4],x="start",y="log2",data=cnr_5,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[4],x="start",y="log2",data=cns_5_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[4],x="end",y="log2",data=cns_5_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g6= sns.scatterplot(ax=axes[5],x="start",y="log2",data=cnr_6,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[5],x="start",y="log2",data=cns_6_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[5],x="end",y="log2",data=cns_6_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g7= sns.scatterplot(ax=axes[6],x="start",y="log2",data=cnr_7,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[6],x="start",y="log2",data=cns_7_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[6],x="end",y="log2",data=cns_7_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g8= sns.scatterplot(ax=axes[7],x="start",y="log2",data=cnr_8,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[7],x="start",y="log2",data=cns_8_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[7],x="end",y="log2",data=cns_8_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g9= sns.scatterplot(ax=axes[8],x="start",y="log2",data=cnr_9,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[8],x="start",y="log2",data=cns_9_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[8],x="end",y="log2",data=cns_9_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g10= sns.scatterplot(ax=axes[9],x="start",y="log2",data=cnr_10,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[9],x="start",y="log2",data=cns_10_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[9],x="end",y="log2",data=cns_10_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g11= sns.scatterplot(ax=axes[10],x="start",y="log2",data=cnr_11,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[10],x="start",y="log2",data=cns_11_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[10],x="end",y="log2",data=cns_11_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g12= sns.scatterplot(ax=axes[11],x="start",y="log2",data=cnr_12,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[11],x="start",y="log2",data=cns_12_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[11],x="end",y="log2",data=cns_12_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g13= sns.scatterplot(ax=axes[12],x="start",y="log2",data=cnr_13,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[12],x="start",y="log2",data=cns_13_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[12],x="end",y="log2",data=cns_13_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g14= sns.scatterplot(ax=axes[13],x="start",y="log2",data=cnr_14,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[13],x="start",y="log2",data=cns_14_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[13],x="end",y="log2",data=cns_14_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g15= sns.scatterplot(ax=axes[14],x="start",y="log2",data=cnr_15,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[14],x="start",y="log2",data=cns_15_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[14],x="end",y="log2",data=cns_15_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g16= sns.scatterplot(ax=axes[15],x="start",y="log2",data=cnr_16,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[15],x="start",y="log2",data=cns_16_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[15],x="end",y="log2",data=cns_16_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g17= sns.scatterplot(ax=axes[16],x="start",y="log2",data=cnr_17,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[16],x="start",y="log2",data=cns_17_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[16],x="end",y="log2",data=cns_17_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g18= sns.scatterplot(ax=axes[17],x="start",y="log2",data=cnr_18,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[17],x="start",y="log2",data=cns_18_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[17],x="end",y="log2",data=cns_18_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g19= sns.scatterplot(ax=axes[18],x="start",y="log2",data=cnr_19,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[18],x="start",y="log2",data=cns_19_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[18],x="end",y="log2",data=cns_19_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g20= sns.scatterplot(ax=axes[19],x="start",y="log2",data=cnr_20,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[19],x="start",y="log2",data=cns_20_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[19],x="end",y="log2",data=cns_20_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g21= sns.scatterplot(ax=axes[20],x="start",y="log2",data=cnr_21,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[20],x="start",y="log2",data=cns_21_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[20],x="end",y="log2",data=cns_21_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g22= sns.scatterplot(ax=axes[21],x="start",y="log2",data=cnr_22,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[21],x="start",y="log2",data=cns_22_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[21],x="end",y="log2",data=cns_22_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g23= sns.scatterplot(ax=axes[22],x="start",y="log2",data=cnr_X,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[22],x="start",y="log2",data=cns_X_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[22],x="end",y="log2",data=cns_X_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            g24= sns.scatterplot(ax=axes[23],x="start",y="log2",data=cnr_Y,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[23],x="start",y="log2",data=cns_Y_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[23],x="end",y="log2",data=cns_Y_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            try:
                g.legend_.remove()
                g.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom1")
            try:
                g2.legend_.remove()
                g2.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom2")
            try:
                g3.legend_.remove()
                g3.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom3")
            try:
                g4.legend_.remove()
                g4.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom4")
            try:
                g5.legend_.remove()
                g5.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom5")
            try:
                g6.legend_.remove()
                g6.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom6")
            try:
                g7.legend_.remove()
                g7.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom7")
            try:
                g8.legend_.remove()
                g8.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom8")
            try:
                g9.legend_.remove()
                g9.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom9")
            try:
                g10.legend_.remove()
                g10.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom10")
            try:
                g11.legend_.remove()
                g11.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom11")
            try:
                g12.legend_.remove()
                g12.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom12")
            try:
                g13.legend_.remove()
                g13.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom13")
            try:
                g14.legend_.remove()
                g14.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom14")
            try:
                g15.legend_.remove()
                g15.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom15")
            try:
                g16.legend_.remove()
                g16.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom16")
            try:
                g17.legend_.remove()
                g17.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom17")
            try:
                g18.legend_.remove()
                g18.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom18")
            try:
                g19.legend_.remove()
                g19.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom19")
            try:
                g20.legend_.remove()
                g20.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom20")
            try:
                g21.legend_.remove()
                g21.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom21")
            try:
                g22.legend_.remove()
                g22.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom22")
            try:
                g23.legend_.remove()
                g23.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chromX")
            try:
                g24.legend_.remove()
                g24.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chromY")
            plt.ylim(-2.5, 1.5)
            st.pyplot(fig)
        else:
            cns_call_filtered = cns_df[cns_df["cn"]!=2]
            #Split cnr file into chromosomes
            cnr_1 = cnr_df[cnr_df["chromosome"]==selected_scatter]
            #Split cns file into chromosomes
            cns_1_call = cns_call_filtered[cns_call_filtered["chromosome"]==selected_scatter]
            redpalette = ["#ff0000","#ff0000"]
            fig, axes = plt.subplots(1, 1, sharey=True, figsize=(15,5))
            fig.suptitle('Scatterplot for: '+selected_scatter)
            g= sns.scatterplot(x="start",y="log2",data=cnr_1,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(x="start",y="log2",data=cns_1_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(x="end",y="log2",data=cns_1_call,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            try:
                g.legend_.remove()
                g.set(xlabel=None,xticklabels=[])
            except Exception:
                    st.write("relevant chromosome missing in cnr file")
            plt.ylim(-2.5, 1.5)
            st.pyplot(fig)


    def plot_baf(
        self,
        cnr_df: pd.DataFrame,
        cns_df: pd.DataFrame,
        vcf: cyvcf2.VCF,
        selected_scatter: str,
    ):
        """
        Function used to plot a b-allele frequency plot based on the alt frequency values of the vcf file and 
        the regions of the cnr file and the segments from the cns file.

        Args:
            cnr_df (pd.DataFrame): DataFrame created from the .cnr file
            cns_df (pd.DataFrame): DataFrame created from the .cns file 
            vcf_file_df (cyvcf2.VCF): VCF created by cyvcf2
        """
        if selected_scatter == "All":
            (
            chrom,
            pos,
            ref,
            alt,
            af,
            ) = [[] for i in range(5)]

            for variant in vcf:
                chrom.append(variant.CHROM)
                pos.append(variant.POS)
                ref.append(variant.REF)
                alt.append(variant.ALT)
                af.append(max(variant.gt_alt_freqs))

            sample_vcf = pd.DataFrame(
            {
                "chromosome": chrom,
                "start": pos,
                "ref": ref,
                "alt": alt,
                "alt freq": af,
            }
            )
            sample_vcf["alt"] = sample_vcf["alt"].str[0]

            merged2_f_1,cns1 =  self.set_up_merged_frame("chr1",sample_vcf,cnr_df,cns_df)
            merged2_f_2,cns2 =  self.set_up_merged_frame("chr2",sample_vcf,cnr_df,cns_df)
            merged2_f_3,cns3 =  self.set_up_merged_frame("chr3",sample_vcf,cnr_df,cns_df)
            merged2_f_4,cns4 =  self.set_up_merged_frame("chr4",sample_vcf,cnr_df,cns_df)
            merged2_f_5,cns5 =  self.set_up_merged_frame("chr5",sample_vcf,cnr_df,cns_df)
            merged2_f_6,cns6 =  self.set_up_merged_frame("chr6",sample_vcf,cnr_df,cns_df)
            merged2_f_7,cns7 =  self.set_up_merged_frame("chr7",sample_vcf,cnr_df,cns_df)
            merged2_f_8,cns8 =  self.set_up_merged_frame("chr8",sample_vcf,cnr_df,cns_df)
            merged2_f_9,cns9 =  self.set_up_merged_frame("chr9",sample_vcf,cnr_df,cns_df)
            merged2_f_10,cns10 =  self.set_up_merged_frame("chr10",sample_vcf,cnr_df,cns_df)
            merged2_f_11,cns11 =  self.set_up_merged_frame("chr11",sample_vcf,cnr_df,cns_df)
            merged2_f_12,cns12 =  self.set_up_merged_frame("chr12",sample_vcf,cnr_df,cns_df)
            merged2_f_13,cns13 =  self.set_up_merged_frame("chr13",sample_vcf,cnr_df,cns_df)
            merged2_f_14,cns14 =  self.set_up_merged_frame("chr14",sample_vcf,cnr_df,cns_df)
            merged2_f_15,cns15 =  self.set_up_merged_frame("chr15",sample_vcf,cnr_df,cns_df)
            merged2_f_16,cns16 =  self.set_up_merged_frame("chr16",sample_vcf,cnr_df,cns_df)
            merged2_f_17,cns17 =  self.set_up_merged_frame("chr17",sample_vcf,cnr_df,cns_df)
            merged2_f_18,cns18 =  self.set_up_merged_frame("chr18",sample_vcf,cnr_df,cns_df)
            merged2_f_19,cns19 =  self.set_up_merged_frame("chr19",sample_vcf,cnr_df,cns_df)
            merged2_f_20,cns20 =  self.set_up_merged_frame("chr20",sample_vcf,cnr_df,cns_df)
            merged2_f_21,cns21 =  self.set_up_merged_frame("chr21",sample_vcf,cnr_df,cns_df)
            merged2_f_22,cns22 =  self.set_up_merged_frame("chr22",sample_vcf,cnr_df,cns_df)
            merged2_f_X,cnsX =  self.set_up_merged_frame("chrX",sample_vcf,cnr_df,cns_df)
            merged2_f_Y,cnsY =  self.set_up_merged_frame("chrY",sample_vcf,cnr_df,cns_df)

            fig, axes = plt.subplots(1, 24, sharey=True, figsize=(15,5))
            fig.suptitle('Genome wide BAF plot')
            axes[0].set_title('chr1')
            axes[1].set_title('chr2')
            axes[2].set_title('chr3')
            axes[3].set_title('chr4')
            axes[4].set_title('chr5')
            axes[5].set_title('chr6')
            axes[6].set_title('chr7')
            axes[7].set_title('chr8')
            axes[8].set_title('chr9')
            axes[9].set_title('chr10')
            axes[10].set_title('chr11')
            axes[11].set_title('chr12')
            axes[12].set_title('chr13')
            axes[13].set_title('chr14')
            axes[14].set_title('chr15')
            axes[15].set_title('chr16')
            axes[16].set_title('chr17')
            axes[17].set_title('chr18')
            axes[18].set_title('chr19')
            axes[19].set_title('chr20')
            axes[20].set_title('chr21')
            axes[21].set_title('chr22')
            axes[22].set_title('chrX')
            axes[23].set_title('chrY')
            
            redpalette = ["#ff0000","#ff0000"]

            g= sns.scatterplot(ax=axes[0],x="coord",y="alt freq",data=merged2_f_1,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[0],x="start",y="baf",data=cns1,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[0],x="end",y="baf",data=cns1,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g2= sns.scatterplot(ax=axes[1],x="coord",y="alt freq",data=merged2_f_2,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[1],x="start",y="baf",data=cns2,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[1],x="end",y="baf",data=cns2,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g3= sns.scatterplot(ax=axes[2],x="coord",y="alt freq",data=merged2_f_3,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[2],x="start",y="baf",data=cns3,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[2],x="end",y="baf",data=cns3,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g4= sns.scatterplot(ax=axes[3],x="coord",y="alt freq",data=merged2_f_4,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[3],x="start",y="baf",data=cns4,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[3],x="end",y="baf",data=cns4,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g5= sns.scatterplot(ax=axes[4],x="coord",y="alt freq",data=merged2_f_5,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[4],x="start",y="baf",data=cns5,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[4],x="end",y="baf",data=cns5,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g6= sns.scatterplot(ax=axes[5],x="coord",y="alt freq",data=merged2_f_6,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[5],x="start",y="baf",data=cns6,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[5],x="end",y="baf",data=cns6,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g7= sns.scatterplot(ax=axes[6],x="coord",y="alt freq",data=merged2_f_7,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[6],x="start",y="baf",data=cns7,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[6],x="end",y="baf",data=cns7,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g8= sns.scatterplot(ax=axes[7],x="coord",y="alt freq",data=merged2_f_8,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[7],x="start",y="baf",data=cns8,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[7],x="end",y="baf",data=cns8,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g9= sns.scatterplot(ax=axes[8],x="coord",y="alt freq",data=merged2_f_9,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[8],x="start",y="baf",data=cns9,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[8],x="end",y="baf",data=cns9,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g10= sns.scatterplot(ax=axes[9],x="coord",y="alt freq",data=merged2_f_10,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[9],x="start",y="baf",data=cns10,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[9],x="end",y="baf",data=cns10,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g11= sns.scatterplot(ax=axes[10],x="coord",y="alt freq",data=merged2_f_11,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[10],x="start",y="baf",data=cns11,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[10],x="end",y="baf",data=cns11,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g12= sns.scatterplot(ax=axes[11],x="coord",y="alt freq",data=merged2_f_12,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[11],x="start",y="baf",data=cns12,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[11],x="end",y="baf",data=cns12,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g13= sns.scatterplot(ax=axes[12],x="coord",y="alt freq",data=merged2_f_13,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[12],x="start",y="baf",data=cns13,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[12],x="end",y="baf",data=cns13,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g14= sns.scatterplot(ax=axes[13],x="coord",y="alt freq",data=merged2_f_14,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[13],x="start",y="baf",data=cns14,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[13],x="end",y="baf",data=cns14,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g15= sns.scatterplot(ax=axes[14],x="coord",y="alt freq",data=merged2_f_15,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[14],x="start",y="baf",data=cns15,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[14],x="end",y="baf",data=cns15,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g16= sns.scatterplot(ax=axes[15],x="coord",y="alt freq",data=merged2_f_16,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[15],x="start",y="baf",data=cns16,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[15],x="end",y="baf",data=cns16,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g17= sns.scatterplot(ax=axes[16],x="coord",y="alt freq",data=merged2_f_17,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[16],x="start",y="baf",data=cns17,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[16],x="end",y="baf",data=cns17,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g18= sns.scatterplot(ax=axes[17],x="coord",y="alt freq",data=merged2_f_18,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[17],x="start",y="baf",data=cns18,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[17],x="end",y="baf",data=cns18,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g19= sns.scatterplot(ax=axes[18],x="coord",y="alt freq",data=merged2_f_19,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[18],x="start",y="baf",data=cns19,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[18],x="end",y="baf",data=cns19,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g20= sns.scatterplot(ax=axes[19],x="coord",y="alt freq",data=merged2_f_20,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[19],x="start",y="baf",data=cns20,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[19],x="end",y="baf",data=cns20,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g21= sns.scatterplot(ax=axes[20],x="coord",y="alt freq",data=merged2_f_21,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[20],x="start",y="baf",data=cns21,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[20],x="end",y="baf",data=cns21,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g22= sns.scatterplot(ax=axes[21],x="coord",y="alt freq",data=merged2_f_22,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[21],x="start",y="baf",data=cns22,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[21],x="end",y="baf",data=cns22,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g23= sns.scatterplot(ax=axes[22],x="coord",y="alt freq",data=merged2_f_X,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[22],x="start",y="baf",data=cnsX,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[22],x="end",y="baf",data=cnsX,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            g24= sns.scatterplot(ax=axes[23],x="coord",y="alt freq",data=merged2_f_Y,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[23],x="start",y="baf",data=cnsY,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)
            sns.scatterplot(ax=axes[23],x="end",y="baf",data=cnsY,hue="chromosome",palette=redpalette,linewidth=0,alpha=0.5)

            try:
                g.legend_.remove()
                g.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom1")
            try:
                g2.legend_.remove()
                g2.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom2")
            try:
                g3.legend_.remove()
                g3.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom3")
            try:
                g4.legend_.remove()
                g4.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom4")
            try:
                g5.legend_.remove()
                g5.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom5")
            try:
                g6.legend_.remove()
                g6.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom6")
            try:
                g7.legend_.remove()
                g7.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom7")
            try:
                g8.legend_.remove()
                g8.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom8")
            try:
                g9.legend_.remove()
                g9.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom9")
            try:
                g10.legend_.remove()
                g10.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom10")
            try:
                g11.legend_.remove()
                g11.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom11")
            try:
                g12.legend_.remove()
                g12.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom12")
            try:
                g13.legend_.remove()
                g13.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom13")
            try:
                g14.legend_.remove()
                g14.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom14")
            try:
                g15.legend_.remove()
                g15.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom15")
            try:
                g16.legend_.remove()
                g16.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom16")
            try:
                g17.legend_.remove()
                g17.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom17")
            try:
                g18.legend_.remove()
                g18.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom18")
            try:
                g19.legend_.remove()
                g19.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom19")
            try:
                g20.legend_.remove()
                g20.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom20")
            try:
                g21.legend_.remove()
                g21.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom21")
            try:
                g22.legend_.remove()
                g22.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chrom22")
            try:
                g23.legend_.remove()
                g23.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chromX")
            try:
                g24.legend_.remove()
                g24.set(xlabel=None,xticklabels=[])
            except Exception:
                print("no chromY")
            plt.ylim(0, 1)
            st.pyplot(fig)
        else:
            try:
                (
                chrom,
                pos,
                ref,
                alt,
                af,
                ) = [[] for i in range(5)]

                for variant in vcf:
                    chrom.append(variant.CHROM)
                    pos.append(variant.POS)
                    ref.append(variant.REF)
                    alt.append(variant.ALT)
                    af.append(max(variant.gt_alt_freqs))

                sample_vcf = pd.DataFrame(
                {
                    "chromosome": chrom,
                    "start": pos,
                    "ref": ref,
                    "alt": alt,
                    "alt freq": af,
                }
                )
                sample_vcf["alt"] = sample_vcf["alt"].str[0]
                
                sample_vcf_chr = sample_vcf[sample_vcf["chromosome"]==selected_scatter]
                sample_vcf_chr = sample_vcf_chr.dropna()

                filter_frame = pd.DataFrame()
                cnr_chrom_filtered = cnr_df[cnr_df["chromosome"]==selected_scatter]
                
                chrom_cnr = list(cnr_chrom_filtered.chromosome)
                start_cnr = list(cnr_chrom_filtered.start)
                end_cnr = list(cnr_chrom_filtered.end)
                filter_frame["chromosome"] = chrom_cnr
                filter_frame["start"] = start_cnr
                filter_frame["end"] = end_cnr

                filter_frame1 = filter_frame[filter_frame["chromosome"]==selected_scatter]


                merged = pd.merge_asof(left=sample_vcf_chr,right=filter_frame1.rename(columns={"start": "coord","chromosome": "chrom"}),
                      left_on="start",right_on="coord",direction="backward")

                merged[merged["end"].notna()]

                merged2 = pd.merge_asof(left=merged,right=filter_frame1.rename(columns={"end": "coord2","chromosome": "chrom2"}),
                                    left_on="start",right_on="coord2",direction="forward")

                merged2 = merged2[merged2["end"].notna()]
                merged2 = merged2[merged2["coord2"].notna()]
                merged2_f = merged2[merged2["coord"]==merged2["start_y"]]

                fig, axes = plt.subplots(1, 1, sharey=True, figsize=(15,5))
                fig.suptitle('BAF plot for: '+selected_scatter)
                g= sns.scatterplot(x="coord",y="alt freq",data=merged2_f,s=8,hue="chromosome",palette="grey",linewidth=0,alpha=0.5)
                try:
                    g.legend_.remove()
                    g.set(xlabel=None,xticklabels=[])
                except Exception:
                    st.write("relevant chromosome missing in cnr file")
                st.pyplot(fig)

            except Exception as e:
                st.write("no b allele frequency found for chrom: "+selected_scatter)
                st.write(e)

