import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np
import datetime
import os
from ensembl_to_GOannotation import ensembl_to_hgnc


def find_significant_differences_df(df, low_potent_columns, high_potent_columns, threshold=25, id='ensembl_transcript'):
    """
    Finds ensembl_transcript rows where the difference between the average values of
    low-potent and high-potent populations is at least a specified threshold percentage,
    returning a DataFrame with these rows and their statistics.

    Parameters:
    - df: DataFrame containing the populations and ensembl_transcript indices.
    - low_potent_columns: List of column names for the low-potent population.
    - high_potent_columns: List of column names for the high-potent population.
    - threshold: Minimum percentage difference between the averages to consider.

    Returns:
    - DataFrame with ensembl_transcript indices, their low-potent and high-potent averages,
      and the percentage difference where the differences meet or exceed the threshold.
    """
    # Calculate the mean across specified columns for low-potent and high-potent populations
    df['low_potent_avg'] = df[low_potent_columns].mean(axis=1)
    df['high_potent_avg'] = df[high_potent_columns].mean(axis=1)

    # Calculate the percentage difference between the two averages
    df['percent_difference'] = np.abs(df['high_potent_avg'] - df['low_potent_avg']) / (
                (df['low_potent_avg'] + df['high_potent_avg']) / 2) * 100

    # Filter for significant differences
    significant_df = df[df['percent_difference'] >= threshold][
        [id, 'low_potent_avg', 'high_potent_avg', 'percent_difference']]

    return significant_df


def separate_populations_dynamic_rowwise(df, low_potent_columns, high_potent_columns, id='ensembl_transcript'):
    """
    Identifies ensembl_transcript rows where there's no overlap between the two populations,
    dynamically adjusting the comparison based on which population has higher values for each row.
    Additionally, calculates max/min ratios where the numerator is always the higher value,
    irrespective of the population.

    Parameters:
    - df: DataFrame containing the populations and ensembl_transcript indices.
    - low_potent_columns: List of column names for the first population.
    - high_potent_columns: List of column names for the second population.

    Returns:
    - DataFrame with ensembl_transcript indices, dynamic percentile comparisons,
      max/min ratios, and sorted by the difference where the specified no-overlap condition is met.
    """
    results = []
    data = df.copy()
    data = data.reset_index()

    for index, row in data.iterrows():
        low_values = row[low_potent_columns].dropna().astype(float)
        high_values = row[high_potent_columns].dropna().astype(float)

        low_mean = low_values.mean()
        high_mean = high_values.mean()

        if low_mean > high_mean:
            low_95 = high_values.quantile(0.95)
            high_05 = low_values.quantile(0.05)
        else:
            low_95 = low_values.quantile(0.95)
            high_05 = high_values.quantile(0.05)

        # Calculate max and min across both populations for each row
        combined_values = pd.concat([low_values, high_values])
        max_value = combined_values.max()
        min_value = combined_values.min()

        # Ensure to use the higher value as the numerator
        # max_min_ratio = max_value / min_value if min_value != 0 else float('inf')  # Handle division by zero
        if low_mean > high_mean:
            mean_ratio = low_mean / high_mean if high_mean != 0 else float('inf')  # Handle division by zero
        else:
            mean_ratio = high_mean / low_mean if low_mean != 0 else float('inf')  # Handle division by zero


        # Check if there's no overlap for this row
        if low_95 < high_05:
            results.append({
                id: row[id],
                'low_potent_95_or_high_05': low_95,
                'high_potent_05_or_low_95': high_05,
                'mean_ratio': mean_ratio,
                'max_value': max_value,
                'min_value': min_value
            })

    # Convert results to DataFrame and sort by the difference for clearer visualization
    no_overlap_df = pd.DataFrame(results)
    no_overlap_df["difference"] = abs(no_overlap_df["high_potent_05_or_low_95"] - no_overlap_df["low_potent_95_or_high_05"])
    no_overlap_df = no_overlap_df.sort_values(by='difference', ascending=False)

    return no_overlap_df

def main():
    data_deseq = False

    # """Potent data from DESeq2"""
    if data_deseq:
        input_dir = ("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/"
                     "20240228R_outputs_with_groups/DESeq2/1/")
        data = pd.read_csv(input_dir + "row_count_data_potency.csv")
        data = data.rename(columns={"Name": "ensembl_transcript"})
        data = data.set_index("ensembl_transcript")
        data = data.drop(["ensemblt", "baseMean", "log2FoldChange", "lfcSE", "stat", "pvalue", "padj", "ensembl_gene_id",
                          "hgnc_symbol"], axis=1)

    # """All data from ComBat-Seq"""
    else:
        potency_test = "IC50"
        potency_df = pd.read_excel(
            "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/characteristics/OptiDonor_potency_table.xlsx", sheet_name="Sheet1")
        potency_df = potency_df.set_index("Batch No.")
        samples_with_potency = potency_df.index.tolist()
        input_dir = ("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/20240228R_outputs_with_groups/")
        today = datetime.date.today().strftime("%Y%m%d")
        method = "SPIDR"
        output_dir = input_dir + "{0}_{1}/".format(today, method)
        os.makedirs(output_dir, exist_ok=True)
        # try:
        #     os.mkdir(output_dir)
        # except OSError:
        #     print("Creation of the directory {0} failed".format(output_dir))
        # else:
        #     print("Successfully created the directory {0}".format(output_dir))
        data = pd.read_csv(input_dir + "final_data_svf_combat_seq_hgnc_grouped_filtered_wo_ens.csv")
        data = data.rename(columns={"Name": "hgnc_symbol"})
        data = data.set_index("hgnc_symbol")

        data = data[samples_with_potency]
        data = data.replace(0, 0.1)
        # no_samples = len(data.columns)
        # data["sum"] = data.sum(axis=1)
        # data = data[data["sum"] > no_samples*5]
        # data = data.drop(["sum"], axis=1)
        # inf = np.log10(0)
        # data = np.log10(data)
        # data = data.loc[(data != inf).all(axis=1), :]
        # data = data.astype(float)
        # data = data.reset_index()
    flag = "low_high"
    green = ["AD377", "AD379", "AD382", "AD383", "BFII.109", "BFII.216"]
    yellow = ['AD373', 'AD376', 'AD380', 'AD384', 'BFII.112',
                          'AD385', 'AD386', 'AD387', 'AD388', 'AD389']
    red = ['AD369', 'AD370', 'AD372', 'BFII.702']
    low_potent_columns = red.copy()# 'AD369', 'AD370', 'AD372', 'BFII.702'
    high_potent_columns = green.copy() #
    if flag == "high":
        data = data.drop(columns=low_potent_columns)
        high_potent_columns = high_potent_columns  # , 'AD373', 'AD386'] , "AD379", "AD386", "BFII.109"
        low_potent_columns = [col for col in data.columns if col not in high_potent_columns]

    if flag == "low":
        low_potent_columns = low_potent_columns#, 'AD373', 'AD386'] , "AD379", "AD386", "BFII.109"
        high_potent_columns = [col for col in data.columns if col not in low_potent_columns]
    # interest_genes = ["ENST00000523622.1", "ENST00000536371.5", "ENST00000607347.1", "ENST00000545555.2",
    #                  "ENST00000527290.1", "ENST00000548933.5", "ENST00000473546.1", "ENST00000648312.1",
    #                  "ENST00000628703.2", "ENST00000397274.6", "ENST00000675418.2", "ENST00000475721.5",
    #                  "ENST00000551952.5", "ENST00000567610.1", "ENST00000515837.7", "ENST00000460701.1",
    #                  "ENST00000678826.1", "ENST00000574431.5", "ENST00000557621.5", "ENST00000443648.6"]
    if flag == "low_high":
        columns_to_calc = low_potent_columns + high_potent_columns
        data = data[columns_to_calc]
    data = data.reset_index()
    # data = data[data["ensembl_transcript"].isin(interest_genes)]
    significant_transcripts = find_significant_differences_df(data, low_potent_columns, high_potent_columns, threshold=50, id="hgnc_symbol")
    # significant_transcripts.to_csv(output_dir + "significant_differences.csv", index=False)
    significant_transcripts = significant_transcripts.sort_values(by='high_potent_avg', ascending=False)
    # significant_transcripts = ensembl_to_hgnc(significant_transcripts, ensembl_id="ensembl_transcript",
    #                                           scopes="ensembl.transcript")
    print(significant_transcripts.head(10).to_string(index=False))
    print(len(significant_transcripts))

    dynamic_no_overlap_df = separate_populations_dynamic_rowwise(data, low_potent_columns, high_potent_columns, id="hgnc_symbol")
    dynamic_no_overlap_df.to_csv(output_dir + "dynamic_no_overlap_batches_{0}_potent.csv".format(flag), index=False)
    print("Transcripts with dynamically adjusted no overlap between populations:")
    print(dynamic_no_overlap_df.to_string())

    dynamic_no_overlap_df_merged = pd.merge(left=data, right=dynamic_no_overlap_df, on="hgnc_symbol", how="right")
    dynamic_no_overlap_df_merged.to_csv(output_dir + "dynamic_no_overlap_count_{0}_potent.csv".format(flag), index=False)

    # dynamic_no_overlap_df = pd.read_csv(input_dir + "dynamic_no_overlap_count_low_potent.csv")
    # dynamic_no_overlap_df_hgnc = ensembl_to_hgnc(dynamic_no_overlap_df, ensembl_id="ensembl_transcript",
    #                                              scopes="ensembl.transcript")
    # dynamic_no_overlap_df_hgnc.to_csv(output_dir + "dynamic_no_overlap_count_low_potent_hgnc.csv", index=False)


if __name__ == '__main__':
    main()
