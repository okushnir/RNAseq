import pandas as pd


def filter_rows_by_series(df, column, series):
    """
    Filters a DataFrame to keep only rows where the specified column's value is in the given Series.

    Parameters:
    - df: pandas.DataFrame - The DataFrame to filter.
    - column: str - The name of the column to check against the Series.
    - series: pandas.Series - The Series containing values to match in the specified DataFrame column.

    Returns:
    - A new DataFrame containing only the rows where the specified column's value is in the Series.
    """
    # Ensure the column exists in the DataFrame
    if column not in df.columns:
        raise ValueError(f"Column '{column}' does not exist in the DataFrame.")

    # Filter the DataFrame
    filtered_df = df[df[column].isin(series)]

    return filtered_df

def main():
    data = pd.read_csv("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/NaiveVSMesenCure/DESeq2_Results_hgnc_unique.csv")
    angiogenesis_series = pd.read_csv("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/20230824_SVF_MesenCure"
                                      "/angiogenesis_Secreted_log2FC2+.csv")["hgnc_symbol"]
    filtered_df = filter_rows_by_series(data, "hgnc_symbol", angiogenesis_series)
    filtered_df.to_csv("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/20230824_Naive_MesenCure/angiogenesis_naive_MesenCure_filtered.csv", index=False)
if __name__ == '__main__':
    main()