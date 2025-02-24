
import pandas as pd
import os
import statsmodels.api as sm
from statsmodels.formula.api import ols


def main():
    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/20230824_SVF_MesenCure/"
    method = "ANOVA"
    output_dir = input_dir + "{0}/".format(method)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))

    data = pd.read_csv(input_dir + "python_final_data_logFC0_mc_743_clusters.csv")
    df = data.T
    new_header = df.iloc[0]  # Grab the first row for the header
    df.columns = new_header
    df = df[1:]  # Take the data less the header row
    df["Cluster"] = df["Cluster"].astype(int)

    results = {}

    # Iterate over each row, treating the row as a variable and the columns as the measurements for different groups
    for row_name, row in df.drop('Cluster', axis=1).iterrows():
        # Create a dataframe for each row
        df_row = pd.DataFrame({
            'Value': row,
            'Group': df.columns[:-1],
            'Cluster': df['Cluster']
        })

        model = ols('Value ~ C(Group) * C(Cluster)', data=df_row).fit()
        anova_table = sm.stats.anova_lm(model, typ=2)

        results[row_name] = anova_table['PR(>F)']['C(Group):C(Cluster)']

    # Print results
    for key, value in results.items():
        print(f"ANOVA results for {key}: p-value = {value}")

if __name__ == '__main__':
        main()