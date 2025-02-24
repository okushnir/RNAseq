import matplotlib.pyplot as plt
from matplotlib_venn import venn2
import datetime
import os
import pandas as pd
import numpy as np


def high_var_venn_diagram(mc_df, svf_df, output_dir):
    mc_df = mc_df.set_index("hgnc_symbol")
    mc_df = mc_df.drop(["Name"], axis=1)
    mc_df["Var"] = mc_df.var(axis=1)
    mc_df["SD"] = mc_df.std(axis=1)
    mc_df["mean"] = mc_df.mean(axis=1)
    mc_df["sum"] = mc_df.sum(axis=1)
    mc_df["cv"] = mc_df["SD"] / mc_df["mean"]*100
    mc_df = mc_df[mc_df["cv"] < 100]
    mc_df = mc_df[mc_df["cv"] >= 30]
    mc_df = mc_df[mc_df["sum"] > 30]
    mc_df = mc_df.sort_values(by="cv", ascending=True)
    mc_df = mc_df.sort_values(by="sum", ascending=False)
    mc_df = mc_df.reset_index()
    mc_genes = mc_df["hgnc_symbol"]
    mc_set = set(mc_genes)
    svf_df = svf_df.set_index("hgnc_symbol")
    svf_df = svf_df.drop(["Name"], axis=1)
    svf_df["Var"] = svf_df.var(axis=1)
    svf_df["SD"] = svf_df.std(axis=1)
    svf_df["mean"] = svf_df.mean(axis=1)
    svf_df["sum"] = svf_df.sum(axis=1)
    svf_df["cv"] = svf_df["SD"] / svf_df["mean"]*100
    svf_df = svf_df[svf_df["cv"] < 100]
    svf_df = svf_df[svf_df["cv"] >= 30]
    svf_df = svf_df[svf_df["sum"] > 30]
    svf_df = svf_df.sort_values(by="cv", ascending=True)
    svf_df = svf_df.sort_values(by="sum", ascending=False)
    svf_df = svf_df.reset_index()
    svf_genes = svf_df["hgnc_symbol"]
    svf_set = set(svf_genes)

    # Create Venn Diagram
    # Create a figure for the diagram
    plt.figure(figsize=(8, 6))

    venn2([mc_set, svf_set], ('MC-Set', 'SVF-Set'), set_colors=('purple', 'skyblue'), alpha=0.7)
    plt.title("Venn diagram")
    plt.savefig(output_dir + "high_var_Venn.png")

    mc_diff = mc_set - svf_set
    svf_diff = svf_set - mc_set
    mc_svf_eqaul = mc_set - mc_diff
    svf_diff_df = pd.DataFrame(svf_diff)
    mc_diff_df = pd.DataFrame(mc_diff)
    mc_svf_eqaul_df = pd.DataFrame(mc_svf_eqaul)
    svf_diff_df.to_csv(output_dir + "high_var_svf_diff_df.csv", index=False)
    mc_diff_df.to_csv(output_dir + "high_var_mc_diff_df.csv", index=False)
    mc_svf_eqaul_df.to_csv(output_dir + "high_var_mc_svf_eqaul_df.csv", index=False)
    return svf_diff_df, mc_diff_df, mc_svf_eqaul_df


def low_var_venn_diagram(mc_df, svf_df, output_dir):
    mc_df = mc_df.set_index("hgnc_symbol")
    mc_df = mc_df.drop(["Name"], axis=1)
    mc_df["Var"] = mc_df.var(axis=1)
    mc_df["SD"] = mc_df.std(axis=1)
    mc_df["mean"] = mc_df.mean(axis=1)
    mc_df["sum"] = mc_df.sum(axis=1)
    mc_df["cv"] = mc_df["SD"] / mc_df["mean"] * 100
    mc_df = mc_df[mc_df["cv"] < 30]
    mc_df = mc_df[mc_df["sum"] > 30]
    mc_df = mc_df.sort_values(by="cv", ascending=True)
    mc_df = mc_df.sort_values(by="sum", ascending=False)
    mc_df = mc_df.reset_index()
    mc_genes = mc_df["hgnc_symbol"]
    mc_set = set(mc_genes)
    svf_df = svf_df.set_index("hgnc_symbol")
    svf_df = svf_df.drop(["Name"], axis=1)
    svf_df["Var"] = svf_df.var(axis=1)
    svf_df["SD"] = svf_df.std(axis=1)
    svf_df["mean"] = svf_df.mean(axis=1)
    svf_df["sum"] = svf_df.sum(axis=1)
    svf_df["cv"] = svf_df["SD"] / svf_df["mean"] * 100
    svf_df = svf_df[svf_df["cv"] < 30]
    svf_df = svf_df[svf_df["sum"] > 30]
    svf_df = svf_df.sort_values(by="cv", ascending=True)
    svf_df = svf_df.sort_values(by="sum", ascending=False)
    svf_df = svf_df.reset_index()
    svf_genes = svf_df["hgnc_symbol"]
    svf_set = set(svf_genes)

    # Create Venn Diagram
    # Create a figure for the diagram
    plt.figure(figsize=(8, 6))

    venn2([mc_set, svf_set], ('MC-Set', 'SVF-Set'), set_colors=('purple', 'skyblue'), alpha=0.7)
    plt.title("Venn diagram")
    plt.savefig(output_dir + "low_var_Venn.png")

    mc_diff = mc_set - svf_set
    svf_diff = svf_set - mc_set
    mc_svf_eqaul = mc_set - mc_diff
    svf_diff_df = pd.DataFrame(svf_diff)
    mc_diff_df = pd.DataFrame(mc_diff)
    mc_svf_eqaul_df = pd.DataFrame(mc_svf_eqaul)
    svf_diff_df.to_csv(output_dir + "low_var_svf_diff_df.csv", index=False)
    mc_diff_df.to_csv(output_dir + "low_var_mc_diff_df.csv", index=False)
    mc_svf_eqaul_df.to_csv(output_dir + "low_var_mc_svf_eqaul_df.csv", index=False)
    return svf_diff_df, mc_diff_df, mc_svf_eqaul_df

def main():
    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/20230824_SVF_MesenCure/"
    today = datetime.date.today().strftime("%Y%m%d")
    method = "{}_Venn".format(today)
    output_dir = input_dir + "{0}/".format(method)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))

    mc_df = pd.read_csv(input_dir + "python_final_data_mc.csv")
    svf_df = pd.read_csv(input_dir + "python_final_data_svf.csv")
    final_data = pd.read_csv(input_dir + "python_final_data.csv")

    high_svf_diff_df, high_mc_diff_df, high_mc_svf_eqaul_df = high_var_venn_diagram(mc_df, svf_df, output_dir)
    low_svf_diff_df, low_mc_diff_df, low_mc_svf_eqaul_df = low_var_venn_diagram(mc_df, svf_df, output_dir)
    """High Var"""
    high_svf_diff_df = high_svf_diff_df.rename(columns={0: "hgnc_symbol"})
    high_mc_diff_df = high_mc_diff_df.rename(columns={0: "hgnc_symbol"})
    high_mc_svf_eqaul_df = high_mc_svf_eqaul_df.rename(columns={0: "hgnc_symbol"})
    expression_high_svf_diff_df = final_data.merge(high_svf_diff_df, on="hgnc_symbol", how="inner")
    expression_high_mc_diff_df = final_data.merge(high_mc_diff_df, on="hgnc_symbol", how="inner")
    expression_high_mc_svf_eqaul_df = final_data.merge(high_mc_svf_eqaul_df, on="hgnc_symbol", how="inner")

    with pd.ExcelWriter(output_dir + "expression_high_svf_diff_df.xlsx") as writer:
        expression_high_svf_diff_df.to_excel(writer, sheet_name="high_svf_diff", index=False)
        expression_high_mc_diff_df.to_excel(writer, sheet_name="high_mc_diff", index=False)
        expression_high_mc_svf_eqaul_df.to_excel(writer, sheet_name="high_mc_svf_eqaul", index=False)

    """Low Var"""
    low_svf_diff_df = low_svf_diff_df.rename(columns={0: "hgnc_symbol"})
    low_mc_diff_df = low_mc_diff_df.rename(columns={0: "hgnc_symbol"})
    low_mc_svf_eqaul_df = low_mc_svf_eqaul_df.rename(columns={0: "hgnc_symbol"})
    expression_low_svf_diff_df = final_data.merge(low_svf_diff_df, on="hgnc_symbol", how="inner")
    expression_low_mc_diff_df = final_data.merge(low_mc_diff_df, on="hgnc_symbol", how="inner")
    expression_low_mc_svf_eqaul_df = final_data.merge(low_mc_svf_eqaul_df, on="hgnc_symbol", how="inner")

    with pd.ExcelWriter(output_dir + "expression_low_svf_diff_df.xlsx") as writer:
        expression_low_svf_diff_df.to_excel(writer, sheet_name="low_svf_diff", index=False)
        expression_low_mc_diff_df.to_excel(writer, sheet_name="low_mc_diff", index=False)
        expression_low_mc_svf_eqaul_df.to_excel(writer, sheet_name="low_mc_svf_eqaul", index=False)




if __name__ == '__main__':
    main()