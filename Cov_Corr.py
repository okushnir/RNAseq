import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import numpy as np

def main():
    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_1/characteristics/"
    output_dir = input_dir + "ScoringPotency_FGF/"
    test_data = pd.read_csv(output_dir + "data.csv")
    test_data = test_data[["Doubling Time (mean)", "%Viable-day 0", "Autophagy",
                           "log(IDO activity)", "IC50", "FGF-7"]]
    # test_data = test_data.rename(columns={"Autophagy (RFU per 1Mio cells)": "Autophagy"})
    cov_data = test_data.cov()
    corr_data = test_data.corr()

    # z_scaled = (cov_data.values - np.mean(cov_data.values, axis=1, keepdims=True)) / np.std(cov_data.values, axis=1, keepdims=True)
    # mat_z_var = pd.DataFrame(z_scaled).T
    # mat_z_var.columns = cov_data.columns
    # mat_z_var.index = cov_data.index
    plt1 = sns.heatmap(cov_data, cmap=sns.cubehelix_palette(as_cmap=True), annot=True)
    plt.tight_layout()
    plt.savefig(output_dir + "cov.png")
    plt.close()

    plt2 = sns.heatmap(corr_data, cmap=sns.diverging_palette(25, 25, l=65, center="light", as_cmap=True),
                                                             vmin=-1, vmax=1, annot=True)
    plt.tight_layout()
    plt.savefig(output_dir + "corr.png")
    plt.close()

    # plt3 = sns.pairplot(data=test_data, height=2, kind="reg", diag_kind='kde')
    # plt.autoscale()
    # plt.savefig(output_dir + "SVF_characteristics.png")
    # plt.close()


if __name__ == '__main__':
    main()