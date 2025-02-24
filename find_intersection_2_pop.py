import numpy as np
from sklearn.mixture import GaussianMixture
from scipy.optimize import fsolve
import matplotlib.pyplot as plt
import pandas as pd
import os


def fit_and_find_intersection(data, output_dir, feature):
    """
    Fits a bimodal distribution using a Gaussian Mixture Model and finds the intersection.

    Parameters:
    - data: Array-like, the dataset to analyze.

    Returns:
    - The intersection point(s) of the two Gaussian distributions.
    """
    # Reshape data for the model
    data = np.array(data).reshape(-1, 1)

    # Fit a Gaussian Mixture Model with 2 components
    gmm = GaussianMixture(n_components=2, random_state=0).fit(data)

    # Extract means and covariances
    means = gmm.means_.flatten()
    stds = np.sqrt(gmm.covariances_.flatten())

    # Define the equation to solve: equating two Gaussian functions
    def eq_to_solve(x):
        return (1 / (stds[0] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[0]) / stds[0]) ** 2) - \
            (1 / (stds[1] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[1]) / stds[1]) ** 2)

    # Initial guess for the intersection point
    initial_guess = np.mean(data)

    # Solve for the x that makes the equation zero
    x_intersection = fsolve(eq_to_solve, initial_guess)

    # Plotting for visual confirmation
    x = np.linspace(min(data) - 1, max(data) + 1, 1000)
    y1 = (1 / (stds[0] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[0]) / stds[0]) ** 2)
    y2 = (1 / (stds[1] * np.sqrt(2 * np.pi))) * np.exp(-0.5 * ((x - means[1]) / stds[1]) ** 2)
    plt.plot(x, y1, label='Gaussian 1')
    plt.plot(x, y2, label='Gaussian 2')
    plt.scatter(x_intersection, eq_to_solve(x_intersection) + (1 / (stds[1] * np.sqrt(2 * np.pi))) * np.exp(
        -0.5 * ((x_intersection - means[1]) / stds[1]) ** 2), c='red', label='Intersection')
    plt.hist(data, density=True, bins=15, alpha=0.5, label='Data Histogram')
    plt.legend()
    plt.savefig(output_dir + "histogram_intersection_{0}.png".format(feature))
    plt.close()

    return x_intersection


def find_intersection_2_pop(input_dir, output_dir, potency_lst, remove_samples=None):
    for potency in potency_lst:
        data = pd.read_excel(input_dir + "OptiDonor_master_table.xlsx")
        data = data.set_index("Batch No.")
        if remove_samples is not None:
            data = data.drop(labels=remove_samples, axis=0)
        data = data[[potency]]
        # data = data.set_index("Batch No.")
        data = data.dropna()

        # Find and print the intersection point
        intersection_point = fit_and_find_intersection(data, output_dir, potency)
        print("Intersection point:", intersection_point)

        median_point = data.T.median(axis=1)
        print("Median point:", median_point)

        data_binomial = data.copy()
        if potency == "IC50":
            data_binomial["{0}_Grade".format(potency)] = np.where(data_binomial[potency] < intersection_point[0], 1, 0)
            data_binomial["{0}_Bool".format(potency)] = np.where(data_binomial[potency] < intersection_point[0], "Yes",
                                                                 "No")
        else:
            data_binomial["{0}_Grade".format(potency)] = np.where(data_binomial[potency] > intersection_point[0], 1, 0)
            data_binomial["{0}_Bool".format(potency)] = np.where(data_binomial[potency] > intersection_point[0], "Yes",
                                                                 "No")
        data_binomial.to_csv(input_dir + "sample_info_{0}.csv".format(potency))



def main():
    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/characteristics/"
    output_dir = input_dir + "pop_intersection/"
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))
    potency_lst = ["IC50", "log(IDO activity)", "FGF7"]
    remove_samples = ["AD369_2", "AD373_2", "AD371"]
    find_intersection_2_pop(input_dir, output_dir, potency_lst, remove_samples=remove_samples)


if __name__ == '__main__':
    main()