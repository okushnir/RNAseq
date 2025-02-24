
import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn import metrics
from sklearn.model_selection import cross_val_score, KFold
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
import shap
from sklearn.tree import plot_tree
import seaborn as sns
import numpy as np
import scipy
from sklearn.tree import DecisionTreeClassifier, plot_tree
from sklearn.tree import export_graphviz
import graphviz
from relationship_plots import assign_grades_highest_best
import datetime
from Potency_genes_linearity_STAR import ensembl_to_hgnc
from relationship_plots import assign_grades_highest_best



def remove_rows_below_threshold_entire_df(dataframe, threshold):
    """
    Removes rows from the DataFrame where any value in the row falls below a given threshold.

    Parameters:
    - dataframe (pd.DataFrame): The DataFrame to filter.
    - threshold (float or int): The threshold value.

    Returns:
    - pd.DataFrame: A DataFrame with rows removed where any value is below the threshold.
    """
    # Applying the condition across all columns and rows, then filtering
    filtered_dataframe = dataframe[dataframe.apply(lambda row: row >= threshold, axis=1).all(axis=1)]
    return filtered_dataframe


def find_best_max_depth(X, y, max_depth_range):
    # Create an empty list to store the cross-validation scores
    cv_scores = []

    # Loop over each max_depth value and compute the cross-validation score
    for max_depth in max_depth_range:
        clf = DecisionTreeClassifier(max_depth=max_depth)
        scores = cross_val_score(clf, X, y, cv=KFold(n_splits=5))
        cv_scores.append(scores.mean())

    # Find the index of the best max_depth value
    best_idx = cv_scores.index(max(cv_scores))

    # Get the best max_depth value and its corresponding score
    best_max_depth = max_depth_range[best_idx]
    best_score = cv_scores[best_idx]

    # Return the best max_depth value and its corresponding score
    return best_max_depth, best_score


def run_the_model(X, y):
    X_trainset, X_testset, y_trainset, y_testset = train_test_split(X, y, test_size=0.2)
    best_max_depth, best_score = find_best_max_depth(X, y, max_depth_range=range(1, 11))
    print("***********************")
    print("best max depth:", best_max_depth)
    print("best score:", best_score)
    clusterTree = DecisionTreeClassifier(criterion="entropy", max_depth=best_max_depth).fit(X_trainset, y_trainset)
    # # Perform 10-fold cross-validation
    # cv = LeaveOneOut()
    # scores = cross_val_score(clusterTree, X, y, cv=cv)
    # print("Accuracy: %0.2f (+/- %0.2f)" % (scores.mean(), scores.std() * 2))
    predTree = clusterTree.predict(X_testset)
    decisionTree_accuracy = metrics.accuracy_score(y_testset, predTree)
    importances = clusterTree.feature_importances_
    count_score = 0
    # for feature_name, importance in zip(data.columns, importances):
    #     print(feature_name, importance)
    #     if importance > 0:
    #         count_score += 1
    return clusterTree, predTree, decisionTree_accuracy, count_score, y_testset

def check_tree_structure(decision_tree, n_classes=4):
    tree = decision_tree.tree_

    def is_leaf(node):
        # A node is a leaf if it doesn't have any children
        return (tree.children_left[node] == -1 and tree.children_right[node] == -1)

    def get_leaves(node):
        # If it's a leaf, return the leaf
        if is_leaf(node):
            return [node]
        # Otherwise, get leaves from both left and right children
        return get_leaves(tree.children_left[node]) + get_leaves(tree.children_right[node])

    # Get all the leaves of the tree
    leaf_nodes = get_leaves(0)

    # For each leaf node, find out the predicted class
    leaf_classes = [np.argmax(tree.value[node]) for node in leaf_nodes]

    # Generate the expected sequence of classes
    expected_classes = np.arange(len(leaf_classes)) % n_classes

    # Check if the sequences match
    return np.array_equal(leaf_classes, expected_classes)


def run_DecisionTree(X, y, test_data, features, output_dir, test_column, accuracy):
    decisionTree_accuracy = 0
    count_score = 0
    features.remove(test_column)
    data = test_data.drop([test_column], axis=1)
    while ((decisionTree_accuracy <= accuracy)):# | (count_score<2)):
        clusterTree, predTree, decisionTree_accuracy, count_score, y_testset = run_the_model(X, y)
        print("DecisionTrees's Accuracy:", decisionTree_accuracy)
    print("***********The algorithm final result:**************")
    print("The predicted-set is: {0}".format(predTree.tolist()))
    print("While the test-set is: {0}".format(y_testset.to_list()))
    print("DecisionTree's Accuracy:", decisionTree_accuracy)


    # Export the decision tree to DOT format
    # dot_data = export_graphviz(clusterTree, out_file=None,
    #                            feature_names=data.columns,
    #                            filled=True, rounded=True,
    #                            special_characters=True, node_ids=True)
    # # Create a graph from the DOT data
    # graph = graphviz.Source(dot_data)
    # # Render the graph
    # graph.render(output_dir + "decision_tree_{0}".format(test_column))

    # tree_structure_ok = False
    # decisionTree_accuracy = 0
    # features.remove(test_column)
    # data = test_data.drop([test_column], axis=1)
    #
    # while not tree_structure_ok or decisionTree_accuracy <= 0.5:
    #     clusterTree, predTree, decisionTree_accuracy, count_score, y_testset = run_the_model(X, y, data)
    #
    #     # Check if the current tree has the desired structure
    #     tree_structure_ok = check_tree_structure(clusterTree)
    #     if tree_structure_ok:
    #         # Check if the tree meets the accuracy threshold
    #         if decisionTree_accuracy > 0.5:
    #             print("Satisfactory Decision Tree found.")
    #             break
    #         else:
    #             print("Decision Tree does not meet the accuracy threshold. Retraining...")
    #     else:
    #         print("Decision Tree does not have the desired structure. Retraining...")

    plt.figure(figsize=(20, 10))
    plot_tree(clusterTree, filled=True, feature_names=data.columns, rounded=True,
              proportion=False, precision=2)
    plt.savefig(output_dir + "decision_tree.png")

    print(test_column)
    return predTree, decisionTree_accuracy


def run_RFC(X, y, data, test_data, features, test_size, output_dir, n_estimators=100, test_column="Cluster"):
    X_trainset, X_testset, y_trainset, y_testset = train_test_split(X, y, test_size=test_size, random_state=52)

    # Create a Random Forest classifier object
    rf_classifier = RandomForestClassifier(n_estimators=n_estimators, random_state=52, criterion="log_loss")
    # Train the Random Forest classifier
    rf_classifier.fit(X_trainset, y_trainset)

    # Make predictions on the test set
    y_pred = rf_classifier.predict(X_testset)

    # Print feature importances
    feature_lst = []
    feature_importances = rf_classifier.feature_importances_
    for feature, importance in zip(test_data.columns, feature_importances):
        print("{0}: {1:.2g}".format(feature, importance))


    # Get SHAP values for feature interaction analysis
    explainer = shap.TreeExplainer(rf_classifier)
    shap_values = explainer.shap_values(X_testset)

    # Plot summary plot of feature importance and interaction effects
    summary_plot = shap.summary_plot(shap_values, X_testset, feature_names=features, plot_type='bar')
    plt.autoscale()
    plt.savefig(output_dir + "summary_plot.png")
    for i in range(0, n_estimators, 1):
        # Get one of the decision trees from the Random Forest
        tree_to_plot = rf_classifier.estimators_[i]

        # Plot the decision tree
        fig, ax = plt.subplots(figsize=(20, 12))
        p1 = plot_tree(tree_to_plot, feature_names=features, filled=True, ax=ax)
        plt.autoscale()
        plt.savefig(output_dir + "RandForest_tree{0}.png".format(i))

    # Evaluate the accuracy of the classifier
    accuracy = accuracy_score(y_testset, y_pred)

    # model_data = test_data.copy()
    # model_data[test_column] = model_data[test_column].fillna(0)
    # model_data = model_data[model_data[test_column] == 0]
    # features.remove(test_column)
    # input_data = model_data[features]
    # input_data = input_data.dropna(axis="rows")
    # predicted_cluster = rf_classifier.predict(input_data)
    # features1 = features.copy()
    # features1.insert(0, "Batch No.")
    # model_data = model_data[features1]
    # model_data = model_data.dropna(axis="rows")
    # model_data = model_data.reset_index()
    # model_data[test_column] = pd.Series(predicted_cluster)
    # model_data = model_data[["Batch No.", test_column]]
    # model_data = pd.merge(data, model_data, on="Batch No.", how="outer")
    # model_data["{}_x".format(test_column)] = model_data["{}_x".format(test_column)].fillna(0)
    # model_data["{}_y".format(test_column)] = model_data["{}_y".format(test_column)].fillna(0)
    # model_data[test_column] = model_data["{}_x".format(test_column)] + model_data["{}_y".format(test_column)]
    # model_data = model_data.drop(["{}_x".format(test_column), "{}_y".format(test_column)], axis=1)
    # features2 = features1.copy()
    # features2.append(test_column)
    # model_data = model_data[features2]
    # model_data.to_csv(output_dir + "model_result.csv")
    return y_pred, accuracy

def main():

    """Gene VS Potency"""
    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/"
    r_input_dir = input_dir + "20240228R_outputs_with_groups/"
    today = datetime.date.today().strftime("%Y%m%d")
    method = "Decision_Tree" #"Random_Forest" #
    filter = "CV_sum" #"linear_reg" #
    top = "top100"
    cv_min_threshold = 0.5
    fc_min_threshold = 5
    count_threshold = 20
    accuracy = 0.7
    drop_rows = True
    output_dir = r_input_dir + "{0}_{1}_{2}_acc{3}_fc{4}_no/".format(today, method, top, accuracy, fc_min_threshold)
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))
    gene_data = pd.read_csv(r_input_dir + "20240325_Linear_reg_Autophagy_IDO activity_IC50/" + "final_data_svf_combat_seq_hgnc_symbol_by_r.csv")
    gene_data = gene_data.drop(columns=["AD374", "BFII.110", "BFII.215", "BFII.216", "MRB003.23", "MRB005.23", "AD388",
                                        "AD389", "AD390", "AD392", "AD393", "BFII.702", "MRB006.23"])
    # gene_data = gene_data.rename(columns={"Unnamed: 0": "hgnc_symbol"})
    # gene_data = ensembl_to_hgnc(gene_data)
    gene_data = gene_data.drop_duplicates(subset=["hgnc_symbol"])
    gene_data = gene_data.set_index("hgnc_symbol")
    gene_data = gene_data.drop(["ensembl_transcript"], axis=1)
    gene_data = 10**(gene_data)

    if filter == "CV_sum":
        """Filtering on CV and sum"""
        # gene_data = gene_data.rename(columns={"Unnamed: 0": "Batch No."})
        # gene_data = gene_data.set_index("Batch No.")
        # gene_data = remove_rows_below_threshold_entire_df(gene_data, 10)
        no_samples = len(gene_data.columns)
        gene_data["sum"] = gene_data.sum(axis=1)
        gene_data = gene_data[gene_data["sum"] > 20 * no_samples]
        gene_data["min"] = gene_data.iloc[:, 0:no_samples].min(axis=1)
        gene_data["max"] = gene_data.iloc[:, 0:no_samples].max(axis=1)
        gene_data["fc"] = gene_data["max"] / gene_data["min"]
        gene_data = gene_data.sort_values(by="fc", ascending=False)
        gene_data["sd"] = gene_data.iloc[:, 0:no_samples].std(axis=1)
        gene_data["mean"] = gene_data.iloc[:, 0:no_samples].mean(axis=1)
        gene_data["cv"] = gene_data["sd"] / gene_data["mean"]
        gene_data["var"] = gene_data.iloc[:, 0:no_samples].var(axis=1)
        gene_data["r"] = gene_data["mean"] ** 2 / gene_data["var"] - gene_data["mean"]

        gene_data = gene_data[gene_data["fc"] >= fc_min_threshold]

        gene_data = gene_data.sort_values(by="r", ascending=False)
        gene_data.to_csv(output_dir + "gene_data.csv")
        gene_data = gene_data.drop(["sum", "min", "max", "fc", "sd", "mean", "cv", "var", "r"], axis=1)
        # gene_data = remove_rows_below_threshold_entire_df(gene_data, count_threshold)
        top_percent = int(1 * len(gene_data))
        gene_data = gene_data.head(top_percent)
        gene_data = gene_data.reset_index()

    elif filter == "linear_reg":
        """top 10% from IDO linear regression model"""
        top10_df = pd.read_csv(r_input_dir + "/20240108_Linear_reg_Autophagy_log(IDO activity)_IC50_FGF-7_all_data/significant_df_{0}.csv".format(top))
        top10_df = top10_df["Gene"]
        gene_data = gene_data.merge(top10_df, left_on="Unnamed: 0", right_on="Gene", how="inner")
        gene_data = gene_data.rename(columns={"Unnamed: 0": "Batch No."})
        gene_data = gene_data.drop(["Batch No."], axis=1)
    gene_data = gene_data.rename(columns={"hgnc_symbol": "Batch No."})
    gene_data = gene_data.T
    if drop_rows:
        gene_data = gene_data.reset_index()
        # Set the first row as column headers
        gene_data.columns = gene_data.iloc[0]
        # Drop the first row from the DataFrame
        gene_data = gene_data.drop(gene_data.index[0])
        features = gene_data.columns.drop("Batch No.").tolist()
    else:
        gene_data = gene_data.reset_index()
        gene_data = gene_data.rename(columns={"index": "Batch No."})
        features = gene_data.columns.drop("Batch No.").tolist()



    grades = ["Grade_FGF-7"]#"Grade_IDO activity","Grade_IC50", "Grade_Autophagy", "Grade_FGF-7", Total grade [%]
    for grade in grades:
        potency_data = pd.read_csv(input_dir + "/characteristics/20240207_ScoringPotency_Autophagy_IDO activity_IC50_FGF-7/data.csv")
        potency_data = potency_data[["Batch No.", grade]]
        if grade == "Total grade [%]":
            potency_data = assign_grades_highest_best(potency_data,"Total grade [%]", quantiles=4)
            potency_data = potency_data[["Batch No.", "Grade"]]
            potency_data = potency_data.rename(columns={"Grade": "Total grade [%]"})
        data = pd.merge(gene_data, potency_data, on="Batch No.")
        data = data.set_index("Batch No.")
        test_data = data.dropna()
        X = test_data.drop(columns=grade, axis=1).values
        y = test_data[grade]

        if method == "Decision_Tree":
            predTree, decisionTree_accuracy = run_DecisionTree(X, y, test_data, grades,
                                                                output_dir=output_dir, test_column=grade, accuracy=accuracy)
            print("***********The algorithm final result:**************")
            print("y_pred:", predTree)
            print("yhat_prob", decisionTree_accuracy)
            print(grade)
        elif method == "Random_Forest":
            y_pred, accuracy = run_RFC(X, y, data, test_data[features], features, test_size=0.25, output_dir=output_dir,
                                       n_estimators=100, test_column="Grade")
            print("***********The algorithm final result:**************")
            print("y_pred:", y_pred)
            print("yhat_prob", accuracy)


if __name__ == '__main__':
    main()

