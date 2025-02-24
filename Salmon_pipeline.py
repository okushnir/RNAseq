

import os
import sys
import glob

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
from Bio import Entrez
import subprocess
import mygene
import contextlib2
import numpy as np
import requests, json
from io import StringIO
import datetime
import time
from RNAseq_pipeline import ensembl_to_GOannotation


def equal_base_name(input_dir):
    dir_lst = glob.glob(input_dir + "/*")
    for dir in dir_lst:
        basename = dir.split("/")[-1]
        file_lst = glob.glob(dir + "/*")
        for file in file_lst:
            temp_base = file.split("/")[-1].split("R")[0]
            temp_base = temp_base[:-1]
            if basename != temp_base:
                new_file_base = basename + "_R" + file.split("/")[-1].split("R")[-1]
                new_file_path = dir + "/" + new_file_base
                commmad = "mv {0} {1}".format(file, new_file_path)
                print(commmad)
                os.system(commmad)

def equal_dir_file_name_single(input_dir):
    file = glob.glob(input_dir + "/*")[0]
    basename = file.split("/")[-1]
    dir_name = file.split("/")[-2]
    temp_base = file.split("/")[-1].split("R1")[0]
    temp_base = temp_base[:-1]
    if basename != temp_base:

        new_file_path = input_dir.replace(dir_name, temp_base)
        try:
            os.mkdir(new_file_path)
        except OSError:
            print("Creation of the directory {0} failed".format(new_file_path))
        else:
            print("Successfully created the directory {0} ".format(new_file_path))
        new_file_path = new_file_path + "/" + basename
        commmad = "mv {0} {1}".format(file, new_file_path)
        print(commmad)
        os.system(commmad)
        commmad = "rmdir {0}".format(input_dir)
        print(commmad)
        os.system(commmad)


def run_salmon(input_dir, salmon_path, output_dir, file_name, paired=False, L=False):
    try:
        os.mkdir(output_dir)
    except OSError:
        print("Creation of the directory {0} failed".format(output_dir))
    else:
        print("Successfully created the directory {0}".format(output_dir))
    if L:
            intra_dir_lst = glob.glob(input_dir + "/*")
            dir_lst = intra_dir_lst
    else:
        dir_lst = [input_dir]
    for dir in dir_lst:
        if paired:
            with open(salmon_path + "{0}".format(file_name), "w") as shell_script:
                shell_script.write('#!/bin/bash\n'
                                    '   for fn in %s/\n'
                                    'do\n'
                                    'echo "${fn}"\n'
                                    'samp=`basename ${fn}`\n'
                                    'echo "Processing sample ${samp}"\n'
                                    'salmon quant -i ~/Salmon/homo_sapiens_index -l A -1 ${fn}/${samp}_R1_001.fastq.gz -2 ${fn}/${samp}_R2_001.fastq.gz -p 8 --validateMappings -o %s/${samp}_quant\n'
                                    'done\n' % (dir, output_dir))
        else:
            # batch_s = glob.glob(dir + "/*")[0].split("/")[-1]
            # batch_s = dir + "/*".split("/")[-1]
            with open(salmon_path+"{0}".format(file_name), "w") as shell_script:
                shell_script.write('#!/bin/bash\n'
                                    '   for fn in %s/\n'
                                    'do\n'
                                    'echo "${fn}"\n'
                                    'samp=`basename ${fn}`\n'
                                    'echo "Processing sample ${samp}"\n'
                                    'salmon quant -i ~/Salmon/homo_sapiens_index -l A -r ${fn}/${samp}_R1_001.fastq.gz -p 8 --validateMappings -o %s/${samp}_quant\n'
                                    'done\n' % (dir, output_dir))
                                    # -2 ${fn}/${samp}_R2_001.fastq.gz \" % (no_generations))
        command = "bash {0}{1}".format(salmon_path, file_name)
        os.system(command)

def salmon_arrange_paired(salmon_results_path):
    dir_lst = glob.glob(salmon_results_path + "/*")
    dir_lst.sort()
    data_all = pd.DataFrame(columns=["Name"])
    for dir in dir_lst:
        batch = dir.split("/")[-1]
        data_batch = pd.DataFrame(columns=["Name"])
        dir_batch = glob.glob(dir + "/*")
        for sub_dir in dir_batch:
            file = glob.glob(sub_dir + "/*.sf")[0]
            line = sub_dir.split("/")[-1].split("_")[-2]
            data = pd.read_table(file, sep="\t")
            data = data[["Name", "TPM"]]
            data = data.rename(columns={"TPM": "{0}".format(line)})
            data_batch = data_batch.merge(data, how="outer", on="Name")
        data_batch = data_batch.set_index("Name")
        data_batch["Sum"] = data_batch.sum(axis=1)
        data_batch = data_batch.drop(labels=["L001", "L002", "L003", "L004"], axis=1)
        data_batch = data_batch.rename(columns={"Sum": "{0}".format(batch)})
        # data_batch.to_csv(dir + "/{0}_quant.sf".format(batch))
        data_all = data_all.merge(data_batch, how="outer", on="Name")
    data_all = data_all.set_index("Name")
    return data_all
def salmon_arrange(salmon_results_path):
    dir_lst = glob.glob(salmon_results_path + "/*")
    dir_lst.sort()
    data_all = pd.DataFrame(columns=["Name"])
    for dir in dir_lst:
        dir = dir + "/" + dir.split("/")[-1] + "_quant"
        file = glob.glob(dir + "/*.sf")[0]
        try:
            batch = dir.split("/")[-1].split("_")[0]  #+ "_" + dir.split("/")[-1].split("_")[1]
        except IndexError:
            batch = dir.split("/")[-1].split("_")[0]
        if batch == "BFII":
            batch = "BFII.702"
        data = pd.read_table(file, sep="\t")
        data = data[["Name", "TPM"]]
        data = data.rename(columns={"TPM": "{0}".format(batch)})
        data_all = data_all.merge(data, how="outer", on="Name")
    data_all = data_all.set_index("Name")
    return data_all


def unique_deseq(results_path, myGO_path, create_GO=False):
    deseq_data = pd.read_csv(results_path + "DESeq2_Results.csv")#, index_col=0)
    if create_GO:
        data = ensembl_to_GOannotation(deseq_data, index='ensembl_gene_id', arrange=False)
        # data = data.drop_duplicates()#subset=['ensembl_gene_id']
        data.to_csv(results_path + "myGo.csv")
        go_data = data.reset_index()
    else:
        go_data = pd.read_csv(myGO_path + "myGo.csv")
    final_data = pd.merge(deseq_data, go_data, left_on="ensembl_gene_id", right_on="ensembl_gene")
    final_data = final_data.set_index("ensembl_gene_id")
    final_data = final_data.drop(["ensembl_gene"], axis=1)
    final_data = final_data.rename(columns={"Unnamed: 0": "ensembl_transcript"})
    final_data["Gene Ontology (biological process)"] = final_data["Gene Ontology (biological process)"].astype(str)
    final_data["Gene Ontology (biological process)"] = final_data["Gene Ontology (biological process)"].apply(lambda x: str(x.split(" [")[0]))

    final_data_unique = final_data.sort_values(by="log2FoldChange", ascending=False)
    final_data_unique = final_data_unique.drop_duplicates(subset=['hgnc_symbol'])
    final_data_unique = final_data_unique.rename(columns={"Unnamed: 0": "ensembl_transcript"})
    final_data_unique = final_data_unique.sort_values(by="ensembl_transcript", ascending=True)
    final_data_unique = final_data_unique.set_index("ensembl_transcript")
    return final_data_unique


def normalize_read_count_ref_gene(table, batch_no):
    table = table.drop(["Unnamed: 0"], axis=1)
    table_annot = table[["ensembl_transcript", "ensembl_gene_id", "hgnc_symbol"]]
    table_ref = table.loc[table["ensembl_transcript"] == "ENST00000392514.9"]
    table_ref = table_ref.set_index("ensembl_transcript")
    table = table.set_index("ensembl_transcript")
    batch_lst = table.columns.values[1:batch_no+1].tolist()
    table_ref = table_ref[batch_lst]
    table_temp = table[batch_lst]
    table_divide = np.divide(table_temp, table_ref)
    table_divide["sum"] = table_divide.sum(axis=1)
    table_divide = table_divide.sort_values(by="sum", ascending=False)
    table_divide = table_divide.reset_index()
    table_merge = pd.merge(table_divide, table_annot, on="ensembl_transcript", how="left")
    table_merge = table_merge.drop_duplicates(subset=['ensembl_gene_id'])
    column_lst = batch_lst
    column_lst.append("ensembl_transcript")
    table_final = table_merge[column_lst]
    table_final = table_final.set_index("ensembl_transcript")
    table_final["median"] = table_final.median(axis=1)
    table_final = table_final.reset_index()
    table_final = pd.merge(table_final, table_merge, on=column_lst, how="left")
    # quantile = table_final["median"].quantile(q=0.5)
    # table_final = table_final[table_final.quantile(axis=1) > quantile]
    table_final = table_final.set_index("ensembl_transcript")
    return table_final

def main():
    """OptiDonor"""
    """1"""
    input_dir ="/media/bonus/Data/RNASeq_files/CB001/"
    batch_dir = glob.glob(input_dir + "*")
    change_name = False
    if change_name:
        for dir in batch_dir:
            equal_dir_file_name_single(dir)
    salmon_path = "/home/bonus/Salmon/"
    # today = datetime.date.today().strftime("%Y%m%d")
    salmon_results_path = "/media/bonus/Data/RNASeq_files/p4_MesenCure_2020-12-08_analysis/"
    try:
        os.mkdir(salmon_results_path)
    except OSError:
        print("Creation of the directory {0} failed".format(salmon_results_path))
    else:
        print("Successfully created the directory {0} ".format(salmon_results_path))

    for dir in batch_dir:
        batch_output = salmon_results_path + "{0}/".format(dir.split("/")[-1])
        seq_dir = input_dir + "{0}".format(dir.split("/")[-1])
        # run_salmon(seq_dir, salmon_path, batch_output, "python_quant_OptiDonor20240306_samples.sh", paired=True, L=True)
    """2"""
    data = salmon_arrange_paired(salmon_results_path)
    data.to_csv(salmon_results_path + "final_data.csv")
    # data = pd.read_csv(salmon_results_path + "transcription_annotation.csv")
    #
    # normalized_to_ref_data = normalize_read_count_ref_gene(data, batch_no=9)
    # normalized_to_ref_data.to_csv(salmon_results_path + "normalized_to_ref_data.csv")

    """3 - copy final_data to OptiDonor_1 GitHub"""
    # batches1_df = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/OptiDonor_1/AD369-BFII109_final_data.csv")
    # batches2_df = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/OptiDonor_1/AD379-BFII215_final_data.csv")
    # batches3_df = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/OptiDonor_1/AD369_X-MRB005_23_final_data.csv")
    # final_df = pd.merge(batches1_df, batches2_df, on="Name", how="outer")
    # final_df = pd.merge(final_df, batches3_df, on="Name", how="outer")
    # final_df = final_df.set_index(keys="Name")
    # final_df.to_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/OptiDonor_1/final_data.csv")

    """Convert ensembl_transcript to ensembl_gene_id using R script"""

    # df = pd.read_csv("C:/Users/odedku/PycharmProjects/RNASeqProject/Results/OptiDonor_1/final_data_transcription_annotation_var.csv")
    # df = ensembl_to_GOannotation(df, index="ensembl_gene_id", arrange=False)
    # df.to_csv("C:/Users/odedku/PycharmProjects/RNASeqProject/Results/OptiDonor_1/updated_myGO.csv")

    """Naive vs Mesencure"""
    # input_dir ="/home/bonus/RNASeq_files/NaiveVSMesenCure/"
    # batch_dir = glob.glob(input_dir + "*")
    # # for dir in batch_dir:
    # #     equal_base_name(dir)
    # salmon_path = "/home/bonus/Salmon/"
    # today = datetime.date.today().strftime("%Y%m%d")
    # # salmon_results_path = "/home/bonus/RNASeq_files/{0}_NaiveVSMesenCure_analysis/".format(today)
    # salmon_results_path = "/home/bonus/PycharmProjects/RNASeqProject/Results/NaiveVSMesencure/"
    # # try:
    # #     os.mkdir(salmon_results_path)
    # # except OSError:
    # #     print("Creation of the directory {0} failed".format(salmon_results_path))
    # # else:
    # #     print("Successfully created the directory {0} ".format(salmon_results_path))
    # # for dir in batch_dir:
    # #     batch_output = salmon_results_path + "{0}/".format(dir.split("/")[-1])
    # #     try:
    # #         os.mkdir(batch_output)
    # #     except OSError:
    # #         print("Creation of the directory {0} failed".format(batch_output))
    # #     else:
    # #         print("Successfully created the directory {0} ".format(batch_output))
    # #     seq_dir = input_dir + "{0}".format(dir.split("/")[-1])
    # #     run_salmon(seq_dir, salmon_path, batch_output, "python_quant_NaiveVSMesenCure_samples.sh", paired=True)

    # data = salmon_arrange_paried(salmon_results_path)
    # data.to_csv(salmon_results_path + "final_data.csv")


    # my_GO_path = "/home/bonus/PycharmProjects/RNASeqProject/Results/20230119_SVF_p4_analysis/"
    # unique_data = unique_deseq(salmon_results_path, my_GO_path, create_GO=True)
    # unique_data.to_csv(salmon_results_path + "DESeq2_GO_Results_hgnc_unique.csv")

    """SVF vs p4"""
    # svf_data = pd.read_csv("/home/bonus/RNASeq_files/20221212Analysis/20221214SalmonResults/final_data.csv")
    # p4_data = pd.read_csv("/home/bonus/RNASeq_files/20230117_NaiveVSMesenCure_analysis/final_data.csv")
    # p4_data = p4_data.drop(labels=["AD331_X", "CB001_X", "CB002_X"], axis=1)
    # merge_data = pd.merge(svf_data, p4_data, on="Name")
    # merge_data = merge_data.set_index("Name")
    # svf_p4_path = "/home/bonus/RNASeq_files/20230119_SVF_p4_analysis/"
    # # merge_data.to_csv(svf_p4_path + "svf_p4_data.csv")
    # my_GO_path = "/home/bonus/RNASeq_files/20221212Analysis/20221214SalmonResults/"
    # unique_data = unique_deseq(svf_p4_path, my_GO_path, create_GO=True)
    # unique_data.to_csv(svf_p4_path + "DESeq2_GO_Results_hgnc_unique.csv")
    # normailze_data = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/SVF_p4/normalized_data_transcription_annotation.csv")
    # normailze_data = normailze_data.set_index("Unnamed: 0")
    # normailze_data = normailze_data[0:199]
    # g = sns.clustermap(normailze_data, cmap="YlGnBu")
    # plt.savefig("python_heatmap_normailzed_top200.png", dpi=300)
    # plt.close()

    """p0 vs p4 vs MesenCure"""
    # svf_data = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/OptiDonor_1/final_data.csv")
    # p4_mesen_data = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/NaiveVSMesencure/final_data.csv")
    # # p4_data = p4_data.drop(labels=["AD331_X", "CB001_X", "CB002_X"], axis=1)
    # merge_data = pd.merge(svf_data, p4_mesen_data, on="Name")
    # merge_data = merge_data.set_index("Name")
    # p0_p4_mesen_path = "/home/bonus/PycharmProjects/RNASeqProject/Results/p0_p4_Mesencure/"
    # merge_data.to_csv(p0_p4_mesen_path + "p0_p4_mesen_data.csv")
    # my_GO_path = "/home/bonus/RNASeq_files/20221212Analysis/20221214SalmonResults/"
    # # unique_data = unique_deseq(svf_p4_path, my_GO_path, create_GO=True)
    # # unique_data.to_csv(svf_p4_path + "DESeq2_GO_Results_hgnc_unique.csv")
    # # normailze_data = pd.read_csv("/home/bonus/PycharmProjects/RNASeqProject/Results/SVF_p4/normalized_data_transcription_annotation.csv")
    # # normailze_data = normailze_data.set_index("Unnamed: 0")
    # # normailze_data = normailze_data[0:199]
    # # g = sns.clustermap(normailze_data, cmap="YlGnBu")
    # # plt.savefig("python_heatmap_normailzed_top200.png", dpi=300)
    # # plt.close()

if __name__ == '__main__':
    main()