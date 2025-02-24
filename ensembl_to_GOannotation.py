

import sys
import pandas as pd
import mygene
import requests
from io import StringIO
import numpy as np


def process_list(my_list):
    website_api = "https://rest.uniprot.org/"
    my_list = ",".join(my_list)
    r = get_url(f"{website_api}/uniprotkb/accessions?accessions={my_list}&fields=id,accession,go_p,go_c,go_f,cc_function,cc_catalytic_activity,cc_subcellular_location&format=tsv")
    return r


def get_url(url, **kwargs):
    response = requests.get(url, **kwargs)
    if not response.ok:
        print(response.text)
        response.raise_for_status()
        sys.exit()
    return response


def count_file_to_df(count_file):
    """
    :param count_file:
    :return:
    """
    data = pd.read_csv(count_file, header=None, sep="\t")
    data = data.rename(columns={0: "ensembl_gene", 1: "count"})
    return data


def ensembl_to_GOannotation(count_file, index, arrange=True, mg_scopes="ensembl.gene", fields='symbol, uniprot'):
    if arrange:
        data = count_file_to_df(count_file)
    else:
        data = count_file
    if mg_scopes == "ensembl.transcript":
        data["ensembl.transcript"] = data["ensembl.transcript"].apply(lambda x: str(x).split(".")[0])
    data = data.set_index(index)


    # Debug
    # data = data[0:2000]
    mg = mygene.MyGeneInfo()
    ens_lst = data.index.tolist()
    ginfo = mg.querymany(ens_lst, scopes=mg_scopes, fields=fields, species='human')
    ginfo_df = pd.DataFrame.from_dict(ginfo)
    ginfo_df["swiss-prot"] = ginfo_df["uniprot"].apply(lambda x: x.get("Swiss-Prot") if type(x) == dict else "RNA")
    ginfo_df["swiss-prot"] = ginfo_df["swiss-prot"].apply(lambda x: x[-1] if type(x) == list else x)
    ginfo_df = ginfo_df.rename(columns={"query": mg_scopes.replace(".", "_")})
    ginfo_df = ginfo_df.rename(columns={"_id": "ncbi_entrez_gene_id"})
    ginfo_df = ginfo_df.drop(["_score", "uniprot"], axis=1)

    data_proteins = ginfo_df[ginfo_df["swiss-prot"] != "RNA"]
    data_proteins = data_proteins[data_proteins["swiss-prot"].notnull()]
    data_RNA = ginfo_df[ginfo_df["swiss-prot"].isnull()]
    data_out = ginfo_df[ginfo_df["swiss-prot"] == "RNA"]
    data_RNA = pd.concat([data_out, data_RNA], axis=0)
    uniprot_lst = data_proteins["swiss-prot"].tolist()
    data_uniprot_csv_all = pd.DataFrame()
    for i in range(0, len(uniprot_lst), 1000):
        data_uniprot = process_list(uniprot_lst[i:i + 1000])
        print("querying {0}...{1}".format(i+1, i+1000))
        data_uniprot_tsv = data_uniprot.text
        data_uniprot_csv = pd.read_table(StringIO(data_uniprot_tsv), sep="\t", header=0)
        data_uniprot_csv_all = data_uniprot_csv_all.append(data_uniprot_csv)
    data_uniprot_csv_all["Function [CC]"] = data_uniprot_csv_all["Function [CC]"].astype(str)
    data_uniprot_csv_all["Catalytic activity"] = data_uniprot_csv_all["Catalytic activity"].astype(str)
    data_uniprot_csv_all["Subcellular location [CC]"] = data_uniprot_csv_all["Subcellular location [CC]"].astype(str)
    data_uniprot_csv_all["Function [CC]"] = data_uniprot_csv_all["Function [CC]"].apply(lambda x: x.replace("FUNCTION: ", ""))
    data_uniprot_csv_all["Catalytic activity"] = data_uniprot_csv_all["Catalytic activity"].apply(lambda x: x.replace("CATALYTIC ACTIVITY: ", ""))
    data_uniprot_csv_all["Subcellular location [CC]"] = data_uniprot_csv_all["Subcellular location [CC]"].apply(
        lambda x: x.replace("SUBCELLULAR LOCATION: ", ""))
    data_uniprot_csv_all = data_uniprot_csv_all.rename(columns={"Entry": "swiss-prot"})
    ginfo_df_final = pd.merge(data_proteins, data_uniprot_csv_all, on="swiss-prot", how="inner")
    ginfo_df_final = pd.concat([ginfo_df_final, data_RNA], axis=0)
    ginfo_df_final = ginfo_df_final.rename(columns={"Function [CC]": "Function", "Subcellular location [CC]": "Subcellular location"})
    ginfo_df_final = ginfo_df_final.sort_values(by=mg_scopes, ascending=True)
    ginfo_df_final = ginfo_df_final.set_index(mg_scopes)
    return ginfo_df_final


def ensembl_to_hgnc(data, ensembl_id="ensembl_gene_id", scopes="ensembl.gene"):
    genes = data[ensembl_id].tolist()
    if (ensembl_id == "ensembl_transcript") | (ensembl_id == "ensembl_gene_id"):
        genes = [id.split(".")[0] for id in genes]
    mg = mygene.MyGeneInfo()
    hgnc = mg.querymany(genes, scopes=scopes, fields='symbol', species='human')
    hgnc_df = pd.DataFrame.from_dict(hgnc)
    hgnc_df = hgnc_df.rename(columns={"query": scopes})
    hgnc_df = hgnc_df.rename(columns={"_id": "ncbi_entrez_gene_id"})
    hgnc_df = hgnc_df.rename(columns={"symbol": "hgnc_symbol"})
    hgnc_df = hgnc_df.drop(["_score"], axis=1)
    # if ensembl_id == "ensembl_transcript":
    data["ensembl_full"] = data[ensembl_id]
    data[ensembl_id] = data[ensembl_id].apply(lambda x: x.split(".")[0])
    final_data_svf = pd.merge(data, hgnc_df, left_on=ensembl_id, right_on=scopes, how="left")
    # else:
    #     final_data_svf = pd.merge(data, hgnc_df, left_on=ensembl_id, right_on=scopes, how="left")
    # final_data_svf = final_data_svf.drop([scopes, "ncbi_entrez_gene_id"], axis=1)
    final_data_svf["hgnc_symbol"] = final_data_svf["hgnc_symbol"].fillna("None")
    final_data_svf["ncbi_entrez_gene_id"] = final_data_svf["ncbi_entrez_gene_id"].fillna("None")
    final_data_svf["ncbi_entrez_gene_id"] = np.where(final_data_svf["ncbi_entrez_gene_id"] == "None",
                                                      final_data_svf["ensembl_full"],
                                                      final_data_svf["ncbi_entrez_gene_id"])
    final_data_svf["hgnc_symbol"] = np.where(final_data_svf["hgnc_symbol"] == "None",
                                             final_data_svf["ncbi_entrez_gene_id"], final_data_svf["hgnc_symbol"])


    return final_data_svf