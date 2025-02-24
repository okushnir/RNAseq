from ensembl_to_GOannotation import ensembl_to_hgnc
import pandas as pd
import requests
import pymysql
import pandas as pd
from sqlalchemy import create_engine


def test_ucsc_api():
    response = requests.get("https://api.genome.ucsc.edu/list/tracks?genome=hg38")
    if response.status_code == 200:
        print("Connection successful, response:", response.json())
    else:
        print("Failed to connect, status code:", response.status_code)


def fetch_gene_symbols_by_ensembl(ensembl_ids):
    results = []
    base_url = "https://api.genome.ucsc.edu/getData/track"
    genome = "hg38"
    track = "ensemblGene"  # This is an example, replace with the correct track from your research

    for ensembl_id in ensembl_ids:
        url = f"{base_url}?genome={genome};track={track};name={ensembl_id}"
        response = requests.get(url)
        if response.status_code == 200:
            data = response.json()
            hgnc_symbol = data.get(ensembl_id, {}).get('geneSymbol')  # Adjust according to the actual data structure
            results.append({'ensembl_id': ensembl_id, 'hgnc_symbol': hgnc_symbol})
        else:
            results.append({'ensembl_id': ensembl_id, 'error': 'No data found', 'status_code': response.status_code, 'response': response.text})
        print(f"Response for {ensembl_id}: {response.status_code}, {response.text}")
    return results


def fetch_hgnc_symbol(transcript_id):
    # Database connection
    connection = pymysql.connect(host='genome-mysql.soe.ucsc.edu',
                                 user='genome',
                                 password='',
                                 database='hg38')

    try:
        with connection.cursor() as cursor:
            # SQL query to fetch HGNC symbol
            sql = f"SELECT geneSymbol FROM kgXref WHERE kgID = '{transcript_id}'"
            cursor.execute(sql)
            result = cursor.fetchone()
            if result:
                return result[0]
            else:
                return "No matching HGNC symbol found."
    except Exception as e:
        print(f"Error: {e}")
    finally:
        connection.close()





def main():
    input_dir = "C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_Salmon_ComBatSeq/20240801R_outputs_with_groups/"
    file_name = "final_data_svf_combat_seq"
    drop_batches = ["AD371", "AD374", "AD371.2", "AD376.2", "AD384.2", "AD387.2"]

    data = pd.read_csv(input_dir + file_name + ".csv")
    columns = data.columns.tolist()
    print(columns[0])
    data = data.rename(columns={f"{columns[0]}": "ensembl_transcript"})
    data = data.drop(drop_batches, axis=1)
    df = ensembl_to_hgnc(data, ensembl_id="ensembl_transcript", scopes="ensembl.transcript")
    df.to_csv(input_dir + file_name + "_hgnc.csv")
    grouped = df.groupby("hgnc_symbol").sum()
    grouped.to_csv(input_dir + file_name + "_hgnc_grouped.csv")

    grouped = pd.read_csv(input_dir + file_name + "_hgnc_grouped.csv")
    filtered_df = grouped[~(grouped['hgnc_symbol'].str.startswith('ENST') | grouped['hgnc_symbol'].str.startswith('ENSG'))]#& (grouped.sum(axis=1) < 100))]
    filtered_df.to_csv(input_dir + file_name + "_hgnc_grouped_filtered_wo_ens.csv", index=False)
    # filtered_df = filtered_df.set_index('hgnc_symbol')
    # filtered_df["sum"] = filtered_df.sum(axis=1)
    # filtered_df = filtered_df[filtered_df['sum'] > 100]
    # filtered_df.reset_index(inplace=True)
    # filtered_df = filtered_df.rename(columns={"hgnc_symbol": "ensembl_transcript"})
    # filtered_df_symbol = ensembl_to_hgnc(filtered_df, ensembl_id="ensembl_transcript", scopes="ensembl.transcript")
    # filtered_df_symbol_grouped = filtered_df_symbol.groupby("hgnc_symbol").sum()
    # filtered_df_symbol_grouped.to_csv(input_dir + file_name + "_hgnc_grouped_ensg.csv")
    # print(filtered_df_symbol.to_string())

    # genes = data["ensembl_transcript"].tolist()
    # # test_ucsc_api()
    # genes_test = ["ENST00000361851.1", "ENST00000510231.1", "ENST00000301905.9", "ENST00000522944.5"]
    # results = fetch_gene_symbols_by_ensembl(genes_test)



    ensembl_transcript_id = 'ENST00000437460.1'  # Replace with your actual transcript ID
    hgnc_symbol = fetch_hgnc_symbol(ensembl_transcript_id)
    print(f"The HGNC symbol for transcript ID {ensembl_transcript_id} is {hgnc_symbol}")

if __name__ == '__main__':
    main()