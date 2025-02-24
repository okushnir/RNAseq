import pandas as pd
import torch
from transformers import BioGptTokenizer, BioGptForCausalLM

def find_genes_function_in_biogpt(genes_list):
    # Choose the model and tokenizer
    model_name = "microsoft/biogpt"
    tokenizer = BioGptTokenizer.from_pretrained(model_name)
    model = BioGptForCausalLM.from_pretrained(model_name)

    # Move the model to the CUDA device if available
    device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
    # device = torch.device("cpu")
    print(f"The device running this function is: {device}")
    model.to(device)

    # Generate gene function descriptions
    for gene in genes_list:
        prompt = f"What is the function of {gene}?"

        # Tokenize the prompt
        inputs = tokenizer.encode_plus(prompt, return_tensors="pt")

        # Move the inputs to the CUDA device
        inputs = {k: v.to(device) for k, v in inputs.items()}

        # Generate model output
        with torch.no_grad():
            outputs = model.generate(inputs['input_ids'],
                                     attention_mask=inputs['attention_mask'],
                                     num_return_sequences=1,
                                     min_length=100, max_length=1024, num_beams=5, early_stopping=True)  # Adjust max_length as needed

        # Decode and print the generated results
        decoded_results = tokenizer.decode(outputs[0], skip_special_tokens=True)
        print(f"Gene: {gene}")
        print(f"Description: {decoded_results}\n")


def main():
    # List of genes for which you want to generate descriptions
    genes = pd.read_csv("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_1/cluster7_genes_448_var.csv")
    genes_list = genes["Row.names"].tolist()

    find_genes_function_in_biogpt(genes_list)

if __name__ == '__main__':
    main()



# # tokenizer = BioGptTokenizer.from_pretrained("microsoft/biogpt")
# # model = BioGptForCausalLM.from_pretrained("microsoft/biogpt")
#
# # Load BioGPT model and tokenizer
# model_name = "gpt2"
# tokenizer = GPT2Tokenizer.from_pretrained(model_name)
# model = GPT2LMHeadModel.from_pretrained(model_name)
#
# device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
# model.to(device)
# inputs = {}
# inputs = {k: v.to(device) for k, v in inputs.items()}
#
# set_seed(42)
# genes = pd.read_csv("C:/Users/odedku/PycharmProjects/RNAseqProject/Results/OptiDonor_1/cluster7_genes_448_var.csv")
# genes = genes["Row.names"].tolist()
#
# sentences = [f"The function of {gene} is" for gene in genes]
# #
# #
# # for sentence in sentences:
# #     inputs = tokenizer(sentence, return_tensors="pt")
# #     with torch.no_grad():
# #         beam_output = model.generate(**inputs,
# #                                     min_length=100,
# #                                     max_length=1024,
# #                                     num_beams=5,
# #                                     early_stopping=True
# #                                     )
# #     print(tokenizer.decode(beam_output[0], skip_special_tokens=True))
#
# # Generate hypotheses for each gene
# for gene in genes:
#     # Prompt BioGPT with a specific question about the gene
#     prompt = f"What is the role of {gene} in disease development and progression?"
#
#     # Tokenize the prompt
#     input_ids = tokenizer.encode(prompt, return_tensors="pt")
#
#     # Generate model output
#     with torch.no_grad():
#         output = model.generate(input_ids, max_length=100, num_return_sequences=1)
#
#     # Decode and print the generated hypotheses
#     hypotheses = [tokenizer.decode(ids, skip_special_tokens=True) for ids in output]
#     print(f"Hypotheses for {gene}:")
#     for hypothesis in hypotheses:
#         print(hypothesis)
#     print()