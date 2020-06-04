import pandas as pd
import numpy as np

# generate table for PCA analysis


def PCA_file(path_in):
    seq_df = pd.read_table(path_1 + "deseq2_results_modified.tsv")
    roary_df = pd.read_csv(path_in + "modified_roary_output.csv")
    #roary_df = roary_df.fillna("None")
    roary_dict = {}
    '''for index, row in roary_df.iterrows():
        row_info = [i for i in row]
        if "None" not in row_info[1::]:
            roary_dict.update({row_info[0]:row_info[1::]})'''
    #ortholog_column = seq_df.loc[:,"ortholog_group"]
    roary_df = roary_df.dropna()
    condition_list = seq_df["coef"].unique()
    strain_list = [i for i in seq_df["strain"].unique()]
    strain_list.remove("Hybrid")
    PCA_df = pd.DataFrame(roary_df.loc[:,"Gene"])
    PCA_dict = {}
    print(roary_df.shape[0])
    for strain in strain_list:
        strain_df = seq_df[seq_df.loc[:,"strain"] == strain]
        for con in condition_list:
            data_list = []
            
            con_strain_df = strain_df[strain_df.loc[:,"coef"] == con]
            temp_con = con.split("_")[1]
            column_l = strain + "_" + temp_con
            data = con_strain_df.loc[:,"log2FoldChange"]
            each_dict = dict(zip([i for i in con_strain_df.loc[:,"id"]],[i for i in data]))
            #print(con_strain_df.shape[0])
            for gene, ID in zip(roary_df["Gene"],roary_df[strain]):
                try:
                    fold_change = each_dict[ID]
                    data_list.append(fold_change)
                    
                    #print(fold_change)
                except:
                    data_list.append(None)
                    #print("Not Found in original file")
                    #continue
            PCA_df[column_l] = data_list
    PCA_df = PCA_df.dropna()
    PCA_df.to_csv(path_in + "PCA_table.csv", index = False)    
                
            


if __name__ == "__main__":
    path_1 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/PCA/"
    PCA_file(path_1)