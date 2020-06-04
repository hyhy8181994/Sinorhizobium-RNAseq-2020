import pandas as pd
import pickle as pk
import os
import numpy as np

#generate table for hpost analysis

def occurance(df, total_df):
    module_list = []
    #df = df[df["KEGG_Module"] != "None"]
    #df = df[df["strain"] != "Hybrid"]
    for module in df["KEGG_Module"]:
        try:
            temp_list = module.split(",")
            for t in temp_list:
                module_list.append(t)
        except:
            module_list.append(module)
    count_df = pd.DataFrame(module_list, columns = ["KEGG"])
    #count_df.dropna("None")
    occurance_df = pd.DataFrame(count_df["KEGG"].value_counts())
    #gene_count = df["id"].shape[0]
    if total_df.shape[0] != 0:
        temp_total_occurance_list = []
        gene_count = []
        for total_module in total_df["KEGG_Module"]:
            try:
                temp_list = total_module.split(",")
                for t in temp_list:
                    temp_total_occurance_list.append(t)
            except:
                temp_total_occurance_list.append(total_module)
        for m in occurance_df.index.tolist():
            t_count = temp_total_occurance_list.count(m)
            gene_count.append(t_count)
        occurance_df["Total_occurance_count"] = gene_count
    else:
        gene_count = df["id"].shape[0]
        occurance_df["Total_gene_count"] = list(np.repeat(gene_count,occurance_df.shape[0]))
    return(occurance_df)

def name_change(path_in):
    file_list = os.listdir(path_in + "eggnog_file/")
    print(file_list)
    loci_file = open(path_in + "loci_id.dat", "rb")
    loci_dict = pk.load(loci_file)
    KEGG_dict = {}
    COG_dict = {}
    for file in file_list:
        #print(file)
        
        df = pd.read_table(path_in + "eggnog_file/" +  file)
        module_dict = {}
        CO_dict = {}
        #print(df.columns)
        small_df = df[["query_name","KEGG_Module","COG Functional cat."]]
        replace_list = []
        for name in small_df["query_name"]:
            try:
                loci_id = loci_dict[name.split(".")[0]][1]
                replace_list.append(loci_id)
            except:
                Warning("ID not found")
                exit(1)
        small_df.loc[:,"query_name"] = replace_list
        small_df.columns = ["id","KEGG_Module","COG"]
        small_df = small_df.fillna("NA")
        small_df.to_csv(path_in + "output/strain/" + file.split("_")[0] + ".csv", index = False)
        for index, row in small_df.iterrows():
            #module_element = row[1]
            '''try:
                temp_module_list = module_element.split(",")
            except:
                temp_module_list = [row[1]]'''
            try:
                split_KEGG = row[1].split(",")
                input_KEGG = split_KEGG[0] # only takes the first KEGG module
            except:
                input_KEGG = row[1]
            module_dict.update({row[0]:input_KEGG})
            CO_dict.update({row[0]:row[2]})
        KEGG_dict.update(module_dict)
        COG_dict.update(CO_dict)
    sig_df = pd.read_table(path_in + "deseq2_results_modified.tsv")
    sig_df = sig_df.fillna("NA")
    module_list_f_df = []
    COG_list_f_df =[]
    for id in sig_df["id"]:
        try:
            temp_KEGG = KEGG_dict[id]
            temp_COG = COG_dict[id]
        except:
            temp_KEGG = "NA"
            temp_COG = "NA"
        module_list_f_df.append(temp_KEGG)
        COG_list_f_df.append(temp_COG)
    #strain-condition-up/down table part
    sig_df["KEGG_Module"] = module_list_f_df
    sig_df["COG"] = COG_list_f_df
    strain_list = ["1021","AK83","BL225C"]
    condition_list = sig_df["coef"].unique()
    #small_sig_df = sig_df[["KEGG_Module","COG","padj"]]
    for strain in strain_list:
        strain_df = sig_df[sig_df["strain"] == strain]
         #occurance table for fisher test
        for con in condition_list:
            case = strain + "_" + con.split("_")[1]
            con_df = strain_df[strain_df["coef"] == con]
            small_con_df = con_df[["id","log2FoldChange","padj","KEGG_Module","COG"]]
            regulation_list = []
            for l in small_con_df["log2FoldChange"]:
                if (l > 1 or l == 1):
                    regulation_list.append(1) #up->1;down->2;other->0
                elif (l < -1 or l == -1):
                    regulation_list.append(2)
                else:
                    regulation_list.append(0)
            small_con_df["exp"] = regulation_list
            small_con_df.to_csv(path_in + "output/hpost/" + case + ".csv",index = False)
        ''' 
            #up_df = con_df[con_df["log2FoldChange"] > 0]
            #down_df = con_df[con_df["log2FoldChange"] < 0]
            
            
            ############occurance count#############3
            up_df_small = up_df[["id","KEGG_Module"]]
            down_df_small = down_df[["id","KEGG_Module"]]
            up_df_small.to_csv(path_in + "output/con_strain/" + case + "_up.csv", index = False)
            down_df_small.to_csv(path_in + "output/con_strain/" + case + "_down.csv", index = False)
            occ_df_up = occurance(up_df,small_df)
            occ_df_down = occurance(down_df,small_df)
            occ_df_up.to_csv(path_in + "output/occurance/"+ case + "_up_occurance_table.csv")
            occ_df_down.to_csv(path_in + "output/occurance/"+ case + "_down_occurance_table.csv")
        occ_small = occurance(small_df,pd.DataFrame())
        occ_small.to_csv(path_in + "output/occurance/" + file.split("_")[0] + "occurance_table.csv")'''
    
            

    

                


if __name__ == "__main__":
    path = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/annotation/"
    name_change(path)