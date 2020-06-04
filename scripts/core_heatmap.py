import pandas as pd 
import pickle as pk


#generate core heatmap table

def ortholog_heatmap_data(path_roary,path_sig):
    roary_df = pd.read_csv(path_roary + "modified_roary_output.csv")#
    sig_df = pd.read_csv(path_sig + "significant_output.csv")
    seq_df = pd.read_table(path_sig + "deseq2_results_modified.tsv")
    sig_df = sig_df.fillna("None")
    conditions = sig_df["coef"].unique()
    roary_df = roary_df.fillna("None")
    roary_dict = {}
    ortho_dict = {}
    for index, row in roary_df.iterrows():
        row_info = [i for i in row]
        roary_dict.update({row_info[0]:row_info[1::]})
    for num in range(1,4):
        for id,gene in zip(roary_df.iloc[:,num],roary_df.iloc[:,0]):
            if id == "None":
                pass
            else:
                ortho_dict.update({id:gene})

    for con in conditions:
        each_con_df = sig_df[sig_df["coef"] == con]
        #each_con_df = each_con_df[each_con_df["strain"] != "AK83"]
        sub_sig_df = sig_df[sig_df["coef"] == con]
        sub_seq_df = seq_df[seq_df["coef"] == con]
        #sub_seq_df = sub_seq_df[sub_seq_df["strain"] != "AK83"]
        con_dict = {}
        for id in each_con_df["id"]:
            try:
                ortho_list = roary_dict[ortho_dict[id]]
            except:
                pass
            value_list = []
            for number in range(3):
                #ortho_row = sub_sig_df[sub_sig_df["id"] == ortho]
                #print(ortho_row)
                #value = [i for i in ortho_row["log2FoldChange"]]
                ortho = ortho_list[number]
                if number == 0:
                    strain = "1021"
                elif number == 1:
                    strain = "AK83" #
                elif number == 2:
                    strain = "BL225C"
                sub_strain_df = sub_seq_df[sub_seq_df["strain"] == strain]
                
                value = sub_strain_df[sub_strain_df["id"] == ortho]["log2FoldChange"]
                if len(value) == 0:    
                    value = ["NA"]
                    
                value_list.append([i for i in value])
                
                #print(value_list)
                '''
                if len(value) == 0:
                    #value = [i for i in sub_seq_df[sub_seq_df["id"] == ortho]["log2FoldChange"]]
                    if len(value) == 0:
                        value = ["NA"]
                    elif len(value) == 2:
                        value = [value[0]]
                if len(value) == 2:
                    value = [value[0]]
                value_list.append(value)
                #count += 1'''
            try:
                three_values = [j for i in value_list for j in i]
                
                #if "NA" not in three_values:
                if three_values.count("NA") < 1:
                    con_dict.update({ortho_dict[id]:[j for i in value_list for j in i]})
            except:
                pass
        with open(path_sig + "roary_dict.dat", "wb") as od:
            pk.dump(ortho_dict,od)
        #print(con_dict)
        con_df = pd.DataFrame.from_dict(con_dict)
        #change column name
        con_df.rename({0:"1021",1:"AK83",2:"BL225"}, inplace = True)
        
        print(con_df.shape[1])
        con_df = con_df.loc[:,~con_df.columns.duplicated()] 
        print(con_df.shape[1])
        con_df.to_csv(path_roary + "ortholog_heatmap_{}.csv".format(con.split("_")[1]))
            #table_dict = {ortho_dict[gene]:sig_df[sig_df["id"]]}

        

if __name__ == "__main__":
    #path_1 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap/"
    path_1 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis/"

    

    path_2 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/"
    ortholog_heatmap_data(path_1,path_2)