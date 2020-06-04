import pandas as pd
import pickle as pk

#replace gene id in roary output to loci id

def change_id(path_in,path_out):
    raw_roary_df = pd.read_csv(path_in + "gene_presence_absence.csv")
    roary_df = raw_roary_df[["Gene","cor_Sinorhizobium_meliloti_1021.ASM696v1.46","cor_Sinorhizobium_meliloti_ak83.ASM14779v3.46","cor_Sinorhizobium_meliloti_bl225c.ASM14777v3.46"]]
    #roary_df = raw_roary_df[["Gene","cor_Sinorhizobium_meliloti_1021.ASM696v1.46","cor_Sinorhizobium_meliloti_2011.swap","cor_Sinorhizobium_meliloti_bl225c.ASM14777v3.46"]]
    loci_file = open(path_in + "loci_id.dat","rb")
    loci_dict = pk.load(loci_file)
    roary_df = roary_df.fillna("None")
    new_roary_df = pd.DataFrame()
    new_roary_df["Gene"] = roary_df["Gene"]
    for num in range(1,4):
        gene_list = roary_df.iloc[:,num]
        replace_id_list = []
        for g in gene_list:
            try:
                replace_id = loci_dict[g][1]
            except:
                if g == "None":
                    replace_id = "NA"
            replace_id_list.append(replace_id)
        temp_df = pd.DataFrame(replace_id_list,columns = [roary_df.columns.values[num]])
        new_roary_df[roary_df.columns.values[num]] = temp_df
    new_roary_df.to_csv(path_out + "modified_roary_output.csv", index = False)
    #new_roary_df.to_csv(path_out + "modified_roary_output_with_Hybrid.csv", index = False)

            
            
            

if __name__ == "__main__":
    path_1 = 
    #path_1 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap/"
    path_2 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis/"
    change_id(path_1,path_2)