import pandas as pd
import os


#change gene names in each core heatmap table to 1021 gene name

def rename(path_in,path_read):
    roary_df = pd.read_csv(path_in + "modified_roary_output.csv")
    roary_df = roary_df.fillna("None")
    roary_dict = {}
    for index, row in roary_df.iterrows():
        #print(row[0]),print(row[1])
        roary_dict.update({row[0]:row[1]})
    for path in os.listdir(path_read):
        file_path = os.path.join(path_read,path)
        core_csv = pd.read_csv(file_path,index_col= 0)
        old_column_name = core_csv.columns
        new_column_name = []
        for col in old_column_name:
            if roary_dict[col] != "None":
                new_column_name.append(roary_dict[col])
            else:
                print(col)
                new_column_name.append(col)
        core_csv.columns = new_column_name
        core_csv.to_csv(path_in + path)

if __name__ == "__main__":
    path_1 = "~/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis/"
    path_2 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis/core_heatmap"
    rename(path_1,path_2)