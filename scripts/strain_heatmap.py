import pandas as pd
import os

#generate strain specific heatmap table

def move(main,input_path,out_path):
    df = pd.read_csv(input_path)
    #out_df = pd.DataFrame(columns= ["id","camporegio","lodi","luteolin","verbena"])
    out_dict = dict()
    for index, row in df.iterrows():
        main_df = open(main,"r")
        fold_dict = dict()
        condition_list = ["condition_camporegio_vs_blank","condition_lodi_vs_blank","condition_luteolin_vs_blank","condition_verbena_vs_blank"]
        df_row_id = row["id"]
        #print(df_row_id)
        df_row_coef = row["coef"]
        df_row_fold = {df_row_coef:row["log2FoldChange"]}
        condition_list.remove(df_row_coef)
        for line in main_df:
            main_row = line.strip("\n").split("\t")
            if df_row_id in main_row[2]:
                if main_row[1] in condition_list:
                    df_row_fold.update({main_row[1]:main_row[5]})
                    fold_dict.update(df_row_fold)
                    condition_list.remove(main_row[1])
        try:
            out_dict.update({df_row_id:fold_dict}) 
        except:
            pass
    #print(out_dict)
    df = pd.DataFrame.from_dict(out_dict)
    out_file_name = os.path.basename(input_path).split(".")[0] + "_heatmap.csv"
    df.to_csv(out_path + out_file_name)


        
        

            
      


                



path = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/"
main_df = path + "deseq2_results_modified.tsv"

for i in ["1021","AK83","BL225","Hybrid"]:
    move(main_df,path + "{}_data.csv".format(i),path)
