import pandas as pd


#generate some stats for significant expressed genes

# total; total significant; significant in each strain and each condition; percentage
def stats(input_file,out_stats):
    df_sig = pd.read_csv(input_file + "significant_output.csv")
    df_main = pd.read_table(input_file + "deseq2_results_modified.tsv")
    conditions = df_main["coef"].unique()
    stats_dict = {}
    strain_list = df_main["strain"].unique()
    for strain in strain_list:
        stats_list = []
        for con in conditions:
            each_strain_df = df_main[df_main["strain"] == strain]
            each_strain_df = each_strain_df[each_strain_df["coef"] == con]
            total_number = each_strain_df.shape[0]
            sig_strain_df = df_sig[df_sig["strain"] == strain]
            sig_strain_df = sig_strain_df[sig_strain_df["coef"] == con]
            sig_number = sig_strain_df.shape[0]
            each_stats = "{}({}%)".format(sig_number,round(sig_number*100/total_number,2))
            stats_list.append(each_stats)
        temp_sig_dict = {strain:stats_list}
        stats_dict.update(temp_sig_dict)
    
    stats_df = pd.DataFrame.from_dict(stats_dict,orient="index",columns=conditions)
    stats_df.to_csv(out_stats + "significant_stats.csv")

        

if __name__ == "__main__":
    path = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/"  
    stats(path,path)