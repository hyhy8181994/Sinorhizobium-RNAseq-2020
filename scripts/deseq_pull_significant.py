import pandas as pd


# extract significant differentiate expressed gene from deseq2 ouput


data = pd.read_table("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/deseq2_results_modified.tsv")
out_file = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/significant_output.csv","w")
out_file.write("strain,coef,id,replicon,baseMean,log2FoldChange,lfcSE,pvalue,padj\n")
out_file.close()
for index, row in data.iterrows():
    if row["log2FoldChange"] > 1 or row["log2FoldChange"] < -1:
        if row["padj"] < 0.01:
            out_file = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/significant_output.csv","a")
            out_file.write("{},{},{},{},{},{},{},{},{}\n".format(row["strain"],row['coef'],row['id'],row['replicon'],row['baseMean'],row['log2FoldChange'],row['lfcSE'],row['pvalue'],row['padj']))
            out_file.close()

sig_data = pd.read_csv("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/significant_output.csv")          
data_1021 = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/1021_data.csv","w")
data_1021.write("id,strain,coef,log2FoldChange\n")
data_AK83 = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/AK83_data.csv","w")
data_AK83.write("id,strain,coef,log2FoldChange\n")
data_BL225 = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/BL225_data.csv","w")
data_BL225.write("id,strain,coef,log2FoldChange\n")
data_hybrid = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/Hybrid_data.csv","w")
data_hybrid.write("id,strain,coef,log2FoldChange\n")
data_1021.close()
data_AK83.close()
data_BL225.close()
data_hybrid.close()


for index, row in sig_data.iterrows():
    if row["strain"] == "1021":
        out_file_1021 = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/1021_data.csv","a")
        out_file_1021.write("{},{},{},{}\n".format(row['id'],row["strain"],row['coef'],row['log2FoldChange']))
        out_file_1021.close()
    elif row["strain"] == "AK83":
        out_file_AK83 = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/AK83_data.csv","a")
        out_file_AK83.write("{},{},{},{}\n".format(row['id'],row["strain"],row['coef'],row['log2FoldChange']))
        out_file_AK83.close()
    elif row["strain"] == "BL225C":
        out_file_BL225 = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/BL225_data.csv","a")
        out_file_BL225.write("{},{},{},{}\n".format(row['id'],row["strain"],row['coef'],row['log2FoldChange']))
        out_file_BL225.close()
    elif row["strain"] == "Hybrid":
        out_file_Hybrid = open("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/Hybrid_data.csv","a")
        out_file_Hybrid.write("{},{},{},{}\n".format(row['id'],row["strain"],row['coef'],row['log2FoldChange']))
        out_file_Hybrid.close()





