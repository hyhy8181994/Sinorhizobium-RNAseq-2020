import os
import pandas as pd
import pickle as pk

# change gene id in deseq2 ouput to loci id

def swap(reference_file_path,target_file):
    file_list = os.listdir(reference_file_path)
    gene_info_dict = {}
    for reference_file in file_list:
        with open(reference_file_path+reference_file,"r") as r_f:
            for r_f_line in r_f:
                if "#" not in r_f_line:
                    gene_row = r_f_line.split("\t")
                    chr = gene_row[0]
                    gene_info = gene_row[8].split(";")
                    gene_id = gene_info[0].strip("ID=")
                    loc_id = gene_info[1].strip("Parent=gene:")
                    gene_dict = {gene_id: [chr,loc_id]}
                    try:
                        gene_info_dict.update(gene_dict)
                    except:
                        pass
                if "FASTA" in r_f_line:
                    break
    in_df = pd.read_table(target_file)
    id_list = in_df["id"]
    replace_data_list = []
    for id in id_list:
        try:
            replace_data = gene_info_dict[id]
        except:
            replace_data = ["Chromosome","SMc05009"]
        replace_data_list.append(replace_data)
    replace_data_df = pd.DataFrame(replace_data_list, columns = ["replicon","id"])
    in_df["id"] = replace_data_df["id"].values
    in_df.insert(3,"replicon",replace_data_df["replicon"].values)
    in_df.to_csv(os.path.split(target_file)[0] + "/deseq2_results_modified.tsv",index = False,na_rep='NA',sep = "\t")
    with open(os.path.split(target_file)[0] + "/loci_id.dat","wb") as data:
        pk.dump(gene_info_dict,data)


            


if __name__ == "__main__":
    path = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/gff_corrected/"
    in_file = '/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/deseq2_results.tsv'
    swap(path,in_file)