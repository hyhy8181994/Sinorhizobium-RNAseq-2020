import pandas as pd
import os
import pickle as pk

#generate overall summary table

def replicon_table(path_in):
    sig_df = pd.read_csv(path_in + "significant_output.csv")
    task_list = ["overall"]
    condition_list = sig_df["coef"].unique()
    strain_list = ["1021","BL225C","AK83"]
    orth_df = pd.read_csv(path_in + "modified_roary_output.csv")
    orth_df = orth_df.fillna("None")
    orth_dict = {}
    pk_file = open(path_in + "roary_dict.dat","rb")
    roary_dict = pk.load(pk_file)
    catergory_list = []
    ortholog_group_list = []
    for index, row in orth_df.iterrows():
        if "None" in list(row):
            orth_dict.update({row[0]:"accessory"})
        else:
            orth_dict.update({row[0]:"core"})
    with open(path_in + "category_dict.dat", "wb") as cf:
        pk.dump(orth_dict,cf)
    for id in sig_df["id"]:
        try:
            ortholog_group = roary_dict[id]
            ortholog_group_list.append(ortholog_group)
            gene_catergory = orth_dict[ortholog_group]
            catergory_list.append(gene_catergory)
        except:
            catergory_list.append("None")
            ortholog_group_list.append("None")
    sig_df["category"] = catergory_list
    sig_df["ortholog_group"] = ortholog_group_list
    
    for i in condition_list:
        task_list.append(i.split("_")[1])
    for task in task_list:
        open(path_in + "{}_summary_table.csv".format(task),"w+")
        with open(path_in + "{}_summary_table.csv".format(task),"a") as table:
            table.write(task + ",,,\n")
            table.write("Category,Rm1021,BL225C,AK83\n")
            
            if task == "overall":
                overall_whole_genome = []
                overall_genome_up = []
                overall_genome_down = []
                overall_chro = []
                overall_chro_up = []
                overall_chro_down = []
                each_chro = []
                each_chro_up = []
                each_chro_down = []
                p = []
                p_up = []
                p_down = []
                core = []
                core_up = []
                core_down = []
                acc = []
                acc_up = []
                acc_down = []
                
                for strain in strain_list:    
                    overall_strain_df = sig_df[sig_df.loc[:,"strain"] == strain]
                    overall_whole_genome.append(overall_strain_df.shape[0])
                    up_overall = overall_strain_df[overall_strain_df.loc[:,"log2FoldChange"] > 0]
                    down_overall = overall_strain_df[overall_strain_df.loc[:,"log2FoldChange"] < 0]
                    if strain != "AK83":
                        overall_chr_df = overall_strain_df[overall_strain_df.loc[:,"replicon"] == "Chromosome"]
                        overall_chro.append(overall_chr_df.shape[0])
                        overall_chro_up.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] > 0].shape[0])
                        overall_chro_down.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] < 0].shape[0])
                        overall_p_df = overall_strain_df[overall_strain_df["replicon"].str.contains('p')]
                        for pl in overall_p_df["replicon"].unique():
                            each_p_df = overall_p_df[overall_p_df["replicon"] == pl]
                            p.append(each_p_df.shape[0])
                            p_up.append(each_p_df[each_p_df["log2FoldChange"] > 0].shape[0])
                            p_down.append(each_p_df[each_p_df["log2FoldChange"] < 0].shape[0])
                        core_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "core"]
                        acc_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "accessory"]
                        #print(core_df)
                        core_df_up = core_df[core_df.loc[:,"log2FoldChange"] > 0]
                        core_df_down = core_df[core_df.loc[:,"log2FoldChange"] < 0]
                        acc_df_up = acc_df[acc_df.loc[:,"log2FoldChange"] > 0]
                        acc_df_down = acc_df[acc_df.loc[:,"log2FoldChange"] < 0]
                        core.append(core_df.shape[0])
                        core_up.append(core_df_up.shape[0])
                        core_down.append(core_df_down.shape[0])
                        acc.append(acc_df.shape[0])
                        acc_up.append(acc_df_up.shape[0])
                        acc_down.append(acc_df_down.shape[0])
                        
                        
                    else:
                        all_count = 0
                        up_count = 0
                        down_count = 0
                        plasmid_list = []
                        plasmid_up = []
                        plasmid_down = []
                        for i in range(1,4):
                            overall_chr_df = overall_strain_df[overall_strain_df.loc[:,"replicon"] == str(i)]
                            each_chro_up.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] > 0].shape[0])
                            each_chro_down.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] < 0].shape[0])
                            core_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "core"]
                            acc_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "accessory"]
                            plasmid_df = overall_strain_df[overall_strain_df.loc[:,"replicon"].str.contains("pSINME")]
                            plasmid_list.append(plasmid_df.shape[0])
                            plasmid_up.append(plasmid_df[plasmid_df.loc[:,"log2FoldChange"] > 0].shape[0])
                            plasmid_down.append(plasmid_df[plasmid_df.loc[:,"log2FoldChange"] < 0].shape[0])

                            core_df_up = core_df[core_df.loc[:,"log2FoldChange"] > 0]
                            core_df_down = core_df[core_df.loc[:,"log2FoldChange"] < 0]
                            acc_df_up = acc_df[acc_df.loc[:,"log2FoldChange"] > 0]
                            acc_df_down = acc_df[acc_df.loc[:,"log2FoldChange"] < 0]
                            core.append(core_df.shape[0])
                            core_up.append(core_df_up.shape[0])
                            core_down.append(core_df_down.shape[1])
                            acc.append(acc_df.shape[0])
                            acc_up.append(acc_df_up.shape[0])
                            acc_down.append(acc_df_down.shape[0])
                            each_chro.append(overall_chr_df.shape[0])
                            all_count += overall_chr_df.shape[0]
                            up_count += overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] > 0].shape[0]
                            down_count += overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] < 0].shape[0]

                        overall_chro.append(all_count)
                        overall_chro_up.append(up_count)
                        overall_chro_down.append(down_count)
                    overall_genome_up.append(up_overall.shape[0])
                    overall_genome_down.append(down_overall.shape[0])
                table.write("DEGs - whole genome,{},{},{}\n".format(overall_whole_genome[0],overall_whole_genome[1],overall_whole_genome[2]))
                table.write("Upregulated,{},{},{}\n".format(overall_genome_up[0],overall_genome_up[1],overall_genome_up[2]))
                table.write("Downregulatd,{},{},{}\n".format(overall_genome_down[0],overall_genome_down[1],overall_genome_down[2]))
                table.write("Chromosome,{},{},{}\n".format(overall_chro[0],overall_chro[1],each_chro[0]))
                table.write("Upregulatd,{},{},{}\n".format(overall_chro_up[0],overall_chro_up[1],each_chro_up[0]))
                table.write("Downregulatd,{},{},{}\n".format(overall_chro_down[0],overall_chro_down[1],each_chro_down[0]))
                table.write("pSymB,{},{},{}\n".format(p[0],p[2],each_chro[1]))
                table.write("Upregulatd,{},{},{}\n".format(p_up[0],p_up[2],each_chro_up[1]))
                table.write("Downregulatd,{},{},{}\n".format(p_down[0],p_down[2],each_chro_down[1]))
                table.write("pSymA,{},{},{}\n".format(p[1],p[3],each_chro[2]))
                table.write("Upregulatd,{},{},{}\n".format(p_up[1],p_up[3],each_chro_up[2]))
                table.write("Downregulatd,{},{},{}\n".format(p_down[1],p_down[3],each_chro_down[2]))
                table.write("Plasmid,NA,NA,{}\n".format(plasmid_list[0]))
                table.write("Upregulatd,NA,NA,{}\n".format(plasmid_up[0]))
                table.write("Downregulatd,NA,NA,{}\n".format(plasmid_down[0]))
                table.write("DEGs - core genome,{},{},{}\n".format(core[0],core[1],core[2]))
                table.write("Upregulated,{},{},{}\n".format(core_up[0],core_up[1],core_up[2]))
                table.write("Downregulated,{},{},{}\n".format(core_down[0],core_down[1],core_down[2]))
                table.write("DEGs - accessory genome,{},{},{}\n".format(acc[0],acc[1],acc[2]))
                table.write("Upregulated,{},{},{}\n".format(acc_up[0],acc_up[1],acc_up[2]))
                table.write("Downregulated,{},{},{}\n".format(acc_down[0],acc_down[1],acc_down[2]))
            else:
                overall_whole_genome = []
                overall_genome_up = []
                overall_genome_down = []
                overall_chro = []
                overall_chro_up = []
                overall_chro_down = []
                each_chro = []
                each_chro_up = []
                each_chro_down = []
                p = []
                p_up = []
                p_down = []
                core = []
                core_up = []
                core_down = []
                acc = []
                acc_up = []
                acc_down = []
                
                for strain in strain_list:    
                    for con in condition_list:
                        if task in con:
                            condition = con
                    condition_df = sig_df[sig_df.loc[:,"coef"] == condition]
                    overall_strain_df = condition_df[condition_df.loc[:,"strain"] == strain]
                    
                    overall_whole_genome.append(overall_strain_df.shape[0])
                    up_overall = overall_strain_df[overall_strain_df.loc[:,"log2FoldChange"] > 0]
                    down_overall = overall_strain_df[overall_strain_df.loc[:,"log2FoldChange"] < 0]
                    if strain != "AK83":
                        overall_chr_df = overall_strain_df[overall_strain_df.loc[:,"replicon"] == "Chromosome"]
                        overall_chro.append(overall_chr_df.shape[0])
                        overall_chro_up.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] > 0].shape[0])
                        overall_chro_down.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] < 0].shape[0])
                        overall_p_df = overall_strain_df[overall_strain_df["replicon"].str.contains('p')]
                    
                        for pl in overall_p_df["replicon"].unique():
                            each_p_df = overall_p_df[overall_p_df["replicon"] == pl]
                           
                            
                            p.append(each_p_df.shape[0])
                            p_up.append(each_p_df[each_p_df["log2FoldChange"] > 0].shape[0])
                            p_down.append(each_p_df[each_p_df["log2FoldChange"] < 0].shape[0])
                        core_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "core"]
                        acc_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "accessory"]
                        #print(core_df)
                        core_df_up = core_df[core_df.loc[:,"log2FoldChange"] > 0]
                        core_df_down = core_df[core_df.loc[:,"log2FoldChange"] < 0]
                        acc_df_up = acc_df[acc_df.loc[:,"log2FoldChange"] > 0]
                        acc_df_down = acc_df[acc_df.loc[:,"log2FoldChange"] < 0]
                        core.append(core_df.shape[0])
                        core_up.append(core_df_up.shape[0])
                        core_down.append(core_df_down.shape[0])
                        acc.append(acc_df.shape[0])
                        acc_up.append(acc_df_up.shape[0])
                        acc_down.append(acc_df_down.shape[0])
                    
                        
                    else:
                        all_count = 0
                        up_count = 0
                        down_count = 0
                        plasmid_list = []
                        plasmid_up = []
                        plasmid_down = []
    
                        for i in range(1,4):
                            overall_chr_df = overall_strain_df[overall_strain_df.loc[:,"replicon"] == str(i)]
                            each_chro_up.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] > 0].shape[0])
                            each_chro_down.append(overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] < 0].shape[0])
                            core_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "core"]
                            acc_df = overall_strain_df[overall_strain_df.loc[:,"category"] == "accessory"]
                            plasmid_df = overall_strain_df[overall_strain_df.loc[:,"replicon"].str.contains("pSINME")]
                            plasmid_list.append(plasmid_df.shape[0])
                            plasmid_up.append(plasmid_df[plasmid_df.loc[:,"log2FoldChange"] > 0].shape[0])
                            plasmid_down.append(plasmid_df[plasmid_df.loc[:,"log2FoldChange"] < 0].shape[0])

                            core_df_up = core_df[core_df.loc[:,"log2FoldChange"] > 0]
                            core_df_down = core_df[core_df.loc[:,"log2FoldChange"] < 0]
                            acc_df_up = acc_df[acc_df.loc[:,"log2FoldChange"] > 0]
                            acc_df_down = acc_df[acc_df.loc[:,"log2FoldChange"] < 0]
                            core.append(core_df.shape[0])
                            core_up.append(core_df_up.shape[0])
                            core_down.append(core_df_down.shape[0])
                            acc.append(acc_df.shape[0])
                            acc_up.append(acc_df_up.shape[0])
                            acc_down.append(acc_df_down.shape[0])
                            each_chro.append(overall_chr_df.shape[0])
                            all_count += overall_chr_df.shape[0]
                            up_count += overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] > 0].shape[0]
                            down_count += overall_chr_df[overall_chr_df.loc[:,"log2FoldChange"] < 0].shape[0]

                        overall_chro.append(all_count)
                        overall_chro_up.append(up_count)
                        overall_chro_down.append(down_count)
                    overall_genome_up.append(up_overall.shape[0])
                    overall_genome_down.append(down_overall.shape[0])
                '''if task == "lodi":
                    p = [0,0,1,1]
                    each_chro = [2,2,1]
                    each_chro_up = [1,2,0]
                    each_chro_down = [1,0,1]
                    p_up = [0,0,0,0]
                    p_down = [0,0,1,1]
                    print(core_df)'''
                   
                table.write("DEGs - whole genome,{},{},{}\n".format(overall_whole_genome[0],overall_whole_genome[1],overall_whole_genome[2]))
                table.write("Upregulated,{},{},{}\n".format(overall_genome_up[0],overall_genome_up[1],overall_genome_up[2]))
                table.write("Downregulatd,{},{},{}\n".format(overall_genome_down[0],overall_genome_down[1],overall_genome_down[2]))
                table.write("Chromosome,{},{},{}\n".format(overall_chro[0],overall_chro[1],each_chro[0]))
                table.write("Upregulatd,{},{},{}\n".format(overall_chro_up[0],overall_chro_up[1],each_chro_up[0]))
                table.write("Downregulatd,{},{},{}\n".format(overall_chro_down[0],overall_chro_down[1],each_chro_down[0]))
                table.write("pSymB,{},{},{}\n".format(p[0],p[2],each_chro[1]))
                table.write("Upregulatd,{},{},{}\n".format(p_up[0],p_up[2],each_chro_up[1]))
                table.write("Downregulatd,{},{},{}\n".format(p_down[0],p_down[2],each_chro_down[1]))
                table.write("pSymA,{},{},{}\n".format(p[1],p[3],each_chro[2]))
                table.write("Upregulatd,{},{},{}\n".format(p_up[1],p_up[3],each_chro_up[2]))
                table.write("Downregulatd,{},{},{}\n".format(p_down[1],p_down[3],each_chro_down[2]))
                table.write("Plasmid,NA,NA,{}\n".format(plasmid_list[0]))
                table.write("Upregulatd,NA,NA,{}\n".format(plasmid_up[0]))
                table.write("Downregulatd,NA,NA,{}\n".format(plasmid_down[0]))
                table.write("DEGs - core genome,{},{},{}\n".format(core[0],core[1],core[2]))
                table.write("Upregulated,{},{},{}\n".format(core_up[0],core_up[1],core_up[2]))
                table.write("Downregulated,{},{},{}\n".format(core_down[0],core_down[1],core_down[2]))
                table.write("DEGs - accessory genome,{},{},{}\n".format(acc[0],acc[1],acc[2]))
                table.write("Upregulated,{},{},{}\n".format(acc_up[0],acc_up[1],acc_up[2]))
                table.write("Downregulated,{},{},{}\n".format(acc_down[0],acc_down[1],acc_down[2]))
    sig_df.to_csv(path_in + "significant_output_with_category.csv", index = False)

        



if __name__ == "__main__":
    path_1 = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/strain_sig_data/"
    replicon_table(path_1)