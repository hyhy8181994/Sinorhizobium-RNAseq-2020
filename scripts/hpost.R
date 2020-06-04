library(hpost)

folder_path = "/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/annotation/output/hpost/"
table_list = list.files(folder_path,pattern = ".csv")
setwd(folder_path)
# 1 upregulate; 2 downregulate; 0 non-significant
for(table_name in table_list){
  input_table = read.table(table_name, sep = ",", header =TRUE, na.strings = "NA")

  logfold = as.numeric(input_table[,2])
  input_matrix = matrix(logfold,ncol = 1)
  cat = as.character("exp")
  grp = as.numeric(input_table[,6])
  KEGG = as.character(input_table[,4])
  COG = as.character(input_table[,5])
  data = buildDataSet(logfold,fct = cat, grp = grp)
  
  output_KEGG = hpostSingle(data, KEGG)
  output_COG = hpostSingle(data, COG)
  
  output_KEGG_up <- output_KEGG[!(output_KEGG$grp==1),]
  output_KEGG_up <- output_KEGG_up[!(output_KEGG_up$adj.p.value>0.05),]
  output_KEGG_down <- output_KEGG[!(output_KEGG$grp==2),]
  output_KEGG_down <- output_KEGG_down[!(output_KEGG_down$adj.p.value>0.05),]
  
  output_COG_up <- output_COG[!(output_COG$grp==1),]
  output_COG_up <- output_COG_up[!(output_COG_up$adj.p.value>0.05),]
  output_COG_down <- output_COG[!(output_COG$grp==2),]
  output_COG_down <- output_COG_down[!(output_COG_down$adj.p.value>0.05),]
  name = strsplit(table_name, "//.")[1]
  
  write.csv(output_KEGG_up,paste("output_KEGG_up_",name,sep = ""))
  write.csv(output_KEGG_down,paste("output_KEGG_down_",name,sep = ""))
  write.csv(output_COG_up,paste("output_COG_up_",name,sep = ""))
  write.csv(output_COG_down,paste("output_COG_down_",name,sep = ""))
  
}

