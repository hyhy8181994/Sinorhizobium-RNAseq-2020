library(heatmaply)

#strain heatmap
###################################
setwd("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap")

file_list = c("1021_data_heatmap.csv","AK83_data_heatmap.csv","BL225_data_heatmap.csv","Hybrid_data_heatmap.csv")


for (file in file_list){  
  df = read.csv(file, header = TRUE)
  
  df_row_name = c("camporegio","lodi","luteolin","verbena")
  df$X = NULL
  
  df_matrix = as.matrix(df)
  rownames(df_matrix) = df_row_name
  
  figure_name_1 = unlist(strsplit(file,"\\."))
  
  strain_name = unlist(strsplit(figure_name_1[1],"\\_"))[1]
  
  if (strain_name == "1021"){
    strain_name = "Rm1021"
  }
  
  figure_name = paste(figure_name_1[1],".html", sep = "")
  
  heatmaply(df_matrix,hclustfun_col = hclust ,hclust_method = "complete", file = figure_name, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"),xlab = strain_name, ylab = "Cultivar")
}
 
setwd("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis")

condition_list = c("ortholog_heatmap_camporegio.csv","ortholog_heatmap_lodi.csv","ortholog_heatmap_luteolin.csv","ortholog_heatmap_verbena.csv")

 
for (name in condition_list){
  
  df_camporeglo = read.csv(name)
  
  row_name = df_camporeglo$X
  
  df_camporeglo$X = NULL
  
  matrix_camproeglo = as.matrix(df_camporeglo)
  
  rownames(matrix_camproeglo) = row_name
  
  
  plot_name_1 = unlist(strsplit(name,"\\."))
  
  plot_name = paste(plot_name_1[1],".html", sep = "")
  
  column_name = unlist(strsplit(plot_name_1,"\\_"))
  
  heatmaply(matrix_camproeglo, hclust_method = "complete", file = plot_name, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"),xlab = column_name[3], ylab = "Strain")
  
  
}

#################Hybrid core


setwd("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap/Heatmap")

condition_list = c("ortholog_heatmap_camporegio.csv","ortholog_heatmap_lodi.csv","ortholog_heatmap_luteolin.csv","ortholog_heatmap_verbena.csv")
for (name in condition_list){
  
  df_camporeglo = read.csv(name)
  
  row_name = df_camporeglo$X
  
  df_camporeglo$X = NULL
  
  matrix_camproeglo = as.matrix(df_camporeglo)
  
  rownames(matrix_camproeglo) = row_name
  
  
  plot_name_1 = unlist(strsplit(name,"\\."))
  
  plot_name = paste(plot_name_1[1],".html", sep = "")
  
  
  #pdf(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap/Heatmap",plot_name,sep = "/"), family = "Helvetica", height = 130, width = 250)
  #png(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap",plot_name,sep = "/"), width = 1500, height = 1000)
  
  
  column_name = unlist(strsplit(plot_name_1,"\\_"))
  
  heatmaply(matrix_camproeglo, hclust_method = "complete", file = plot_name, scale_fill_gradient_fun = ggplot2::scale_fill_gradient2(low = "blue", high = "red"),xlab = column_name[3], ylab = "Strain")
  
  
}
