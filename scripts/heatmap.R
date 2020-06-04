#source("https://bioconductor.org/biocLite.R")
#BiocManager::install("ComplexHeatmap")
#library("BiocManager")
library("ComplexHeatmap")
library(cluster)
library("heatmaply")
library("grDevices")

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
  
  figure_name = paste(figure_name_1[1],".pdf", sep = "")
  
  pdf(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap",figure_name,sep = "/"),family = "Helvetica", width = 10, height = 6)
  #png(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap",figure_name,sep = "/"), width = 1500, height = 1000)
  

  hm = Heatmap(df_matrix,
          name = "log2 Fold Change",
          clustering_method_columns = "complete",
          row_title = "Cultivar",
          column_title = strain_name,
          show_column_names = FALSE,
          row_names_gp = gpar(fontsize = 15),
          column_title_gp = gpar(fontsize = 15, fontface = "bold"),
          row_title_gp = gpar(fontsize = 15, fontface = "bold"),
          
          column_names_gp = gpar(fontsize = 12)
          #width = unit(170, "inch"), height = unit(90, "inch")
          
  )
  
  draw(hm)

  dev.off()
}

#core_gene_heatmap
####################################

setwd("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis")

condition_list = c("ortholog_heatmap_camporegio.csv","ortholog_heatmap_lodi.csv","ortholog_heatmap_luteolin.csv","ortholog_heatmap_verbena.csv")
for (name in condition_list){

df_camporeglo = read.csv(name)

row_name = df_camporeglo$X

df_camporeglo$X = NULL

matrix_camproeglo = as.matrix(df_camporeglo)

rownames(matrix_camproeglo) = row_name


plot_name_1 = unlist(strsplit(name,"\\."))

plot_name = paste(plot_name_1[1],".pdf", sep = "")



pdf(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis",plot_name,sep = "/"),family = "Helvetica", width = 10, height = 6)
#png(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/roary_analysis",plot_name,sep = "/"), width = 1500, height = 1000)


column_name = unlist(strsplit(plot_name_1,"\\_"))

ht = Heatmap(matrix_camproeglo,
        name = "log2 Fold Change",
        column_title = column_name[3],
        clustering_method_rows = "complete",
        show_column_names = FALSE,
        row_title = "Strain",
        row_names_gp = gpar(fontsize = 15),
        column_title_gp = gpar(fontsize = 15, fontface = "bold"),
        row_title_gp = gpar(fontsize = 15, fontface = "bold"),
        column_names_gp = gpar(fontsize = 12)
        #width = unit(20, "inch"), height = unit(8, "inch")
        )
draw(ht)
dev.off()
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
  
  plot_name = paste(plot_name_1[1],".pdf", sep = "")
  
  
  pdf(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap/Heatmap",plot_name,sep = "/"), family = "Helvetica", height = 130, width = 250)
  #png(paste("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/heatmap/for_heatmap",plot_name,sep = "/"), width = 1500, height = 1000)
  
  
  column_name = unlist(strsplit(plot_name_1,"\\_"))
  
  ht = Heatmap(matrix_camproeglo,
               name = "log2 Fold Change",
               column_title = column_name[3],
               clustering_method_rows = "complete",
               show_column_names = FALSE,
               row_title = "Strain",
               row_names_gp = gpar(fontsize = 70),
               column_title_gp = gpar(fontsize = 70, fontface = "bold"),
               row_title_gp = gpar(fontsize = 70, fontface = "bold"),
               
               column_names_gp = gpar(fontsize = 12)
  )
  draw(ht)
  dev.off()

}
        