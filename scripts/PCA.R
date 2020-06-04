library(ggplot2)

setwd("/home/rhuang06/Documents/Rhizobia/RNAseq_data/output_stats/PCA")

PCA_table = read.csv("PCA_table.csv", header = TRUE)

row_name = PCA_table$Gene

PCA_table$Gene = NULL

row.names(PCA_table) = row_name

PCA_data = prcomp(t(as.matrix(PCA_table)),center = TRUE)
PCA_data

summary_table = summary(PCA_data)$importance

write.csv(summary_table,"PCA_summary_table.csv")
write.csv(PCA_data$x,"PCA_result.csv")


plot(PCA_data$x[,1],PCA_data$x[,2])
pca_var = PCA_data$sdev^2

pca_var_per = round(pca_var/sum(pca_var)*100,1)
 
barplot(pca_var_per,xlab = "PC", ylab = "Percent Variation")

pca_data = data.frame(conditions = rownames(PCA_data$x),
                      X = PCA_data$x[,1],
                      Y = PCA_data$x[,2])
pca_data

ggplot(pca_data, aes(x = X, y = Y, label = conditions)) +
  xlab(paste("PC1 -", pca_var_per[1], "%", sep = "")) +
  ylab(paste("PC2 -", pca_var_per[2], "%", sep = "")) +
  theme_bw() +
  geom_text(size = 4)

pca_data = data.frame(conditions = rownames(PCA_data$x),
                      X = PCA_data$x[,3],
                      Y = PCA_data$x[,4])
pca_data

ggplot(pca_data, aes(x = X, y = Y, label = conditions)) +
  xlab(paste("PC3 -", pca_var_per[3], "%", sep = "")) +
  ylab(paste("PC4 -", pca_var_per[4], "%", sep = "")) +
  theme_bw() +
  geom_text(size = 4)