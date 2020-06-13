#Creamos un directorio para los datos, otro para los resultados y otro par las figuras
dir.create("data")
dir.create("results")
dir.create("figures")
#leemos el archivo targets
muestras <- read.csv(file.path("C:/ADO_PEC2/data/targets.csv"), header = TRUE)
#ordenamos el archivo según la columna Grupo_analisis, que contiene la información sobre el grupo 
#al que pertenece cada muestra (1=NIT, 2=SFI, 3=ELI)
muestras <- muestras[order(muestras$Grupo_analisis),]
#obtenemos los nombres de las muestras que han sido ordenadas según el grupo al que pertenecen
misMuestras <- as.character (muestras$Sample)
#elegimos al azar 10 muestras de cada grupo, teniendo en cuenta que 236 pertenecen al grupo 1, 42 al grupo 2
#y 14 al grupo 3
muestrasAzar <- c(sample(misMuestras[1:236], 10), 
                      sample(misMuestras[237:278], 10), 
                      sample(misMuestras[279:292], 10))
#creamos un nuevo data frame a partir de los datos contenidos en el archivo targets pero únicamente con las muestras
#que hemos seleccionado al azar
muestras2 <-muestras[muestras$Sample %in% muestrasAzar,]
#creamos un vector para renombrar las muestras
nuevosnombres <- c(paste0("NIT", 1:10), 
                   paste0("SFI", 1:10), 
                   paste0("ELI", 1:10))
#anadimos una columna a muestras2 con los nuevos nombres
muestras2 <- cbind(muestras2, nuevosnombres)
#convertimos el data frame "muestras2" en un archivo csv
write.table(muestras2, file.path("C:/ADO_PEC2/data/targets2.csv"), sep=";", row.names=FALSE, quote=FALSE)

#leemos el archivo que contiene los datos de expresión
contajes <- read.csv(file.path("C:/ADO_PEC2/data/counts.csv"), sep=";", row.names = 1 , header = TRUE, check.names = FALSE)
#seleccionamos las 30 muestras que hemos elegido previamente
contajes2 <- subset(contajes, select=muestrasAzar)
#ordenamos las columnas de la tabla de contaje en función de la columna Sample_Name de muestras2
muestras2orden <- as.character (muestras2$Sample)
contajes2 <- contajes2[, muestras2orden]
#ahora las filas de muestras2 (columna Sample_Name) están en el mismo orden que las columnas de contajes2

#convertimos el data frame con los contajes de las 30 muestras en un archivo csv.
write.table(contajes2, file.path("C:/ADO_PEC2/data/counts2.csv"), sep=";", row.names=TRUE, quote=FALSE, col.names = NA)

#transformamos los contajes en pseudocontajes
pseudocontajes = log2(contajes2 + 1)
#renombramos las columnas de pseudocontajes según el grupo al que pertenece cada muestra
colnames(pseudocontajes)[1:10] = paste0("NIT", 1:10)
colnames(pseudocontajes)[11:20] = paste0("SFI", 1:10)
colnames(pseudocontajes)[21:30] = paste0("ELI", 1:10)

#Exploración de los datos crudos y evaluación de la calidad

#boxplot raw
library(reshape2)
pseudocontajesmelt = melt(pseudocontajes)
condition = substr(pseudocontajesmelt$variable, 1, 3)
pseudocontajesmelt = data.frame(pseudocontajesmelt, Condition = substr(pseudocontajesmelt$variable, 1, 3))
library(ggplot2)
ggplot(pseudocontajesmelt, aes(x = variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))

#guardar boxplot raw
tiff("figures/Bloxplot_RawData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
ggplot(pseudocontajesmelt, aes(x = variable, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
dev.off()

#heatmap o CIM
library(mixOmics)
library(stats)
mat.dist = pseudocontajes
mat.dist = as.matrix(dist(t(mat.dist)))
mat.dist = mat.dist/max(mat.dist)
cim(mat.dist, symkey = FALSE, margins = c(5, 5))

#guardar heatmap raw
png("figures/Heatmap_RawData.tiff", width = 10, height = 10, units = 'in', res = 600)
cim(mat.dist, symkey = FALSE, margins = c(5, 5))
dev.off()

#PCA
library(DESeq)
library(DESeq2)
#convertimos el data frame con pseudocontajes en una matriz numérica que pueda servir de input a la función rowVars
pseudocontajes2 <- data.matrix(pseudocontajes)
#calculamos la varianza de cada fila (gen)
rv = rowVars(pseudocontajes2)
#seleccionamos los 500 genes con mayor varianza
select = order(rv, decreasing = TRUE)[1:500]
pca = prcomp(t(pseudocontajes[select, ]))
type <- factor(condition, levels = c("NIT", "SFI", "ELI"))
colType <- c("forestgreen","royalblue", "red2")[type]
plot(
  pca$x,
  col = colType,
  pch = 16,
  cex = 1.5)
legend(
  "right",
  bty = "n",
  c("NIT", "SFI", "ELI"),
  fill = c("forestgreen", "royalblue", "red2"),
  cex = 1.5)

#guardar PCA
png("figures/PCA_RawData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
plot(
  pca$x,
  col = colType,
  pch = 16,
  cex = 1.5)
legend(
  "right",
  bty = "n",
  c("NIT", "SFI", "ELI"),
  fill = c("forestgreen", "royalblue", "red2"),
  cex = 1.5)
dev.off()

#Creamos nuestro DESeqDataSet
#le damos a las filas de muestras2 el mismo nombre que las columnas de contajes2
row.names(muestras2) <- muestras2orden
countdata <- contajes2
coldata <- muestras2

dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = countdata, 
  colData   = coldata, 
  design    = ~Group
)

#filtraje
nrow(dds)
dds <- dds[ rowSums(counts(dds)) > 100, ]
nrow(dds)

#normalización
dds1 <- estimateSizeFactors(dds)
#factor de normalización
sizeFactors(dds1)
#ver valores filtrados vs filtrados y normalizados
head(counts(dds1))
head(counts(dds1, normalized=TRUE))
#guardamos los datos normalizados
normalized_counts <- counts(dds1, normalized=TRUE)

#Exploración de los datos normalizados y evaluación de la calidad
#transformamos los contajes normalizados en pseudocontajes
normalized_pseudocontajes = log2(normalized_counts + 1)
#renombramos las columnas de pseudocontajes según el grupo al que pertenece cada muestra
colnames(normalized_pseudocontajes)[1:10] = paste0("NIT", 1:10)
colnames(normalized_pseudocontajes)[11:20] = paste0("SFI", 1:10)
colnames(normalized_pseudocontajes)[21:30] = paste0("ELI", 1:10)

#boxplot datos normalizados
library(reshape2)
norm_pseudoc_melt = melt(normalized_pseudocontajes)
norm_pseudoc_melt = data.frame(norm_pseudoc_melt, Condition = substr(norm_pseudoc_melt$Var2, 1, 3))
library(ggplot2)
ggplot(norm_pseudoc_melt, aes(x = Var2, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))

#guardar boxplot norm
tiff("figures/Bloxplot_NormData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
ggplot(norm_pseudoc_melt, aes(x = Var2, y = value, fill = Condition)) + geom_boxplot() + xlab("") + ylab(expression(log[2](count + 1))) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0))
dev.off()

#heatmap o CIM datos normalizados
library(mixOmics)
library(stats)
mat.dist_norm = normalized_pseudocontajes
mat.dist_norm = as.matrix(dist(t(mat.dist_norm)))
mat.dist_norm = mat.dist_norm/max(mat.dist_norm)
cim(mat.dist_norm, symkey = FALSE, margins = c(9, 9))

#guardar heatmap norm
png("figures/Heatmap_NormData.tiff", width = 10, height = 10, units = 'in', res = 600)
cim(mat.dist_norm, symkey = FALSE, margins = c(5, 5))
dev.off()

#PCA datos normalizados
library(DESeq)
library(DESeq2)
#calculamos la varianza de cada fila (gen)
rv_norm = rowVars(normalized_pseudocontajes)
#seleccionamos los 500 genes con mayor varianza
select_norm = order(rv_norm, decreasing = TRUE)[1:500]
pca_norm = prcomp(t(normalized_pseudocontajes[select_norm, ]))
plot(
  pca_norm$x,
  col = colType,
  pch = 16,
  cex = 2.0)
legend(
  "right",
  bty = "n",
  c("NIT", "SFI", "ELI"),
  fill = c("forestgreen", "royalblue", "red2"),
  cex = 1.5)

#guardar PCA datos normalizados
png("figures/PCA_NormData.tiff", res = 200, width = 4.5, height = 4, units = 'in')
plot(
  pca_norm$x,
  col = colType,
  pch = 16,
  cex = 1.5)
legend(
  "right",
  bty = "n",
  c("NIT", "SFI", "ELI"),
  fill = c("forestgreen", "royalblue", "red2"),
  cex = 1.5)
dev.off()

#pipeline expresión diferencial
dds <- DESeq(dds)

#tabla de resultados
resultados <- results(dds, contrast=c("Group","ELI","NIT"), alpha = 0.05, lfcThreshold=1)
resultados
#resumen
summary(resultados)

#lista genes diferencialmente expresados
resSig <- subset(resultados, padj < 0.05, log2FoldChange > 1)
head(resSig[ order(resSig$padj, decreasing = FALSE), ])

#MA plot
DESeq2::plotMA(resultados, main = "ELI vs NIT", alpha = 0.05, ylim = c(-8, 8))

#guardar MA plot
png("figures/MA.tiff", width = 10, height = 10, units = 'in', res = 600)
DESeq2::plotMA(resultados, main = "ELI vs NIT", alpha = 0.05, ylim = c(-8, 8))
dev.off()

#Volcano plot
alpha <- 0.05 # Threshold on the adjusted p-value
cols <- densCols(resultados$log2FoldChange, -log10(resultados$pvalue))
plot(resultados$log2FoldChange, -log10(resultados$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~p~value),
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")

#guardar Volcano plot
png("figures/Volcano.tiff", width = 10, height = 10, units = 'in', res = 600)
plot(resultados$log2FoldChange, -log10(resultados$padj), col=cols, panel.first=grid(),
     main="Volcano plot", xlab=bquote(~Log[2]~fold~change), ylab=bquote(~-log[10]~p~value),
     pch=20, cex=0.6)
abline(v=0)
abline(v=c(-1,1), col="brown")
abline(h=-log10(alpha), col="brown")
dev.off()

#histograma de p valores
hist(resultados$pvalue[resultados$baseMean > 1], breaks = 20,
     col = "grey50", border = "white", xlab = "p-valor", ylab = "Frecuencia", main = "Histograma de p-valores")

#guardar histograma de p-valores
png("figures/histo_p_val.tiff")
hist(resultados$pvalue[resultados$baseMean > 1], breaks = 20,
     col = "grey50", border = "white", xlab = "p-valor", ylab = "Frecuencia", main = "Histograma de p-valores")
dev.off()

#agrupamiento de genes FDR < 0.01
rld = rlog(dds, blind = FALSE)
assay <- assay(rld)
colnames(assay) <- nuevosnombres
nuevosnombres2 <- c(paste0("NIT", 1:10), 
                   paste0("ELI", 1:10))
assay=subset(assay, select = nuevosnombres2)

de = (resultados$padj < 0.01)
de[is.na(de)] = FALSE

cim(t(assay[de, ]), ylab = "Samples", symkey = FALSE)

#guardar agrupamiento de genes
png("figures/agrup_genes.tiff", width = 15, height = 10, units = 'in', res = 600)
cim(t(assay[de, ]), ylab = "Samples", symkey = FALSE)
dev.off()

#anotación de resultados
library("AnnotationDbi")
library("org.Hs.eg.db")

tmp=gsub("\\..*","",row.names(resultados))

resultados$symbol <- mapIds(org.Hs.eg.db,
                     keys=tmp,
                     column="SYMBOL",
                     keytype="ENSEMBL",
                     multiVals="first")
resultados$entrez <- mapIds(org.Hs.eg.db,
                     keys=tmp,
                     column="ENTREZID",
                     keytype="ENSEMBL",
                     multiVals="first")
#mostrar tabla sin NA
head(resultados[complete.cases(resultados),])

#Análisis de significación biológica
#creamos una lista de DEGs
resultados$ENSEMBL <- tmp
genes <- resultados$padj < 0.05 & !is.na(resultados$padj)
names(genes) <- resultados$ENSEMBL
#ajustamos la Probability Weighting Function (PWF)
library(goseq)
pwf <- nullp(genes, "hg19","ensGene")
#análisis de significación biológica
goseq_resultados <- goseq(pwf, "hg19","ensGene",test.cats="GO:BP")
#análisis con clusterProfiler
library(clusterProfiler)
library(dplyr)
resSig2 <- subset(resultados, padj < 0.05)
sigGenes=gsub("\\..*","",row.names(resSig2))

enrich_go_BP <- enrichGO(
  gene= sigGenes,
  OrgDb = org.Hs.eg.db,
  keyType = "ENSEMBL",
  ont = "BP",
  universe = tmp,
  qvalueCutoff = 0.05,
  readable=TRUE
)

#Visualizar clusterProfiler resultados
dotplot(enrich_go_BP, showCategory=30)

#guardard dotplot
png("figures/dotplot.tiff", width = 15, height = 10, units = 'in', res = 600)
dotplot(enrich_go_BP, showCategory=30)
dev.off()

