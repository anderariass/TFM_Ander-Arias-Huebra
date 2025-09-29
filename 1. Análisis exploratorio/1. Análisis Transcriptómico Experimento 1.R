# 01. Análisis exploratorio

# Se importa la matriz de conteos inicial
data=read.delim("raw_counts_matrix_names.txt")

# Se modifican los nombres de los genes
names=useMart("plants_mart", dataset="athaliana_eg_gene", 
              host= "https://plants.ensembl.org")
gene_ids=rownames(data)
change=getBM(attributes=c("tair_locus", "external_gene_name"),
             filters="tair_locus", values=gene_ids, mart=names)
position=match(gene_ids, change$tair_locus)
final_names=ifelse(is.na(position), gene_ids, change$external_gene_name[position])
rownames(data)=make.unique(final_names)

# Se calcula el número total de lecturas por muestra
reads=colSums(data)
barplot(reads, las=2, cex.names=1, cex.axis=0.9, cex.lab=1.4, 
        mgp=c(5,0.7,0),
        font.lab=2, col="skyblue", xlab="Muestra", ylab="Número de lecturas")
par(mar = c(7, 7, 2, 5))
barplot(reads, las=2, cex.names=0.7, cex.axis=0.7, mgp=c(5,0.7,0), col="skyblue"
        , xlab="Muestra", ylab="Número de lecturas")

# Se visualiza el número total de genes sin expresión
no_expression=apply(data, 2, function (x) sum(x==0))
barplot(no_expression, las=2, col="skyblue", cex.names=1, cex.axis=0.9, 
        cex.lab=1.4, 
        font.lab=2, mgp=c(5,0.7,0), xlab="Muestra", ylab="Número de genes con 0 lecturas")

# Se visualiza el número total de genes con expresión
expression=apply(data, 2, function(x) sum(x>0))
barplot(expression, ylim=c(0, 25000), las=2, col="skyblue", cex.names=1, 
        cex.axis=0.9, cex.lab=1.4, font.lab=2, mgp=c(5,0.7,0), 
        xlab="Muestra", ylab="Número de genes con > 0 lecturas")

# Se crea el objeto DGElist
metadata <- factor(c(
  "WTD0", "WTD0", 
  "WTWD6", "WTWD6", 
  "WTDD6", "WTDD6",
  "OxD0", "OxD0", 
  "OxWD6", "OxWD6", 
  "OxDD6", "OxDD6",
  "KOD0", "KOD0", "KOD0",
  "KOWD6", "KOWD6", "KOWD6",
  "KODD6", "KODD6", "KODD6"
))
table(metadata)
dge=DGEList(counts=data, group=metadata)
View(dge)
dge$counts
dge$samples

# Se filtran los genes con baja expresión
keep=filterByExpr(dge)
summary(keep)
dge=dge[keep, keep.lib.sizes=FALSE]

# Se realiza la normalización Trimmed Mean of M-values (TMM)
dge=calcNormFactors(dge)
dge$samples

# Se calcula el logCPM
logcpm=cpm(dge, log=TRUE)
par(font.lab=2)
boxplot(logcpm, las=2, col="skyblue", cex.lab=1.4, 
        font.lab=2, cex.names=1, cex.axis=0.9, ylab="logCPM")

dge_sin_normalizar=DGEList(counts=data, group=metadata)
logcpm_sin_filtrar_normalizar=cpm(dge_sin_normalizar, log=TRUE)
boxplot(logcpm_sin_filtrar_normalizar, las=2, col="skyblue", 
        cex.names=1, cex.axis=0.9, cex.lab=1.4, ylab="logCPM")
keep=filterByExpr(dge_sin_normalizar)
dge_sin_normalizar=dge_sin_normalizar[keep, keep.lib.sizes=FALSE]
logcpm_sin_normalizar=cpm(dge_sin_normalizar, log=TRUE)
boxplot(logcpm_sin_normalizar, las=2, col="skyblue", 
        cex.names=1, cex.axis=0.9, cex.lab=1.4, ylab="logCPM")

# Se representa la variabilidad

# Escalado multidimensional (EMD)
colors=c("forestgreen", "gold", "orange", "deepskyblue", 
         "dodgerblue4", "purple", "firebrick", "grey", "black")
pch <- c(0, 1, 2, 3, 4, 5, 15, 16, 17)
group <- factor(c(rep("WTD0", 2), rep("WTWD6", 2), rep("WTDD6", 2), 
                  rep("OxD0", 2), rep("OxWD6", 2), rep("OxDD6", 2), 
                  rep("KOD0", 3), rep("KOWD6", 3), rep("KODD6", 3)))
plotMDS(dge, col=colors[group], pch=pch[group], cex.lab=1.4, font.lab=2)
legend("bottomright", legend=levels(group), pch=pch, 
       col=colors, ncol=2, cex = 0.6)

# Análisis de componentes principales (ACP)
pca=prcomp(t(logcpm), scale.=TRUE)
summary(pca)
var=pca$sdev^2/sum(pca$sdev^2)
porcentaje=round(var*100, 0)
pca_dataframe=data.frame(pca$x)
pca_dataframe$group=group
pca_dataframe$sample=colnames(logcpm)
ggplot(pca_dataframe, aes(x=PC1, y=PC2, color=group, label=sample))+
  geom_point(size=3)+geom_text(vjust=-1.2, size= 2.5)+
  labs(x=paste0("PC1 (", porcentaje[1], "%)"),
       y=paste0("PC2 (", porcentaje[2], "%)"))+
  theme_minimal()+theme(legend.position="none", 
                        axis.title.x=element_text(size=15, face="bold"),
                        axis.title.y=element_text(size=15, face="bold") )

# Gráfico de violín
logcpm_dataframe=as.data.frame(logcpm)
logcpm_dataframe$Gene=rownames(logcpm_dataframe)
logcpm_long=logcpm_dataframe %>% 
  pivot_longer(-Gene, names_to="Muestra", values_to="logcpm")
ggplot(logcpm_long, aes(x=Muestra, y=logcpm))+geom_violin(fill="lightblue")+
  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1), 
                        axis.title.x=element_text(size=15, face="bold"),
                        axis.title.y=element_text(size=15, face="bold"))+
  ylab("logCPM")