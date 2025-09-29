# 01. Análisis exploratorio

# Se importa el archivo inicial
data=read_excel("ROOT target results_formato.xlsx")
data_dataframe=as.data.frame(data)
rownames(data_dataframe)=data_dataframe[[1]]
data_dataframe=data_dataframe[, -1]

# Se calcula la cuantificación total por muestra
data_numeric=as.data.frame(lapply(data_dataframe, 
                                  function(x) as.numeric(as.character(x))))
rownames(data_numeric)=data[[1]]
total=colSums(data_numeric)
par(mar = c(7, 7, 2, 5))
barplot(total, las=2, col="skyblue", cex.axis=0.9, cex.names=1, cex.lab=1.4,
        font.lab=2, ylab="Cuantificación total", ylim=c(0, 2000))

# Se visualiza el número total de metabolitos sin cuantificación
no_cuantification=colSums(data_numeric==0)
barplot(no_cuantification, las=2, col="skyblue", 
        ylab="Número de metabolitos sin cuantificación", 
        cex.axis=0.9, cex.names=1, cex.lab=1.4, font.lab=2)

# Se visualiza el número total de metabolitos con cuantificación
cuantification=colSums(data_numeric!=0)
barplot(cuantification, las=2, col="skyblue", 
        ylab="Número de metabolitos con cuantificación", ylim=c(0, 30),
        cex.axis=0.9, cex.names=1, cex.lab=1.4, font.lab=2)

# Se realiza la transformación logarítmica
data_numeric[data_numeric<0] = 0
log_data=log2(data_numeric + 1)

# Se representa la variabilidad

# Escalado multidimensional (EMD)
colors=c("forestgreen", "gold", "orange", "deepskyblue", 
         "dodgerblue4", "purple", "firebrick", "grey", "black")
pch <- c(0, 1, 2, 3, 4, 5, 15, 16, 17)
group <- factor(c(rep("WTD0", 5), rep("WTDD6", 6), rep("WTWD6", 6),
                  rep("KOD0", 4), rep("KODD6", 6), rep("KOWD6", 6),
                  rep("OxD0", 5), rep("OxDD6", 6), rep("OxWD6", 6)))
plotMDS(log_data, col=colors[group], pch=pch[group], cex.lab=1.4, font.lab=2)
legend("topright", legend=levels(group), pch=pch, 
       col=colors, ncol=2, cex = 0.8)

# Análisis de componentes principales (ACP)
pca=prcomp(t(log_data), scale.=TRUE)
summary(pca)
var=pca$sdev^2/sum(pca$sdev^2)
porcentaje=round(var*100, 0)
pca_dataframe=data.frame(pca$x)
pca_dataframe$group=group
pca_dataframe$sample=colnames(log_data)
sd_PC1=sd(pca_dataframe$PC1)
sd_PC2=sd(pca_dataframe$PC2)
ggplot(pca_dataframe, aes(x=PC1, y=PC2, color=group, label=sample))+
  geom_point(size=3)+geom_text(vjust=-1.2, size= 2.5)+
  labs(x=paste0("PC1 (", porcentaje[1], "%)"),
       y=paste0("PC2 (", porcentaje[2], "%)"))+
  theme_minimal()+theme(legend.position="none", 
                        axis.title.x=element_text(size=15, face="bold"),
                        axis.title.y=element_text(size=15, face="bold") )
sd_PC1
sd_PC2

# Gráfico de violín
log_data_dataframe=as.data.frame(log_data)
log_data_dataframe$Metabolite=rownames(log_data_dataframe)
log_data_long=log_data_dataframe %>% 
  pivot_longer(-Metabolite, names_to="Muestra", values_to="logCuantification")
ggplot(log_data_long, aes(x=Muestra, y=logCuantification))+geom_violin(fill="lightblue")+
  theme_minimal()+theme(axis.text.x=element_text(angle=45, hjust=1), 
                        axis.title.x=element_text(size=15, face="bold"),
                        axis.title.y=element_text(size=15, face="bold"))+
  ylab("logCuantificación")


# Se eliminan los valores atípicos
outlayers=c("KOD0R5", "WTWD6R6")
log_data_filtered=log_data[, !(colnames(log_data) %in% outlayers)]
group2=factor(c(rep("WTD0", 5), rep("WTDD6", 6), rep("WTWD6", 5),
                rep("KOD0", 3), rep("KODD6", 6), rep("KOWD6", 6),
                rep("OxD0", 5), rep("OxDD6", 6), rep("OxWD6", 6)))