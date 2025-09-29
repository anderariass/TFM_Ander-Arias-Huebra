# 02. Análisis de expresión diferencial

# Se crea la matriz de diseño
design=model.matrix(~0+group)
colnames(design)=levels(group)

# Se calculan las dispersiones y se ajustan los datos por cuasiverosimilitud
dge=estimateDisp(dge, design, robust = TRUE)
fit=glmQLFit(dge, design, robust = TRUE)

# Se crea el contraste
WTD0vsKOD0=makeContrasts(WTD0-KOD0, levels = design)
res6=glmQLFTest(fit, contrast=WTD0vsKOD0)
head(res6$table)

# Se corrige el testeo múltiple mediante el método de Benjamini-Hochberg
res6_corrected=topTags(res6, n=Inf)

# Se filtran los genes por significancia estadística (FDR) y tasa de cambio (logFC)
nrow(res6_corrected$table[res6_corrected$table$FDR<=0.05 
                          & abs(res6_corrected$table$logFC)>=1, ])
is.de6=decideTests(res6, adjust.method="BH", p.value=0.05, lfc = 1)
summary(is.de6)

# Visualización de resultados

# MD-plot
par(mar = c(5, 7, 2, 3))
plotMD(res6, status=is.de6, font.lab=2, cex.lab=1.4, xlab="logCPM promedio", 
       ylab="logFC")

# Gráfico de volcano
data=res6_corrected$table
ggplot(data, aes(x=logFC, y=-log10(FDR)))+geom_point()
data$DE="NO"
data$DE[data$logFC > 1 & data$FDR<0.05]="up"
data$DE[data$logFC < -1 & data$FDR<0.05]="down"
data$gene=rownames(data)
label=data$FDR <= 0.05 & abs(data$logFC) >= 1
ggplot(data, aes(x=logFC, y = -log10(FDR), col=DE))+
  geom_point(size=1)+geom_text_repel(data=data[label, ], 
                                     aes(label=gene), size=3, max.overlaps=10)+
  theme_minimal()+theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
                        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
                        legend.text=element_text(size=14), 
                        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))


# 03. Análisis de ontología génica (GO)

# Se buscan los genes sobreexpresados en el término ontológico "Proceso Biológico"
res6_table=res6_corrected$table
up6=rownames(res6_table[res6_table$FDR <= 0.05 
                        & res6_table$logFC >= 1, ])
go_up6BP=enrichGO(gene=up6, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="BP", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up6BP_1=pairwise_termsim(go_up6BP)
head(go_up6BP)

# Visualización de resultados
dotplot(go_up6BP)
barplot(go_up6BP)
emapplot(go_up6BP_1)

# Se buscan los genes sobreexpresados en el término ontológico "Función Molecular"
go_up6MF=enrichGO(gene=up6, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="MF", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up6MF_1=pairwise_termsim(go_up6MF)
head(go_up3)

# Visualización de resultados
dotplot(go_up6MF)
barplot(go_up6MF)
emapplot(go_up6MF_1)

# Se buscan los genes sobreexpresados en el término ontológico "Localización Celular"
go_up6CC=enrichGO(gene=up6, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="CC", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up6CC_1=pairwise_termsim(go_up6CC)
head(go_up3)

# Visualización de resultados
dotplot(go_up6CC)
barplot(go_up6CC)
emapplot(go_up6CC_1)

# Se buscan los genes infraexpresados de forma general
down6=rownames(res6_table[res6_table$FDR <= 0.05 
                          & res6_table$logFC <= -1, ])
go_down6=enrichGO(gene=down6, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="ALL", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_down6_1=pairwise_termsim(go_down6)
head(go_down3)

# Visualización de resultados
dotplot(go_down6)
barplot(go_down6)
emapplot(go_down6_1)

# Se crea el contraste
WTD0vsOxD0=makeContrasts(WTD0-OxD0, levels = design)
res7=glmQLFTest(fit, contrast=WTD0vsOxD0)
head(res7$table)

# Se corrige el testeo múltiple mediante el método de Benjamini-Hochberg
res7_corrected=topTags(res7, n=Inf)

# Se filtran los genes por significancia estadística (FDR) y tasa de cambio (logFC)
nrow(res7_corrected$table[res7_corrected$table$FDR<=0.05 
                          & abs(res7_corrected$table$logFC)>=1, ])
is.de7=decideTests(res7, adjust.method="BH", p.value=0.05, lfc = 1)
summary(is.de7)

# Visualización de resultados

# MD-plot
par(mar = c(5, 7, 2, 3))
plotMD(res7, status=is.de7, font.lab=2, cex.lab=1.4, xlab="logCPM promedio", 
       ylab="logFC")


# Gráfico de volcano
data=res7_corrected$table
ggplot(data, aes(x=logFC, y=-log10(FDR)))+geom_point()
data$DE="NO"
data$DE[data$logFC > 1 & data$FDR<0.05]="up"
data$DE[data$logFC < -1 & data$FDR<0.05]="down"
data$gene=rownames(data)
label=data$FDR <= 0.05 & abs(data$logFC) >= 1
ggplot(data, aes(x=logFC, y = -log10(FDR), col=DE))+
  geom_point(size=1)+geom_text_repel(data=data[label, ], 
                                     aes(label=gene), size=3, max.overlaps=10)+
  theme_minimal()+theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
                        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
                        legend.text=element_text(size=14), 
                        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))


# 03. Análisis de ontología génica (GO)

# Se buscan los genes sobreexpresados en el término ontológico "Proceso Biológico"
res7_table=res7_corrected$table
up7=rownames(res7_table[res7_table$FDR <= 0.05 
                        & res7_table$logFC >= 1, ])
go_up7BP=enrichGO(gene=up7, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="BP", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up7BP_1=pairwise_termsim(go_up7BP)
head(go_up7BP)

# Visualización de resultados
dotplot(go_up7BP)
barplot(go_up7BP)
emapplot(go_up7BP_1)

# Se buscan los genes sobreexpresados en el término ontológico "Función Molecular"
go_up7MF=enrichGO(gene=up7, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="MF", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up7MF_1=pairwise_termsim(go_up7MF)
head(go_up7)

# Visualización de resultados
dotplot(go_up7MF)
barplot(go_up7MF)
emapplot(go_up7MF_1)

# Se buscan los genes sobreexpresados en el término ontológico "Localización Celular"
go_up7CC=enrichGO(gene=up7, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="CC", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up7CC_1=pairwise_termsim(go_up7CC)
head(go_up7)

# Visualización de resultados
dotplot(go_up7CC)
barplot(go_up7CC)
emapplot(go_up7CC_1)

# Se buscan los genes infraexpresados de forma general
down7=rownames(res7_table[res7_table$FDR <= 0.05 
                          & res7_table$logFC <= -1, ])
go_down7=enrichGO(gene=down7, OrgDb=org.At.tair.db, keyType="TAIR", 
                  ont="ALL", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_down7_1=pairwise_termsim(go_down7)
head(go_down7)

# Visualización de resultados
dotplot(go_down7)
barplot(go_down7)
emapplot(go_down7_1)

# Se crea el contraste
WTDD6vsOxDD6=makeContrasts(WTDD6-OxDD6, levels = design)
res18=glmQLFTest(fit, contrast=WTDD6vsOxDD6)
head(res18$table)

# Se corrige el testeo múltiple mediante el método de Benjamini-Hochberg
res18_corrected=topTags(res18, n=Inf)

# Se filtran los genes por significancia estadística (FDR) y tasa de cambio (logFC)
nrow(res18_corrected$table[res18_corrected$table$FDR<=0.05 
                           & abs(res18_corrected$table$logFC)>=1, ])
is.de18=decideTests(res18, adjust.method="BH", p.value=0.05, lfc = 1)
summary(is.de18)

# Visualización de resultados

# MD-plot
par(mar = c(5, 7, 2, 3))
plotMD(res18, status=is.de18, font.lab=2, cex.lab=1.4, xlab="logCPM promedio", 
       ylab="logFC")

# Gráfico de volcano
data=res18_corrected$table
ggplot(data, aes(x=logFC, y=-log10(FDR)))+geom_point()
data$DE="NO"
data$DE[data$logFC > 1 & data$FDR<0.05]="up"
data$DE[data$logFC < -1 & data$FDR<0.05]="down"
data$gene=rownames(data)
label=data$FDR <= 0.05 & abs(data$logFC) >= 1
ggplot(data, aes(x=logFC, y = -log10(FDR), col=DE))+
  geom_point(size=1)+geom_text_repel(data=data[label, ], 
                                     aes(label=gene), size=3, max.overlaps=10)+
  theme_minimal()+theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
                        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
                        legend.text=element_text(size=14), 
                        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))


# 03. Análisis de ontología génica (GO)

# Se buscan los genes sobreexpresados en el término ontológico "Proceso Biológico"
res18_table=res18_corrected$table
up18=rownames(res18_table[res18_table$FDR <= 0.05 
                          & res18_table$logFC >= 1, ])
go_up18BP=enrichGO(gene=up18, OrgDb=org.At.tair.db, keyType="TAIR", 
                   ont="BP", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up18BP_1=pairwise_termsim(go_up18BP)
head(go_up6BP)

# Visualización de resultados
dotplot(go_up18BP)
barplot(go_up18BP)
emapplot(go_up18BP_1)

# Se buscan los genes sobreexpresados en el término ontológico "Función Molecular"
go_up18MF=enrichGO(gene=up18, OrgDb=org.At.tair.db, keyType="TAIR", 
                   ont="MF", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up18MF_1=pairwise_termsim(go_up18MF)
head(go_up3)

# Visualización de resultados
dotplot(go_up18MF)
barplot(go_up18MF)
emapplot(go_up18MF_1)

# Se buscan los genes sobreexpresados en el término ontológico "Localización Celular"
go_up18CC=enrichGO(gene=up18, OrgDb=org.At.tair.db, keyType="TAIR", 
                   ont="CC", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up18CC_1=pairwise_termsim(go_up18CC)
head(go_up3)

# Visualización de resultados
dotplot(go_up18CC)
barplot(go_up18CC)
emapplot(go_up18CC_1)

# Se buscan los genes infraexpresados en el término ontológico "Proceso Biológico"
down18=rownames(res18_table[res18_table$FDR <= 0.05 
                            & res18_table$logFC <= -1, ])
go_down18BP=enrichGO(gene=down18, OrgDb=org.At.tair.db, keyType="TAIR", 
                     ont="BP", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_down18BP_1=pairwise_termsim(go_down18BP)
head(go_down3)

# Visualización de resultados
dotplot(go_down18BP)
barplot(go_down18BP)
emapplot(go_down18BP_1)

# Se buscan los genes infraexpresados en el término ontológico "Función Molecular"
go_down18MF=enrichGO(gene=down18, OrgDb=org.At.tair.db, keyType="TAIR", 
                     ont="MF", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_down18MF_1=pairwise_termsim(go_down18MF)
head(go_down3)

# Visualización de resultados
dotplot(go_down18MF)
barplot(go_down18MF)
emapplot(go_down18MF_1)

# Se crea el contraste
WTDD6vsKODD6=makeContrasts(WTDD6-KODD6, levels = design)
res21=glmQLFTest(fit, contrast=WTDD6vsKODD6)
head(res21$table)

# Se corrige el testeo múltiple mediante el método de Benjamini-Hochberg
res21_corrected=topTags(res21, n=Inf)

# Se filtran los genes por significancia estadística (FDR) y tasa de cambio (logFC)
nrow(res21_corrected$table[res21_corrected$table$FDR<=0.05 
                           & abs(res21_corrected$table$logFC)>=1, ])
is.de21=decideTests(res21, adjust.method="BH", p.value=0.05, lfc = 1)
summary(is.de21)

# Visualización de resultados

# MD-plot
par(mar = c(5, 7, 2, 3))
plotMD(res21, status=is.de21, font.lab=2, cex.lab=1.4, xlab="logCPM promedio", 
       ylab="logFC")

# Gráfico de volcano
data=res21_corrected$table
ggplot(data, aes(x=logFC, y=-log10(FDR)))+geom_point()
data$DE="NO"
data$DE[data$logFC > 1 & data$FDR<0.05]="up"
data$DE[data$logFC < -1 & data$FDR<0.05]="down"
data$gene=rownames(data)
label=data$FDR <= 0.05 & abs(data$logFC) >= 1
ggplot(data, aes(x=logFC, y = -log10(FDR), col=DE))+
  geom_point(size=1)+geom_text_repel(data=data[label, ], 
                                     aes(label=gene), size=3, max.overlaps=10)+
  theme_minimal()+theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
                        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
                        legend.text=element_text(size=14), 
                        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))


# 03. Análisis de ontología génica (GO)

# Se buscan los genes sobreexpresados en el término ontológico "Proceso Biológico"
res21_table=res21_corrected$table
up21=rownames(res21_table[res21_table$FDR <= 0.05 
                          & res21_table$logFC >= 1, ])
go_up21BP=enrichGO(gene=up21, OrgDb=org.At.tair.db, keyType="TAIR", 
                   ont="BP", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up21BP_1=pairwise_termsim(go_up21BP)
head(go_up6BP)

# Visualización de resultados
dotplot(go_up21BP)
barplot(go_up21BP)
emapplot(go_up21BP_1)

# Se buscan los genes sobreexpresados en el término ontológico "Función Molecular"
go_up21MF=enrichGO(gene=up21, OrgDb=org.At.tair.db, keyType="TAIR", 
                   ont="MF", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up21MF_1=pairwise_termsim(go_up21MF)
head(go_up3)

# Visualización de resultados
dotplot(go_up21MF)
barplot(go_up21MF)
emapplot(go_up21MF_1)

# Se buscan los genes sobreexpresados en el término ontológico "Localización Celular"
go_up21CC=enrichGO(gene=up21, OrgDb=org.At.tair.db, keyType="TAIR", 
                   ont="CC", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_up21CC_1=pairwise_termsim(go_up21CC)
head(go_up3)

# Visualización de resultados
dotplot(go_up21CC)
barplot(go_up21CC)
emapplot(go_up21CC_1)

# Se buscan los genes infraexpresados en el término ontológico "Proceso Biológico"
down21=rownames(res21_table[res21_table$FDR <= 0.05 
                            & res21_table$logFC <= -1, ])
go_down21BP=enrichGO(gene=down21, OrgDb=org.At.tair.db, keyType="TAIR", 
                     ont="BP", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_down21BP_1=pairwise_termsim(go_down21BP)
head(go_down3)

# Visualización de resultados
dotplot(go_down21BP)
barplot(go_down21BP)
emapplot(go_down21BP_1)

# Se buscan los genes infraexpresados en el término ontológico "Función Molecular"
go_down21MF=enrichGO(gene=down21, OrgDb=org.At.tair.db, keyType="TAIR", 
                     ont="MF", pAdjustMethod ="BH", pvalueCutoff=0.1, qvalueCutoff=0.1)
go_down21MF_1=pairwise_termsim(go_down21MF)
head(go_down3)

# Visualización de resultados
dotplot(go_down21MF)
barplot(go_down21MF)
emapplot(go_down21MF_1)

# Se crea el contraste
OxDD6vsKODD6=makeContrasts(OxDD6-KODD6, levels = design)
res33=glmQLFTest(fit, contrast=OxDD6vsKODD6)
head(res33$table)

# Se corrige el testeo múltilple mediante el método de Benjamini-Hochberg
res33_corrected=topTags(res33, n=Inf)

# Se filtran los genes por significancia estadística (FDR) y tasa de cambio (logFC)
nrow(res33_corrected$table[res33_corrected$table$FDR<=0.05 
                           & abs(res33_corrected$table$logFC)>=1, ])
is.de33=decideTests(res33, adjust.method="BH", p.value=0.05, lfc = 1)
summary(is.de33)

# Visualización de resultados
# MD-plot
par(mar = c(5, 7, 2, 3))
plotMD(res33, status=is.de33, font.lab=2, cex.lab=1.4, xlab="logCPM promedio", 
       ylab="logFC")

# Gráfico de volcano
data=res33_corrected$table
ggplot(data, aes(x=logFC, y=-log10(FDR)))+geom_point()
data$DE="NO"
data$DE[data$logFC > 1 & data$FDR<0.05]="up"
data$DE[data$logFC < -1 & data$FDR<0.05]="down"
data$gene=rownames(data)
label=data$FDR <= 0.05 & abs(data$logFC) >= 1
ggplot(data, aes(x=logFC, y = -log10(FDR), col=DE))+
  geom_point(size=1)+geom_text_repel(data=data[label, ], 
                                     aes(label=gene), size=3, max.overlaps=10)+
  theme_minimal()+theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
                        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
                        legend.text=element_text(size=14), 
                        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))