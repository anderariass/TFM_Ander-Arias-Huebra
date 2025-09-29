# 02. Análisis de expresión diferencial

# Se crea la matriz de diseño
design=model.matrix(~0+group2)
colnames(design)=levels(group2)

# Se ajustan los datos
fit=lmFit(log_data_filtered, design)

# Se crea el contraste
WT22vsOx22=makeContrasts(WT22-Ox22, levels = design)
fit1=contrasts.fit(fit, WT22vsOx22)
fit1=eBayes(fit1)
results1=topTable(fit1, adjust.method="none", number=Inf)
head(results1)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant1=results1[results1$P.Val <= 0.05 & abs(results1$logFC) >= 1, ]
nrow(significant1)

# Visualización de resultados

# Gráfico de barras
results1$DE="NO"
results1$DE[results1$logFC > 1 & results1$P.Val < 0.05]="up"
results1$DE[results1$logFC < -1 & results1$P.Val < 0.05]="down"
results1$metabolite=rownames(results1)
results1=results1[order(results1$logFC), ]

ggplot(results1, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
WT22vsKO22=makeContrasts(WT22-KO22, levels = design)
fit2=contrasts.fit(fit, WT22vsKO22)
fit2=eBayes(fit2)
results2=topTable(fit2, adjust.method="none", number=Inf)
head(results2)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant2=results2[results2$P.Val <= 0.05 & abs(results2$logFC) >= 1, ]
nrow(significant2)

# Visualización de resultados

# Gráfico de barras
results2$DE="NO"
results2$DE[results2$logFC > 1 & results2$P.Val < 0.05]="up"
results2$DE[results2$logFC < -1 & results2$P.Val < 0.05]="down"
results2$metabolite=rownames(results2)
results2=results2[order(results2$logFC), ]

ggplot(results2, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
Ox22vsKO22=makeContrasts(Ox22-KO22, levels = design)
fit3=contrasts.fit(fit, Ox22vsKO22)
fit3=eBayes(fit3)
results3=topTable(fit3, adjust.method="none", number=Inf)
head(results3)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant3=results3[results3$P.Val <= 0.05 & abs(results3$logFC) >= 1, ]
nrow(significant3)

# Gráfico de barras
results3$DE="NO"
results3$DE[results3$logFC > 1 & results3$P.Val < 0.05]="up"
results3$DE[results3$logFC < -1 & results3$P.Val < 0.05]="down"
results3$metabolite=rownames(results3)
results3=results3[order(results3$logFC), ]

ggplot(results3, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
WT28vsKO28=makeContrasts(WT28-KO28, levels = design)
fit4=contrasts.fit(fit, WT28vsKO28)
fit4=eBayes(fit4)
results4=topTable(fit4, adjust.method="none", number=Inf)
head(results4)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant4=results4[results4$P.Val <= 0.05 & abs(results4$logFC) >= 1, ]
nrow(significant4)

# Se crea el contraste
WT28vsOx28=makeContrasts(WT28-Ox28, levels = design)
fit5=contrasts.fit(fit, WT28vsOx28)
fit5=eBayes(fit5)
results5=topTable(fit5, adjust.method="none", number=Inf)
head(results5)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant5=results5[results5$P.Val <= 0.05 & abs(results5$logFC) >= 1, ]
nrow(significant5)

# Visualización de resultados

# Gráfico de barras
results5$DE="NO"
results5$DE[results5$logFC > 1 & results5$P.Val < 0.05]="up"
results5$DE[results5$logFC < -1 & results5$P.Val < 0.05]="down"
results5$metabolite=rownames(results5)
results5=results5[order(results5$logFC), ]

ggplot(results5, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
Ox28vsKO28=makeContrasts(Ox28-KO28, levels = design)
fit6=contrasts.fit(fit, Ox28vsKO28)
fit6=eBayes(fit6)
results6=topTable(fit6, adjust.method="none", number=Inf)
head(results6)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant6=results6[results6$P.Val <= 0.05 & abs(results6$logFC) >= 1, ]
nrow(significant6)

# Visualización de resultados

# Gráfico de barras
results6$DE="NO"
results6$DE[results6$logFC > 1 & results6$P.Val < 0.05]="up"
results6$DE[results6$logFC < -1 & results6$P.Val < 0.05]="down"
results6$metabolite=rownames(results6)
results6=results6[order(results6$logFC), ]

ggplot(results6, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
Psuc22vsKO22=makeContrasts(Psuc22-KO22, levels = design)
fit7=contrasts.fit(fit, Psuc22vsKO22)
fit7=eBayes(fit7)
results7=topTable(fit7, adjust.method="none", number=Inf)
head(results7)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant7=results7[results7$P.Val <= 0.05 & abs(results7$logFC) >= 1, ]
nrow(significant7)

# Visualización de resultados

# Gráfico de barras
results7$DE="NO"
results7$DE[results7$logFC > 1 & results7$P.Val < 0.05]="up"
results7$DE[results7$logFC < -1 & results7$P.Val < 0.05]="down"
results7$metabolite=rownames(results7)
results7=results7[order(results7$logFC), ]

ggplot(results7, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
Psuc22vsOx22=makeContrasts(Psuc22-Ox22, levels = design)
fit8=contrasts.fit(fit, Psuc22vsOx22)
fit8=eBayes(fit8)
results8=topTable(fit8, adjust.method="none", number=Inf)
head(results8)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant8=results8[results8$P.Val <= 0.05 & abs(results8$logFC) >= 1, ]
nrow(significant8)

# Visualización de resultados

# Gráfico de barras
results8$DE="NO"
results8$DE[results8$logFC > 1 & results8$P.Val < 0.05]="up"
results8$DE[results8$logFC < -1 & results8$P.Val < 0.05]="down"
results8$metabolite=rownames(results8)
results8=results8[order(results8$logFC), ]

ggplot(results8, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
Psuc28vsKO28=makeContrasts(Psuc28-KO28, levels = design)
fit9=contrasts.fit(fit, Psuc28vsKO28)
fit9=eBayes(fit9)
results9=topTable(fit9, adjust.method="none", number=Inf)
head(results9)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant9=results9[results9$P.Val <= 0.05 & abs(results9$logFC) >= 1, ]
nrow(significant9)

# Visualización de resultados

# Gráfico de barras
results9$DE="NO"
results9$DE[results9$logFC > 1 & results9$P.Val < 0.05]="up"
results9$DE[results9$logFC < -1 & results9$P.Val < 0.05]="down"
results9$metabolite=rownames(results9)
results9=results9[order(results9$logFC), ]

ggplot(results9, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))

# Se crea el contraste
Psuc28vsOx28=makeContrasts(Psuc28-Ox28, levels = design)
fit10=contrasts.fit(fit, Psuc28vsOx28)
fit10=eBayes(fit10)
results10=topTable(fit10, adjust.method="none", number=Inf)
head(results10)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant10=results10[results10$P.Val <= 0.05 & abs(results10$logFC) >= 1, ]
nrow(significant10)

# Visualización de resultados

# Gráfico de barras
results10$DE="NO"
results10$DE[results10$logFC > 1 & results10$P.Val < 0.05]="up"
results10$DE[results10$logFC < -1 & results10$P.Val < 0.05]="down"
results10$metabolite=rownames(results10)
results10=results10[order(results10$logFC), ]

ggplot(results10, aes(x=reorder(metabolite, logFC), y=logFC, fill=DE))+
  geom_bar(stat="identity")+coord_flip()+theme_minimal()+
  labs(x = "Metabolitos", y = "Log2 FC")+
  scale_fill_manual(values=c("up"="firebrick", "down"="steelblue", "NO"="grey70"))+
  theme(axis.title.x=element_text(size=15, face="bold", margin=margin(t=20)),
        axis.title.y=element_text(size=15, face="bold", margin=margin(r=20)),
        legend.text=element_text(size=14), 
        legend.title=element_text(size=14, face="bold"))+
  theme(plot.margin=margin(t=10, r=20, b=40, l=40))