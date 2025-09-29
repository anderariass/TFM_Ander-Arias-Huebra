# 02. Análisis de expresión diferencial

# Se crea la matriz de diseño
design=model.matrix(~0+group2)
colnames(design)=levels(group2)

# Se ajustan los datos
fit=lmFit(log_data_filtered, design)

# Se crea el contraste
WTD0vsOxD0=makeContrasts(WTD0-OxD0, levels = design)
fit1=contrasts.fit(fit, WTD0vsOxD0)
fit1=eBayes(fit1)
results1=topTable(fit1, adjust.method="none", number=Inf)
head(results1)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant1=results1[results1$P.Value <= 0.05 & abs(results1$logFC) >= 1, ]
nrow(significant1)

# Visualización de resultados

# Gráfico de barras
results1$DE="NO"
results1$DE[results1$logFC > 1 & results1$P.Value < 0.05]="up"
results1$DE[results1$logFC < -1 & results1$P.Value < 0.05]="down"
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
WTD0vsKOD0=makeContrasts(WTD0-KOD0, levels = design)
fit2=contrasts.fit(fit, WTD0vsKOD0)
fit2=eBayes(fit2)
results2=topTable(fit2, adjust.method="none", number=Inf)
head(results2)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant2=results2[results2$P.Value <= 0.05 & abs(results2$logFC) >= 1, ]
nrow(significant2)

# Se crea el contraste
WTDD6vsOxDD6=makeContrasts(WTDD6-OxDD6, levels = design)
fit3=contrasts.fit(fit, WTDD6vsOxDD6)
fit3=eBayes(fit3)
results3=topTable(fit3, adjust.method="none", number=Inf)
head(results3)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant3=results3[results3$P.Value <= 0.05 & abs(results3$logFC) >= 1, ]
nrow(significant3)

# Visualización de resultados

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
WTDD6vsKODD6=makeContrasts(WTDD6-KODD6, levels = design)
fit4=contrasts.fit(fit, WTDD6vsKODD6)
fit4=eBayes(fit4)
results4=topTable(fit4, adjust.method="none", number=Inf)
head(results4)

# Se filtran los metabolitos por significancia estadística (valor-p) y tasa de cambio (logFC)
significant4=results4[results4$P.Val <= 0.05 & abs(results4$logFC) >= 1, ]
nrow(significant4)

# Se crea el contraste
OxDD6vsKODD6=makeContrasts(OxDD6-KODD6, levels = design)
fit5=contrasts.fit(fit, OxDD6vsKODD6)
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