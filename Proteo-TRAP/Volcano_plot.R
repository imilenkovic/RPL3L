#Read the data with gene names, lfc and pvalues
data <- read.csv("ip_foldchange_pvalue.csv")

#Remove the rows with missing gene names
data2 <- subset(data, data$Gene_names!="")
rownames <- data2$Gene_names
row.names(data2) <- rownames
data3 <- subset(data2, select = -c(Gene_names))


### Volcano plot with ggplot

library(ggplot2)
library(ggrepel)

data3$diffexpressed <- "NO"
# if Fold.Change > 0.6 and pvalue < 0.05, set as "UP" 
data3$diffexpressed[data3$Fold.Change > 2 & data3$pvalue < 0.05] <- "UP"
# if Fold.Change < -0.6 and pvalue < 0.05, set as "DOWN"
data3$diffexpressed[data3$Fold.Change < -2 & data3$pvalue < 0.05] <- "DOWN"
mycolors <- c("green", "red", "gray")
names(mycolors) <- c("DOWN", "UP", "NO")

data3$delabel <- NA
data3$gene_symbol <- row.names(data3)
data3$delabel[data3$diffexpressed != "NO"] <- data3$gene_symbol[data3$diffexpressed != "NO"]

ggplot(data=data3, aes(x=Fold.Change, y=-log10(pvalue), col=diffexpressed, label=delabel)) +
  geom_point() + 
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values=c("#02C20E", "#D4DAD4", "#C92306")) +
  geom_vline(xintercept=c(-0.6, 0.6), col="#D4DAD4", linetype="dashed", size=1) +
  geom_hline(yintercept=-log10(0.05), col="#D4DAD4", linetype="dashed", size=1)

### Volcano plot with EnhancedVolcano

library(EnhancedVolcano)


### Make a dataframe with only differentially expressed proteins (pValue<0.05 and lfc>|0.6|)
data_diff <- subset(data3, data3$diffexpressed != "NO")
data_non_significant <- subset(data3, data3$diffexpressed == "NO") ### number of non-significant
data_diff_lfc <- subset(data_diff, data_diff$Fold.Change > 0.6 | data_diff$Fold.Change < -0.6) ### number of dif expressed
data_sig_low_FC <- subset(data3, data3$Fold.Change < 0.6 & data3$Fold.Change > -0.6)
data_sig_low_FC2 <- subset (data_sig_low_FC, data_sig_low_FC$pvalue < 0.05) ### number of pval<0.05 but low LFC


KO_diff <- subset(data_diff, data_diff$diffexpressed == "UP")
WT_diff <- subset(data_diff, data_diff$diffexpressed == "DOWN")

WT_diff_gene_names <- row.names(WT_diff)
WT_diff_gene_names <- as.data.frame(WT_diff_gene_names)
KO_diff_gene_names <- row.names(KO_diff)
KO_diff_gene_names <- as.data.frame(KO_diff_gene_names)



red_dots <- subset(data3, pvalue < 0.05 & (Fold.Change < -2 | Fold.Change > 2))
blue_dots <- subset(data3, pvalue<0.05 & Fold.Change>-2 & Fold.Change<2)
green_dots <- subset(data3, pvalue > 0.05 & (Fold.Change <= -2 | Fold.Change >= 2))
gray_dots <- subset(data3, (pvalue >= 0.05 & Fold.Change >= -2 & Fold.Change <= 2))

nrow(red_dots)
nrow(blue_dots)
nrow(green_dots)
nrow(gray_dots)


pdf("IP_volcano.pdf",height=9, width=11)
EnhancedVolcano(data3,
                title = 'Rpl22HA-IP',
                subtitle = "Protein abundances",
                lab = rownames(data3),
                x = 'Fold.Change',
                y = 'pvalue',
                pCutoff = 0.05,
                FCcutoff = 2,
                pointSize = 3.0,
                labSize = 4.0,
                legendLabels=c(paste('NS, n =', nrow(gray_dots)), bquote(log[2] ~ FC < "|2|" ~ "," ~ "n" == .(nrow(green_dots))), paste('p-value < 0.05, n =', nrow(blue_dots)),
                               bquote(log[2] ~ FC > "|2|" ~ "," ~ "pvalue<0.05" ~ "," ~ "n" == .(nrow(red_dots)))),
                
                legendPosition = 'right',
                legendLabSize = 16,
                legendIconSize = 5.0,)
dev.off()
