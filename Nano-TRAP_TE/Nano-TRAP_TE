## 1. LOAD DATA
counts<-read.csv("Nano_TRAP_full_counts_renormalized.csv", row.names=1)
colnames(counts)<-c(
  "S1","S2","S3","S4","S5","S6",
  "S7","S8","S9","S10","S11","S12")  
# Exclude those that are not expressed in all conditions at RNA level
subset1<-subset(counts, S2 >= 10 & S4 >=10,) 
subset1<-subset(subset1, S6 >= 10 & S8 >=10,) 
subset1<-subset(subset1, 10 >= 10 & S2 >=10,) 
filtered_counts<-subset1

colnames(filtered_counts)<-c(
  "KO.HP.rep1","KO.tot.rep1","KO.HP.rep2","KO.tot.rep2","KO.HP.rep3","KO.tot.rep3",
  "WT.HP.rep1","WT.tot.rep1","WT.HP.rep2","WT.tot.rep2","WT.HP.rep3","WT.tot.rep3")

# Check
dim(counts) # 17873
dim(filtered_counts) # 11231
counts<-filtered_counts ##Â USE FILTERED COUNTS

# Add pseudocount
counts<-counts+1


# Check for replicability
plot_denscols<-function(pdfname,my_x,my_y,xlab, ylab) {
  pdf(file=pdfname, height=6, width=6)
  dcols<-densCols(my_x,my_y, colramp=colorRampPalette(blues9[-(1:3)]))
  plot(my_x,my_y,col=dcols,cex=1, cex.lab=1,cex.main=3,lwd=5,pch=20,xlab=xlab,ylab=ylab)
  title(main=pdfname, col.main="black", font.main=4)
  abline(0,1, lty=2)
  # Correlation
  test<-cor.test(my_x,my_y, method="pearson")
  print(test)
  cor2<-cor(my_x,my_y, method="pearson")
  cor22<-format(round(cor2, 3), nsmall = 2)   #to produce value with two digits
  cor222<-paste("Pearson corr =",cor22)
  pval<-paste("Pvalue =",test$p.value)
  mtext(cor222) #Print the subtitle with the dataset correlation
  dev.off()
}
## 2. ADD TE
################

# add TE
counts$TE.KO.rep1<-counts$KO.HP.rep1/counts$KO.tot.rep1
counts$TE.KO.rep2<-counts$KO.HP.rep2/counts$KO.tot.rep2
counts$TE.KO.rep3<-counts$KO.HP.rep3/counts$KO.tot.rep3
counts$TE.WT.rep1<-counts$WT.HP.rep1/counts$WT.tot.rep1
counts$TE.WT.rep2<-counts$WT.HP.rep2/counts$WT.tot.rep2
counts$TE.WT.rep3<-counts$WT.HP.rep3/counts$WT.tot.rep3

#add_deltaTE
counts$delta_TE.rep1<-counts$TE.KO.rep1-counts$TE.WT.rep1
counts$delta_TE.rep2<-counts$TE.KO.rep2-counts$TE.WT.rep2
counts$delta_TE.rep3<-counts$TE.KO.rep3-counts$TE.WT.rep3

# add mean
counts$mean_deltaTE<-rowMeans(counts[c("delta_TE.rep1","delta_TE.rep2","delta_TE.rep3")])

# add Ratio
counts$meanTE_KO<-rowMeans(counts[c("TE.KO.rep1","TE.KO.rep2","TE.KO.rep3")])
counts$meanTE_WT<-rowMeans(counts[c("TE.WT.rep1","TE.WT.rep2","TE.WT.rep3")])
counts$meanTE_RATIOFINAL<-counts$meanTE_KO/counts$meanTE_WT
counts$log2_foldchangeTE<-log2(counts$meanTE_RATIOFINAL)

# Inspect
head(counts)

# Plot
plot_denscols("PDFs/3.TEcomparison/KO1_KO2_TE.pdf",log(counts$TE.KO.rep1),log(counts$TE.KO.rep2),"TE.KO.rep1","TE.KO.rep2")# rho=0.9624806; p-value < 2.2e-16

plot_denscols("PDFs/3.TEcomparison/KO1_KO2_TE.pdf",log(counts$TE.KO.rep1),log(counts$TE.KO.rep2),"TE.KO.rep1","TE.KO.rep2")# rho=0.9624806; p-value < 2.2e-16

plot_denscols("PDFs/3.TEcomparison/WT1_WT2_TE.pdf",log(counts$TE.WT.rep1),log(counts$TE.WT.rep2),"TE.WT.rep1","TE.WT.rep2")# rho=0.9148412; p-value < 2.2e-16

plot_denscols("PDFs/3.TEcomparison/WT1_KO1_TE.pdf",log(counts$TE.WT.rep1),log(counts$TE.KO.rep1),"TE.WT.rep1","TE.KO.rep1")# rho=0.9155115 ; p-value < 2.2e-16

plot_denscols("PDFs/3.TEcomparison/WT2_KO2_TE.pdf",log(counts$TE.WT.rep2),log(counts$TE.KO.rep2),"TE.WT.rep2","TE.KO.rep2")# rho=0.9581443; p-value < 2.2e-16

plot_denscols("PDFs/3.TEcomparison/WT3_KO3_TE.pdf",log(counts$TE.WT.rep3),log(counts$TE.KO.rep3),"TE.WT.rep3","TE.KO.rep3")# rho=0.9682729; p-value < 2.2e-16


plot_denscols("PDFs/3.TEcomparison/WT_vs_KO_TE.pdf",log(counts$meanTE_WT),log(counts$meanTE_KO),"TE.WT","TE.KO")# rho=0.9682729; p-value < 2.2e-16

write.table(counts,file="Nano_TRAP_TE_analysis.txt",quote=F,row.names=T)
