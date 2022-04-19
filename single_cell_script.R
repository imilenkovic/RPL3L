require("devtools")
library("stringr")

#devtools::install_github("tallulandrews/M3Drop")


dat = read.delim("GSE120064_TAC_raw_umi_matrix.csv", sep=",", header=TRUE)
dat[1:5,1:5]

annotation = read.delim("GSE120064_TAC_clean_cell_info_summary.txt", header=TRUE)
#annotation <- annotation[match(cellIDs, annotation[,1]),]


dim(dat)

rownames(dat) <- dat[,1]
dat <- dat[,-1]

cellIDs <- colnames(dat)
cell_info <- strsplit(cellIDs, "\\_")


Plate <- annotation$plate_id
Plate <- unlist(Plate)

Mouse <- annotation$sample
Mouse <- unlist(Mouse)

Well <- annotation$CellID
Well <- as.data.frame(Well)
#Well <- word(Well$Well,3,4,sep = "\\_")
Well <- unlist(Well)


summary(factor(Mouse))

table(Mouse, Plate)

celltype <- annotation[,8]

library("SingleCellExperiment")
library("scater")
library("data.table")
library("tidyverse")

cell_anns <- data.frame(mouse = Mouse, well=Well, type=celltype)

week0 <- cell_anns[cell_anns$mouse %like% "0w", ]
week2 <- cell_anns[cell_anns$mouse %like% "2w", ]
week5 <- cell_anns[cell_anns$mouse %like% "5w", ]
week8 <- cell_anns[cell_anns$mouse %like% "8w", ]
week11 <- cell_anns[cell_anns$mouse %like% "11w", ]


#dat_week0 <- dat[, c(names(dat) %in% week0$well)]
#dat_week2 <- dat[, (names(dat) %in% week2$well)]
#dat_week5 <- dat[, (names(dat) %in% week5$well)]
#dat_week8 <- dat[, (names(dat) %in% week8$well)]
#dat_week11 <- dat[, (names(dat) %in% week11$well)]


dat_week0 <- dat %>% select(week0$well)
dat_week2 <- dat %>% select(week2$well)
dat_week5 <- dat %>% select(week5$well)
dat_week8 <- dat %>% select(week8$well)
dat_week11 <- dat %>% select(week11$well)


rownames(cell_anns) <- colnames(dat)

rownames(week0) <- week0$well
rownames(week2) <- week2$well
rownames(week5) <- week5$well
rownames(week8) <- week8$well
rownames(week11) <- week11$well 

sceset <- SingleCellExperiment(assays = list(counts = as.matrix(dat)), colData=cell_anns)

sceset_week0 <- SingleCellExperiment(assays = list(counts = as.matrix(dat_week0)), colData=week0)
sceset_week2 <- SingleCellExperiment(assays = list(counts = as.matrix(dat_week2)), colData=week2)
sceset_week5 <- SingleCellExperiment(assays = list(counts = as.matrix(dat_week5)), colData=week5)
sceset_week8 <- SingleCellExperiment(assays = list(counts = as.matrix(dat_week8)), colData=week8)
sceset_week11 <- SingleCellExperiment(assays = list(counts = as.matrix(dat_week11)), colData=week11)



example_sce <- addPerCellQC(sceset, 
                            subsets=list(Mito=grep("mt-", rownames(sceset))))

example_sce_w0 <- addPerCellQC(sceset_week0, 
                            subsets=list(Mito=grep("mt-", rownames(sceset_week0))))
example_sce_w2 <- addPerCellQC(sceset_week2, 
                               subsets=list(Mito=grep("mt-", rownames(sceset_week2))))
example_sce_w5 <- addPerCellQC(sceset_week5, 
                               subsets=list(Mito=grep("mt-", rownames(sceset_week5))))
example_sce_w8 <- addPerCellQC(sceset_week8, 
                               subsets=list(Mito=grep("mt-", rownames(sceset_week8))))
example_sce_w11 <- addPerCellQC(sceset_week11, 
                               subsets=list(Mito=grep("mt-", rownames(sceset_week11))))




#plotColData(example_sce, x = "sum", y="detected", colour_by="mouse") 
#plotHighestExprs(example_sce, exprs_values = "counts")

###logNormCounts

example_sce <- logNormCounts(example_sce)
example_sce_w0 <- logNormCounts(example_sce_w0)
example_sce_w2 <- logNormCounts(example_sce_w2)
example_sce_w5 <- logNormCounts(example_sce_w5)
example_sce_w8 <- logNormCounts(example_sce_w8)
example_sce_w11 <- logNormCounts(example_sce_w11)

###Run PCA

example_sce <- runPCA(example_sce)
example_sce_w0 <- runPCA(example_sce_w0)
example_sce_w2 <- runPCA(example_sce_w2)
example_sce_w5 <- runPCA(example_sce_w5)
example_sce_w8 <- runPCA(example_sce_w8)
example_sce_w11 <- runPCA(example_sce_w11)

### Reduce dimensions

str(reducedDim(example_sce, "PCA"))
str(reducedDim(example_sce_w0, "PCA"))
str(reducedDim(example_sce_w2, "PCA"))
str(reducedDim(example_sce_w5, "PCA"))
str(reducedDim(example_sce_w8, "PCA"))
str(reducedDim(example_sce_w11, "PCA"))


### Run PCA2

example_sce <- runPCA(example_sce, name="PCA2",
                      subset_row=rownames(example_sce)[1:1000],
                      ncomponents=25)
example_sce_w0 <- runPCA(example_sce_w0, name="PCA2",
                      subset_row=rownames(example_sce_w0)[1:1000],
                      ncomponents=25)
example_sce_w2 <- runPCA(example_sce_w2, name="PCA2",
                      subset_row=rownames(example_sce_w2)[1:1000],
                      ncomponents=25)
example_sce_w5 <- runPCA(example_sce_w5, name="PCA2",
                      subset_row=rownames(example_sce_w5)[1:1000],
                      ncomponents=25)
example_sce_w8 <- runPCA(example_sce_w8, name="PCA2",
                      subset_row=rownames(example_sce_w8)[1:1000],
                      ncomponents=25)
example_sce_w11 <- runPCA(example_sce_w11, name="PCA2",
                         subset_row=rownames(example_sce_w11)[1:1000],
                         ncomponents=25)

### Reduce dimensions in PCA2

str(reducedDim(example_sce, "PCA2"))
str(reducedDim(example_sce_w0, "PCA2"))
str(reducedDim(example_sce_w2, "PCA2"))
str(reducedDim(example_sce_w5, "PCA2"))
str(reducedDim(example_sce_w8, "PCA2"))
str(reducedDim(example_sce_w11, "PCA2"))

### TSNE

set.seed(1000)
example_sce <- runTSNE(example_sce, perplexity=10)
example_sce_w0 <- runTSNE(example_sce_w0, perplexity=10)
example_sce_w2 <- runTSNE(example_sce_w2, perplexity=10)
example_sce_w5 <- runTSNE(example_sce_w5, perplexity=10)
example_sce_w8 <- runTSNE(example_sce_w8, perplexity=10)
example_sce_w11 <- runTSNE(example_sce_w11, perplexity=10)


set.seed(1000)
example_sce <- runTSNE(example_sce, perplexity=50, 
                       dimred="PCA", n_dimred=10)
example_sce_w0 <- runTSNE(example_sce_w0, perplexity=50, 
                       dimred="PCA", n_dimred=10)
example_sce_w2 <- runTSNE(example_sce_w2, perplexity=50, 
                       dimred="PCA", n_dimred=10)
example_sce_w5 <- runTSNE(example_sce_w5, perplexity=50, 
                       dimred="PCA", n_dimred=10)
example_sce_w8 <- runTSNE(example_sce_w8, perplexity=50, 
                       dimred="PCA", n_dimred=10)
example_sce_w11 <- runTSNE(example_sce_w11, perplexity=50, 
                       dimred="PCA", n_dimred=10)

### UMAP

example_sce <- runUMAP(example_sce)
example_sce_w0 <- runUMAP(example_sce_w0)
example_sce_w2 <- runUMAP(example_sce_w2)
example_sce_w5 <- runUMAP(example_sce_w5)
example_sce_w8 <- runUMAP(example_sce_w8)
example_sce_w11 <- runUMAP(example_sce_w11)


### Plot PCA

pdf("PCA_full_by_cell_type.pdf",height=9, width=11)
plotReducedDim(example_sce, dimred = "PCA", colour_by = "type")
dev.off()

pdf("PCA_full_by_animal.pdf",height=9, width=11)
plotReducedDim(example_sce, dimred = "PCA", colour_by = "mouse")
dev.off()

pdf("PCA_w0.pdf",height=9, width=11)
plotReducedDim(example_sce_w0, dimred = "PCA", colour_by = "type")
dev.off()

pdf("PCA_w2.pdf",height=9, width=11)
plotReducedDim(example_sce_w2, dimred = "PCA", colour_by = "type")
dev.off()

pdf("PCA_w5.pdf",height=9, width=11)
plotReducedDim(example_sce_w5, dimred = "PCA", colour_by = "type")
dev.off()

pdf("PCA_w8.pdf",height=9, width=11)
plotReducedDim(example_sce_w8, dimred = "PCA", colour_by = "type")
dev.off()

pdf("PCA_w11.pdf",height=9, width=11)
plotReducedDim(example_sce_w11, dimred = "PCA", colour_by = "type")
dev.off()

### Plot TSNE - by cell type (and mouse for full)

pdf("TSNE_full_by_type.pdf",height=9, width=11)
plotTSNE(example_sce, colour_by = "type")
dev.off()

pdf("TSNE_full_by_mouse.pdf",height=9, width=11)
plotTSNE(example_sce, colour_by = "mouse")
dev.off()

pdf("TSNE_w0_cell_type.pdf",height=9, width=11)
plotTSNE(example_sce_w0, colour_by = "type")
dev.off()

pdf("TSNE_w2_cell_type.pdf",height=9, width=11)
plotTSNE(example_sce_w2, colour_by = "type")
dev.off()

pdf("TSNE_w5_cell_type.pdf",height=9, width=11)
plotTSNE(example_sce_w5, colour_by = "type")
dev.off()

pdf("TSNE_w8_cell_type.pdf",height=9, width=11)
plotTSNE(example_sce_w8, colour_by = "type")
dev.off()

pdf("TSNE_w11_cell_type.pdf",height=9, width=11)
plotTSNE(example_sce_w11, colour_by = "type")
dev.off()

### Plot TSNE - by Rpl3l

pdf("TSNE_full_Rpl3l.pdf",height=9, width=11)
plotTSNE(example_sce, colour_by = "Rpl3l") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w0_Rpl3l.pdf",height=9, width=11)
plotTSNE(example_sce_w0, colour_by = "Rpl3l") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w2_Rpl3l.pdf",height=9, width=11)
plotTSNE(example_sce_w2, colour_by = "Rpl3l") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w5_Rpl3l.pdf",height=9, width=11)
plotTSNE(example_sce_w5, colour_by = "Rpl3l") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w8_Rpl3l.pdf",height=9, width=11)
plotTSNE(example_sce_w8, colour_by = "Rpl3l") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w11_Rpl3l.pdf",height=9, width=11)
plotTSNE(example_sce_w11, colour_by = "Rpl3l") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

#3-colour gradient:
#plotTSNE(example_sce_w0, colour_by = "Rpl3l") + scale_color_gradient2(low = "#FFACF1", mid = "#86A8E7", high = "#00DEFF",midpoint = 3, limits = c(0, 6))



#### Plot TSNE - by Rpl3


pdf("TSNE_full_Rpl3.pdf",height=9, width=11) + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
plotTSNE(example_sce, colour_by = "Rpl3")
dev.off()

pdf("TSNE_w0_Rpl3.pdf",height=9, width=11)
plotTSNE(example_sce_w0, colour_by = "Rpl3") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w2_Rpl3.pdf",height=9, width=11)
plotTSNE(example_sce_w2, colour_by = "Rpl3") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w5_Rpl3.pdf",height=9, width=11)
plotTSNE(example_sce_w5, colour_by = "Rpl3") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w8_Rpl3.pdf",height=9, width=11)
plotTSNE(example_sce_w8, colour_by = "Rpl3") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()

pdf("TSNE_w11_Rpl3.pdf",height=9, width=11)
plotTSNE(example_sce_w11, colour_by = "Rpl3") + scale_color_gradient2(low = "#B5B2B2", mid = "#D63334", high = "#900000",midpoint = 3, limits = c(0, 6))
dev.off()
