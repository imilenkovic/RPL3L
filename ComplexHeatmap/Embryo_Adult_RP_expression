library(ComplexHeatmap)
library(circlize)
test <- read.csv("C:\\Users\\imilenkovic\\Desktop\\PhD\\Projects\\Embryo\\RNA-Seq datasets\\embryo_ribosomal_proteins_median_for_heatmap_removed_rt_rrp.csv",row.names = 1)
input <- as.matrix(test)

annotation<- read.delim("annotation.txt")
stage<-annotation[,2]
stage2<- as.data.frame(stage)
stage2 <- t(stage2)
stage2<- as.data.frame(stage2)

row_ha<- rowAnnotation(df=stage2)
column_ha <- columnAnnotation(df=stage2)


tissue<- colnames(input)
tissue2<- gsub("\\..*","", tissue)

input_l<-log(input+1)


pdf("RPs_full_unscaled.pdf",height=7,width=6)
Heatmap(input_l, name = "log(RPKM)", 
        #col = colorRamp2(c(-3,0,4), c("cadetblue3","floralwhite", "maroon4"),space = "RGB"), 
        #cluster_rows = TRUE, 
        col = colorRamp2(c(0, 1, 3, 5 ,7), c("#2c7bb6","#abd9e9","floralwhite","#fdae61", "#d7191c"),space = "RGB"),
        #clustering_distance_columns = "kendall",
        cluster_columns = FALSE,
        column_title = "Developmental stages", 
        column_title_gp = gpar(fontsize = 10, fontface = "bold"),
        column_names_gp = gpar(fontsize = 7, fontface = "bold"),
        row_title = "Ribosomal proteins", row_title_rot = 90,
        row_title_gp = gpar(fontsize = 8, fontface = "bold"),
        cluster_rows = TRUE,
        show_row_names = TRUE,
        row_names_gp = gpar(fontsize = 5), #row names size
        #column_order = 1:dim(data4)[2],#Keep the column order, make clustering FALSE for this
        #row_dend_side = "right", #Dendogram on the right side
        #row_order = 1:dim(data4)[1], #Keep the row order, make clustering FALSE for this
        #show_column_dend = TRUE, #
        column_dend_side = "top",
        column_names_side = "bottom",
        column_split = tissue2, #Splitting by Class
        row_gap = unit(0, "mm"), #Gap
        #top_annotation=column_ha,
)
dev.off()
