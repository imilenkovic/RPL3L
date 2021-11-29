library("devtools")
install_github(repo = "lcalviell/Ribo-seQC")
library("RiboseQC")

suppressPackageStartupMessages(library("RiboseQC"))


prepare_annotation_files(annotation_directory = ".",
                         twobit_file = "~/references/GRCm39.genome.2bit",
                         gtf_file = "/users/enovoa/imilenkovic/RiboseQC/gencode.vM27.annotation.gtf", scientific_name = "Mus.musculus",
                         annotation_name = "GRCm38",export_bed_tables_TxDb = F,forge_BSgenome = T,create_TxDb = T)

load_annotation("gencode.vM27.annotation.gtf_Rannot")
getSeq(genome_seq,GTF_annotation$cds_txs[[4]])

GTF_annotation$start_stop_codons
GTF_annotation$cds_txs_coords

GTF_annotation$genome_package

RiboseQC_analysis(annotation_file="gencode.vM27.annotation.gtf_Rannot",bam_files = c("/users/enovoa/imilenkovic/RiboseQC2/star/bams/unaligned/KO1Aligned.sortedByCoord.out.bam", "/users/enovoa/imilenkovic/RiboseQC2/star/bams/unaligned/KO2Aligned.sortedByCoord.out.bam", "/users/enovoa/imilenkovic/RiboseQC2/star/bams/unaligned/KO3Aligned.sortedByCoord.out.bam", "/users/enovoa/imilenkovic/RiboseQC2/star/bams/unaligned/WT1Aligned.sortedByCoord.out.bam", "/users/enovoa/imilenkovic/RiboseQC2/star/bams/unaligned/WT2Aligned.sortedByCoord.out.bam", "/users/enovoa/imilenkovic/RiboseQC2/star/bams/unaligned/WT3Aligned.sortedByCoord.out.bam", report_file = "All_RiboseQC_riboseq.html",write_tmp_files = F) #, readlength_choice_method = "all")
