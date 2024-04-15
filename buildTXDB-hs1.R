###################
library(BSgenome.Hsapiens.UCSC.hs1)
#repeat.coords





library(Gviz)
# make plot of methylation files.
# annotation in C:/Users/boelle/nextcloud_SU\pape/TxDb.T2T-CHM13v2.0
library(GenomicFeatures)
txdb <- loadDb("C:/Users/boelle/nextcloud_SU/paper/TxDb.T2T-CHM13v2.0")
GeneRegionTrack(txdb)


# construction of the TxDb file :
library(rtracklayer)

## Takes < 1 min, consumes about 7Gb of RAM
gff <- import("C:/Users/boelle/nextcloud_SU/paper/GCF_009914755.1_T2T-CHM13v2.0_genomic.gff.gz")


library(GenomeInfoDb)
chrominfo <- getChromInfoFromNCBI("T2T-CHM13v2.0")
seqlevels(gff) <- setNames(chrominfo$SequenceName, chrominfo$RefSeqAccn)

seqinfo(gff) <- Seqinfo(genome="T2T-CHM13v2.0")


gff = gff[gff$type!= "pseudogene",]
gff = gff[gff$type!= "miRNA",]
gff = gff[gff$type!= "snoRNA",]
gff = gff[gff$type!= "snRNA",]
gff = gff[gff$type!= "tRNA",]
gff = gff[gff$type!= "V_gene_segment",]
gff = gff[gff$type!= "C_gene_segment",]
gff = gff[gff$type!= "J_gene_segment",]
gff = gff[gff$type!= "D_gene_segment",]


#library(txdbmaker)
txdb <- makeTxDbFromGRanges(gff, taxonomyId=9606)

seqlevelsStyle(txdb) <- "NCBI"

saveDb(txdb, "C:/Users/boelle/nextcloud_SU/paper/TxDb.T2T-CHM13v2.0")

save(repeat.coords, file = "C:/Users/boelle/nextcloud_SU/paper/repeats.rda")