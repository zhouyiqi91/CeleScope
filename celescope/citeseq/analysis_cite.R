library(Seurat)
library(tidyverse)
library(argparser)

argv <- arg_parser('')
argv <- add_argument(argv,"--citeseq_mtx", help="citeseq_mtx file")
argv <- add_argument(argv,"--outdir", help="outdir")
argv <- add_argument(argv,"--sample", help="sample")
argv <- add_argument(argv,"--rds", help="rds")
argv <- parse_args(argv)

#args
citeseq_mtx = argv$citeseq_mtx
outdir = argv$outdir
sample = argv$sample
rds = argv$rds

rds = readRDS(rds)
rds.adt <- read.csv(citeseq_mtx, sep = "\t", header = TRUE, row.names = 1)
rds.citeseq <- rds.adt
rownames(rds.citeseq) <- paste0("CITE_", rownames(rds.adt))
rds <- SetAssayData(rds, assay.type = "CITE", slot = "raw.data", new.data = rds.citeseq)
rds <- NormalizeData(rds, assay.type = "CITE", normalization.method = "genesCLR")
rds <- ScaleData(rds, assay.type = "CITE", display.progress = FALSE)
proteins = rownames(rds@assay$CITE@raw.data)
n_proteins = length(proteins)

pdf.out = str_glue('{outdir}/{sample}_CITESeq_featurePlot.pdf')
pdf(pdf.out, height=n_proteins, width=6)
print(FeaturePlot(rds, features.plot = proteins, min.cutoff = "q05", max.cutoff = "q95", 
    nCol = 3, cols.use = c("lightgrey", "blue"), pt.size = 0.5))
dev.off()