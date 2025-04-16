#' @import GenomicRanges
#' @import gUtils
#' @import rtracklayer
#' @import GenomeInfoDb
#' @importFrom GenomeInfoDb seqlengths
#' @importFrom data.table data.table fread rbindlist set setkey setkeyv setnames transpose as.data.table
#' @importFrom stats cor loess predict quantile
#' @importFrom skidb read_gencode
#' @importFrom Biostrings alphabetFrequency
#' @importFrom parallel mclapply
#' @importFrom Rsamtools BamFile
#' @importFrom grDevices col2rgb rgb png dev.off
#' @importFrom IRanges IRanges
#' @importFrom DNAcopy CNA segment smooth.CNA
#' @importFrom graphics plot lines
#' @importFrom utils globalVariables
globalVariables(c(".", ".N", ".SD", ":=", "V1", "V2", "V3", "bam", "bin", "black_list_pct", "blacklisted", "child", "chr", "count", "counts", "gc.wig", "index", "ix.start", "lev0", "lev1", "map.wig", "newcount", "numlevs", "pair", "parent", "reads", "reads.corrected", "rowid", "sortSeqlevels", "tmp"))
"_PACKAGE"