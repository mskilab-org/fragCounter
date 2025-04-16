#!/usr/bin/env Rscript
library(optparse)

dr.str = "

███████╗██████╗  █████╗  ██████╗  ██████╗ ██████╗ ██╗   ██╗███╗   ██╗████████╗███████╗██████╗
██╔════╝██╔══██╗██╔══██╗██╔════╝ ██╔════╝██╔═══██╗██║   ██║████╗  ██║╚══██╔══╝██╔════╝██╔══██╗
█████╗  ██████╔╝███████║██║  ███╗██║     ██║   ██║██║   ██║██╔██╗ ██║   ██║   █████╗  ██████╔╝
██╔══╝  ██╔══██╗██╔══██║██║   ██║██║     ██║   ██║██║   ██║██║╚██╗██║   ██║   ██╔══╝  ██╔══██╗
██║     ██║  ██║██║  ██║╚██████╔╝╚██████╗╚██████╔╝╚██████╔╝██║ ╚████║   ██║   ███████╗██║  ██║
╚═╝     ╚═╝  ╚═╝╚═╝  ╚═╝ ╚═════╝  ╚═════╝ ╚═════╝  ╚═════╝ ╚═╝  ╚═══╝   ╚═╝   ╚══════╝╚═╝  ╚═╝


"


library(optparse)
options(bitmapType='cairo')

if (!exists('opt'))
    {
      option_list = list(
        make_option(c("-b", "--bam"), type = "character", help = "Path to .bam or .cram file"),
        make_option(c("-c", "--cov"), type = "character", help = "Path to existing coverage rds or bedgraph"),
        make_option(c("-m", "--midpoint"), type = "character", default = "TRUE", help = "If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval"),
        make_option(c("-w", "--window"), type = "integer", default = 1000, help = "Window / bin size"),
        make_option(c("-d", "--gcmapdir"), type = "character", help = "Mappability / GC content dir"),
        make_option(c("-q", "--minmapq"), type = "integer", default = 20, help = "Minimal map quality"),
        make_option(c("-r", "--reference"), type = "character", default = NULL, help = "Reference FASTA path (required for CRAM)"),
        make_option(c("-p", "--paired"), type = "logical", default = TRUE, help = "Is paired"),
        make_option(c("-e", "--exome"), type = "logical", default = FALSE, help = "Use exons as bins instead of fixed window"),
        make_option(c("-u", "--use.skel"), type = "logical", default = FALSE, help = "Use user defined regions instead of default exome skeleton"),
        make_option(c("-s", "--skeleton"), type = "character", help = "Path to skeleton file"),
        make_option(c("-o", "--outdir"), type = "character", default = './', help = "Directory to dump output into"),
        make_option("--require_flags", type = "character", default = '2', help = "-f flag for samtools used in fragCounter"),
        make_option("--exclude_flags", type = "character", default = '3868', help = "-f flag for samtools used in fragCounter"),
        make_option(c("-l", "--libdir"), type = "character", default = paste(Sys.getenv('GIT_HOME'), 'isva', sep = '/'), help = "Directory containing this R file")
      )

      parseobj = OptionParser(option_list=option_list)
      opt = parse_args(parseobj)

      if (is.null(opt$libdir) | (is.null(opt$bam) & is.null(opt$cov)))
        stop(print_help(parseobj))
      
      ## keep record of run
      writeLines(paste(paste('--', names(opt), ' ', sapply(opt, function(x) paste(x, collapse = ',')), sep = '', collapse = ' '), sep = ''), paste(opt$outdir, 'cmd.args', sep = '/'))
      saveRDS(opt, paste(opt$outdir, 'cmd.args.rds', sep = '/'))
    }


#############################
suppressWarnings(suppressPackageStartupMessages(library(fragCounter)))
suppressWarnings(suppressPackageStartupMessages(library(skidb)))
suppressWarnings(suppressPackageStartupMessages(library(data.table)))
message(dr.str)

# Samtools flags
# 
# Samtools flags for fragCounter bam query
# -F exclude ANY bits set
# 3868 sets bits
# 4 - Read unmapped
# 8 - Mate unmapped
# 16 - Read reverse strand
# 256 - Not primary alignment (secondary alignment)
# 512 - Read fails platform/vendor quality checks
# 1024 - Read is PCR or optical duplicate
# 2048 - Supplementary alignment (chimeric, non-contiguous alignment)
# -f requires read have ALL bits set
# 2 - Read mapped in proper pair

is_null = is.null(opt$require_flags)
is_character = is.character(opt$require_flags)
is_numeric = is.numeric(opt$require_flags)
is_len_one = NROW(opt$require_flags) == 1
is_na = is_len_one && (is.na(opt$require_flags) || (is_character && opt$require_flags %in% "NA"))
is_valid = is_len_one && (is_character || is_numeric)
require_flags = "-f 2"
if (!is_null && !is_na && is_valid) {
  require_flags = paste("-f", opt$require_flags)
}

is_null = is.null(opt$exclude_flags)
is_character = is.character(opt$exclude_flags)
is_numeric = is.numeric(opt$exclude_flags)
is_len_one = NROW(opt$exclude_flags) == 1
is_na = is_len_one && (is.na(opt$exclude_flags) || (is_character && opt$exclude_flags %in% "NA"))
is_valid = is_len_one && (is_character || is_numeric)
exclude_flags = "-F 3868"
if (!is_null && !is_na && is_valid) {
  exclude_flags = paste("-F", opt$exclude_flags)
}

st.flag = paste(require_flags, exclude_flags)

out = fragCounter(bam = opt$bam, cov = opt$cov, window = opt$window, reference = opt$reference, gc.rds.dir = opt$gcmapdir, map.rds.dir = opt$gcmapdir, minmapq = opt$minmapq, paired = opt$paired, outdir = opt$outdir, exome = opt$exome, use.skel = opt$use.skel, skeleton = opt$skeleton, st.flag = st.flag)

saveRDS(out, paste(opt$outdir, 'cov.rds', sep = '/'))
