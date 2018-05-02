#' @import GenomicRanges
#' @import gUtils
#' @import data.table
#' @import rtracklayer
#' @import bamUtils
#' @import ffTrack


#' @title Mappability calculator
#' @description Calculates mappability as fraction of bases in a tile with mappability of 1
#' @name MAP.fun
#' @param win.size Window size, in basepairs, to calculate mappability for (should match window
#' size of your coverage file
#' @param twobitURL URL to twobit genome file. Default is hg19 from UCSC
#' @param bw.path Path to .bigWig mappability file
#' @param twobit.win How many windows of the twobit file to load into memory on each core
#' @param mc.cores How many cores to use
#' @author Marcin Imielinski
#' @export

MAP.fun = function(win.size = 200, twobitURL = '~/DB/UCSC/hg19.2bit', bw.path = '~/DB/UCSC/wgEncodeCrgMapabilityAlign100mer.bigWig', twobit.win = 1e3, mc.cores = 10) {
  tiles = gr.tile(seqlengths(TwoBitFile(twobitURL)),win.size)
  values(tiles) = NULL
  tiles = tiles %Q% (seqnames %in% c(paste0("chr",seq(1,22)),"chrX","chrY")) # removes all chromomes except 1:22 and X/Y
  seqlevels(tiles) = seqlevelsInUse(tiles) # Get rid of empty seqname factor levels
  x = seq(1,(length(tiles) / twobit.win)) 
  map.score = function(x) {
    this.ix = seq(1 + ((x-1) * twobit.win), x * twobit.win)
    out = data.table(mappability = tiles[this.ix] %O% (rtracklayer::import(bw.path, selection = tiles[this.ix]) %Q% (score==1)), index = this.ix)
    return(out)
  }
  map.out = mclapply(x, map.score, mc.cores = mc.cores)
  map.out = rbindlist(map.out)
  if (!is.integer(length(tiles)/twobit.win)){
    edge.num = seq(twobit.win*max(x)+1,length(tiles))
    edge.out = data.table(mappability = tiles[edge.num] %O% (rtracklayer::import(bw.path, selection = tiles[edge.num]) %Q% (score==1)), index = edge.num)
    map.out = rbind(map.out, edge.out)
  }
  setkey(map.out,index)
  tiles$score = map.out$mappability
  tiles = gr.sub(tiles)
  tiles = gr.stripstrand(tiles)
  return(tiles)
}


#' @title GC content calculator
#' @description Calculates GC content across the genome for a given sized window
#' @name GC.fun
#' @param win.size Window size, in basepairs, to calculate GC content for (should match window
#' size of your coverage file
#' @param twobitURL URL to twobit genome file. Default is hg19 from UCSC
#' @param twobit.win How many windows of the twobit file to load into memory on each core
#' @param mc.cores How many cores to use
#' @author Marcin Imielinski
#' @export

GC.fun = function(win.size = 200, twobitURL = '~/DB/UCSC/hg19.2bit', twobit.win = 1e3, mc.cores = 10) {
  tiles = gr.tile(seqlengths(TwoBitFile(twobitURL)),win.size)
  values(tiles) = NULL
  tiles = tiles %Q% (seqnames %in% c(paste0("chr",seq(1,22)),"chrX","chrY")) # removes all chromomes except 1:22 and X/Y
  x = seq(1,(length(tiles) / twobit.win))
  gc.con = function(x) {
    this.ix = seq(1 + ((x-1) * twobit.win), x * twobit.win)
    # print(this.ix[1])
    tmp = alphabetFrequency(get_seq(twobitURL, tiles[this.ix]))
    out = data.table(gc = rowSums(tmp[, c('C', 'G')])/rowSums(tmp), index = this.ix)
    return(out)
  }
  gc.out = mclapply(x, gc.con, mc.cores = mc.cores)
  gc.out = rbindlist(gc.out)
  if (!is.integer(length(tiles)/twobit.win)){
    edge.num = seq(twobit.win*max(x)+1,length(tiles))
    tmp.edge = alphabetFrequency(get_seq(twobitURL, tiles[edge.num]))
    edge.out = data.table(gc = rowSums(tmp.edge[, c('C', 'G')])/rowSums(tmp.edge), index = edge.num)
    gc.out = rbind(gc.out, edge.out)
  }
  setkey(gc.out,index)
  tiles$score = gc.out$gc
  tiles = gr.sub(tiles)
  tiles = gr.stripstrand(tiles)
  return(tiles)
}


#' @name PrepareCov
#' @title PrepareCov
#' @description Load BAM or coverage file and prepare for use in correctcov_stub
#' @author Marcin Imielinski
#' @param bam path to .bam file
#' @param cov Path to existing coverage rds or bedgraph 
#' @param midpoint If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval
#' @param window window / bin size
#' @param minmapq Minimal map quality
#' @param paired wether or not paired
#' @param outdir Directory to dump output into
#' @export

PrepareCov = function(bam, cov = NULL, midpoint = FALSE, window = 200, minmapq = 1, paired = TRUE, outdir = NULL) {
  if (is.null(bam)) {
    bam = ''
  }
  midpoint = grepl("(T)|(TRUE)", midpoint, ignore.case = T)
  if (file.exists(bam)) { # & is.null(cov))
    if (!midpoint) {
      cat("Running without midpoint!!!\n")
    }
    print('Doing it!')
    if (is.null(paired)) {
      paired = TRUE
    }
    if (paired) {
      cov = bam.cov.tile(bam, window = 200, chunksize = 1e6, midpoint = FALSE, min.mapq = 1)  ## counts midpoints of fragments
    }
    else {
      sl = seqlengths(BamFile(bam))
      tiles = gr.tile(sl, window)
      cov = bam.cov.gr(bam, intervals = tiles, isPaired = NA, isProperPair = NA, hasUnmappedMate = NA, chunksize = 1e7)  ## counts midpoints of fragments    # Can we increase chunksize?
      cov$count = cov$records/width(cov)
    }
  }
  else if (!is.null(cov)) {
    cov = readRDS(cov)
  }
  else {
    stop("Can't locate either bam or coverage input file")
  }
  gc()
  cat('Finished acquiring coverage\n')
  return(cov)
  cat('done\n')
}



#' @title Correct coverage
#' @description Prepares GC, mappability, and coverage files for multicoco
#' @name correctcov_stub
#' @param cov.wig wig file of coverage tiles of width W or pointer to rds file
#' of sorted GRanges object or GRanges object
#' @param mappability threshold for mappability score
#' @param samplesize size of sub-sample
#' @param verbose Wether to print log to console
#' @param gc.rds.dir for tiles of width W, will look here for a file named gc{W}.rds in this directory
#' @param map.rds.dir for tiles of width W will look here for a file named map{W}.rds in this directory
#' @author Marcin Imielinski
#' @export

correctcov_stub = function(cov.wig, mappability = 0.9, samplesize = 5e4, verbose = T, gc.rds.dir, map.rds.dir) {
  if (is.character(cov.wig)) {
    if (grepl('(\\.bedgraph$)|(\\.wig$)|(\\.bw$)', cov.wig)) {
      cov = import.ucsc(cov.wig)
    }
    else if (grepl('(\\.rds)', cov.wig)) {
      cov = gr.stripstrand(readRDS(cov.wig))
      names(values(cov))[1] = 'score' ## will take the first value column as the (uncorrected) read count
    }
    else {
      stop("Unsupported coverage format (should be UCSC bedgraph, wig, bw, or R Bioconductor GRanges .rds file")
    }
  }
  else { ## assume it is a sorted GRanges
    cov = gr.stripstrand(cov.wig)
  }
  # cov = cov %Q% (seqnames == 21) #' twalradt Monday, Apr 23, 2018 10:33:19 AM Done for unit testing
  n = length(cov)
  wid = as.numeric(names(sort(-table(width(cov))))[1])
  gc.rds = paste(gc.rds.dir,'/gc', wid, '.rds', sep = '')
  map.rds = paste(map.rds.dir,'/map', wid, '.rds', sep = '')
  cat('Loaded GC and mappability\n')
  if (!file.exists(gc.rds) | !file.exists(map.rds)) {
    stop(sprintf('GC rds file %s not found, either directory is wrong or file needs to be generated for this particular window width', gc.wig))
    stop(sprintf('mappability rds file %s not found, either directory is wrong or file needs to be generated for this particular window width', map.wig))
  }
  else { ## if we have rds files then let's use these to avoid using rtracklayer
    gc = readRDS(gc.rds)
    map = readRDS(map.rds)
  }
  if (is.null(cov$score)) { ## if $score field is empty then just assume that the first column of coverage is the "score" i.e. read count
    names(values(cov))[1] = 'score'
  }
  map = gr.sub(map)
  gc = gr.sub(gc)
  gc.str = gr.string(gc)
  map.str = gr.string(map)
  cov.str = gr.string(cov)
  all.str = intersect(intersect(map.str, cov.str), gc.str) ## in case we have mismatches in the ordering / genome definition       
  cat(sprintf('length cov is %s, length gc is %s, length map is %s\n',
              length(cov),
              length(gc),
              length(gc)))
  map = map[match(all.str, map.str)]
  cov = cov[match(all.str, cov.str)]
  gc = gc[match(all.str, gc.str)]        
  if (length(cov) != length(gc) | length(gc) != length(map)) {
    stop('Mismatch / problem in cov, gc, or map definition.  Check if they come from the same width tiling')
  }  
  cov$reads = cov$score
  cov$gc = gc$score
  cov$gc[cov$gc<0] = NA
  cov$map = map$score
  cov$score = NULL
  rm(gc)
  rm(map)
  gc()  
  cat('Synced coverage, GC, and mappability\n')
  cov = sort(gr.fix(cov))        
  cat('Modified gc / mappability correction\n')
  return(cov)
}


#' @title Multi-scale coverage correction
#' @description Given gc and mappability coverage correction at k "nested" scales finds the coverage
#' assignment at the finest scale that yields the best correction at every scale
#' @name multicoco
#' @param cov constant with GRanges of coverage samples with (by default) fields $reads, $map, $gc
#' @param numlevs numbers of scales at which to correct
#' @param fields fields of gc to use as covariates
#' @param interative whether to iterate
#' @param presegment whether to presegment
#' @param min.segwidth when presegmenting, minimum segment width
#' @param mono Wether to only do single iteration at finest scale
#' @param verbose Wether to print log to console
#' @param FUN function with which to correct coverage (by default loess
#' correction modified from HMMcopy that takes in granges with fields
#' $reads and other fields specified in "fields"
#' @param ... additional args to FUN
#' @author Marcin Imielinski
#' @usage multicoco(cov, numlevs = 1, fiedls = c("gc", "map"), iterative = TRUE,
#' presegment = "TRUE", min.segwidth = 5e6, mono = FALSE, verbose = TRUE,
#' FUN = NULL, mc.cores = 1)
#' @export

multicoco = function(cov, numlevs = 1, base = max(10,1e5/max(width(cov))), fields = c("gc", "map"), iterative = TRUE,
                     presegment = TRUE, min.segwidth = 5e6, mono = FALSE, verbose = T, FUN = NULL, ..., mc.cores = 1) {
  if (verbose) {
    cat('Converting to data.table\n')
  }
  WID = max(width(cov))
  library(data.table)
  cov.dt = gr2dt(cov)        
  sl = structure(as.numeric(1:length(seqlevels(cov))), names = seqlevels(cov))       
  if (verbose) {
    cat('Grouping intervals\n')
  }
  ## compute level means
  ## lev 0 = raw data
  ## lev 1 = base-fold collapsed
  ## lev 2 = base^2-fold collapsed
  ## etc
  parentmap= list() ## data.tables that map lev labels at level k  to parent lev labels
  cov.dt[, lev0:=as.character(1:length(seqnames))]
  for (k in 1:numlevs) {
    if (verbose) {
      cat('Defining', base^k, 'fold collapsed ranges\n')
    }
    cov.dt[, eval(paste("lev", k, sep = '')) := as.character(sl[seqnames] + as.numeric(Rle(as.character(1:length(start)), rep(base^k, length(start)))[1:length(start)])/length(start)), by = seqnames]
    parentmap[[k]] = data.table(parent = cov.dt[, get(paste("lev", k, sep = ''))], child = cov.dt[, get(paste("lev", k-1, sep = ''))], key = 'child')[!duplicated(child), ]
  }
  if (presegment) { ## perform rough segmentation at highest level
    seg = NULL
    sl = seqlengths(cov)
    if (verbose) {
      cat('Presegmenting at ', as.integer(WID*base^(numlevs)), ' bp scale \n')
    }
    require(DNAcopy)
    set.seed(42) #' twalradt Friday, Apr 20, 2018 01:07:28 PM
    tmp.cov = seg2gr(cov.dt[,list(chr = seqnames[1], start = min(start), end = max(end), strand = strand[1], reads = mean(reads, na.rm = T)), by = get(paste("lev", numlevs, sep = ''))][end>start, ], seqlengths = sl)
    ix = which(!is.na(values(tmp.cov)[, 'reads']))
    tmp = data.frame()
    tryCatch({
      cna = CNA(log(values(tmp.cov)[, 'reads'])[ix], as.character(seqnames(tmp.cov))[ix], start(tmp.cov)[ix], data.type = 'logratio')
      tmp = print(segment(smooth.CNA(cna), alpha = 1e-5, verbose = T))
      tmp = tmp[!is.na(tmp$loc.start) & !is.na(tmp$chrom) & !is.na(tmp$loc.end), , drop = F]
    }, error = function(e) warning('DNACopy error moving on without segmenting'))
    if (nrow(tmp)>0) {
      seg = sort(seg2gr(tmp, seqlengths = sl))
      seg = seg[width(seg)>min.segwidth] ## remove small segments
      seg.dt = gr2dt(seg)
      if (nrow(seg.dt)>0) {
        seg = seg2gr(seg.dt[, list(seqnames = seqnames,
                                   start = ifelse(c(FALSE, seqnames[-length(seqnames)]==seqnames[-1]), c(1, start[-1]), 1),
                                   end = ifelse(c(seqnames[-length(seqnames)]==seqnames[-1], FALSE), c(start[-1]-1, Inf), seqlengths(seg)[as.character(seqnames)]))], seqlengths = sl)
        Seg = gr.val(seg, tmp.cov, 'reads') ## populate with mean coverage
        seg$reads = seg$reads/sum(as.numeric(seg$reads*width(seg))/sum(as.numeric(width(seg)))) ## normalize segs by weigthed mean (so these become a correction factor)
      }
      else {
        seg = NULL
      }
    }
    else {
      seg = NULL
    }
  }
  else {
    seg = NULL
  }
  if (verbose) {
    cat('Aggregating coverage within levels \n')
  }
  ## list of data frames showing scales of increasing collapse
  cov.dt[, ix := 1:nrow(cov.dt)]
  cmd1 = paste('list(ix.start = ix[1], ix.end = ix[length(ix)], reads = mean(reads, na.rm = T),', paste(sapply(fields, function(f) sprintf("%s = mean(%s, na.rm = T)", f, f)), collapse = ','), ')', sep = '')
  cmd2 = paste('list(lab = lev0, reads,', paste(fields, collapse = ','), ', seqnames, start, end)', sep = '')
  if (mono) {
    if (verbose) {
      cat('Mono scale correction \n')
    }
    grs = list(cov.dt[, eval(parse(text=cmd2))])
    numlevs = 1
  }
  else {
    grs = c( list(cov.dt[, eval(parse(text=cmd2))]), lapply(1:numlevs, function(x) {
      if (verbose) {
        cat('Aggregating coverage in level', x,  '\n')
      }
      out = cov.dt[, eval(parse(text=cmd1)), keyby = list(lab = get(paste('lev', x, sep = '')))]
      out[, ":="(seqnames = cov.dt$seqnames[ix.start], end = cov.dt$end[ix.start], start = cov.dt$start[ix.start])]
      out[, ":="(ix.start= NULL, ix.end = NULL)]
      return(out)
    }))
  }
  setkey(grs[[1]], 'lab')
  ## modified from HMMCopy to
  ## (1) take arbitrary set of covariates, specified by fields vector
  ## (2) employ as input an optional preliminary (coarse) segmentation with which to adjust signal immediately prior to loess
  ## NOTE: (this only impacts the loess fitting, does not impose any segmentation on the data)
  ##
  if (is.null(FUN)) {
    FUN = function(x, fields = fields, samplesize = 5e4, seg = NULL, verbose = T, doutlier = 0.001, routlier = 0.01) {
      if (!all(fields %in% names(x))) {
        stop(paste('Missing columns:', paste(fields[!(fields %in% names(x))], collapse = ',')))
      }
      x$valid <- TRUE
      for (f in fields) {
        x$valid[is.na(x[, f])] = FALSE
        x$valid[which(is.infinite(x[, f]))] = FALSE
      }
      if (verbose) {
        cat('Quantile filtering response and covariates\n')
      }
      range <- quantile(x$reads[x$valid], prob = c(routlier, 1 - routlier), na.rm = TRUE)
      if (verbose) {
        cat(sprintf("Response min quantile: %s max quantile: %s\n", round(range[1],2), round(range[2],2)))
      }
      domains = lapply(fields, function(f) quantile(x[x$valid, f], prob = c(doutlier, 1 - doutlier), na.rm = TRUE))
      names(domains) = fields
      x$ideal <- x$valid
      x$ideal[x$reads<=range[1] | x$reads>range[2]] = FALSE
      for (f in fields)
      {
        x$ideal[x[, f] < domains[[f]][1] | x[, f] > domains[[f]][2]] = FALSE
      }
      if (verbose) {
        cat(sprintf('Nominated %s of %s data points for loess fitting\n', sum(x$ideal), nrow(x)))
      }
      set <- which(x$ideal)
      if (length(set)<=10) {
        warning("Not enough samples for loess fitting - check to see if missing or truncated data?")
        return(x$reads)
      }
      for (f in fields) {
        if (verbose) {
          message(sprintf("Correcting for %s bias...", f))
        }
        set.seed(42)
        select <- sample(set, min(length(set), samplesize))
        x2 = x[, c('reads', f)]
        x2$covariate = x2[, f]
        x2s = x2[select, ]                    
        if (!is.null(seg)) {  ## here we apply a prelmiinary segmentation to correct for large scale copy variation allow more power to reveal the covariate signal
          if (verbose) {
            message('Applying preliminary segmentation prior to loess fitting')
          }
          x.grs = gr.val(seg2gr(x[select, ], seqlengths = NULL), seg, 'reads')
          x2s$reads = x2s$reads/x.grs$reads
        }
        fit = tryCatch(loess(reads ~ covariate, data = x2s, span = 0.3), error = function(e) NULL)
        x$reads = NA
        if (!is.null(fit)) {
          if (is.na(fit$s)) {
            warning("Using all points since initial loess failed")
            fit = loess(reads ~ covariate, data = x2[select, ], span = 1)
          }
        }
        tryCatch(
        {
          if (!is.na(fit$s)) {
            domain = domains[[f]]
            yrange <- quantile(x2s$reads, prob = c(routlier, 1 - routlier), na.rm = TRUE)
            df = data.frame(covariate = seq(domain[1], domain[2], 0.001))
            x$reads = x2$reads/predict(fit, x2) ## apply correction
          }
          else {
            print("Loess failed, yielding NA loess object, continuing without transforming data")
          }
        }, error = function(e) print("Unspecified loess or figure output error"))
      }
      return(x$reads)
    }
  }
  if (verbose) {
    cat('Correcting coverage at individual scales\n')
  }
  ## level 1,2, ..., k corrections
  ## these are the computed corrected values that will be input into the objective function
  correction = NULL
  for (i in rev(1:length(grs))) {
    cat('Correcting coverage at ', WID*base^(i-1), 'bp scale, with', nrow(grs[[i]]), 'intervals\n')
    if (i != length(grs)) {
      grs[[i]]$reads = grs[[i]]$reads/correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
    }
    if (WID*base^(i-1) > 1e5) { ## for very large intervals do not quantile trim, only remove NA
      grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, doutlier = 0, seg = seg)
    }
    else {
      grs[[i]]$reads.corrected = FUN(as.data.frame(grs[[i]]), fields, seg = seg);
    }
    if (is.null(correction)) {
      correction = data.table(lab = grs[[i]]$lab, cor = grs[[i]]$reads / grs[[i]]$reads.corrected, key = 'lab')
    }
    else {
      ## multiply new correction and old correction
      old.cor = correction[parentmap[[i]][grs[[i]]$lab, parent], cor]
      new.cor = grs[[i]]$reads / grs[[i]]$reads.corrected                                     
      correction = data.table(lab = grs[[i]]$lab,  cor = old.cor * new.cor, key = 'lab') ## relabel with new nodes
    }
  }
  cov.dt$reads.corrected = grs[[1]][cov.dt$lev0, ]$reads.corrected
  rm(grs)
  gc()        
  if (verbose) {
    cat('Converting to GRanges\n')
  }
  gc()      
  out = seg2gr(as.data.frame(cov.dt), seqlengths = seqlengths(cov)) 
  if (verbose) {
    cat('Made GRanges\n')
  }
  gc()
  return(out)
}


#' @name fragCounter
#' @title fragCounter
#' @description Runs entire fragCounter pipeline
#' @author Marcin Imielinski
#' @param bam path to .bam file
#' @param cov Path to existing coverage rds or bedgraph 
#' @param midpoint If TRUE only count midpoint if FALSE then count bin footprint of every fragment interval
#' @param window window / bin size
#' @param minmapq Minimal map quality
#' @param paired wether or not paired
#' @param outdir Directory to dump output into
#' @param gc.rds.dir for tiles of width W, will look here for a file named gc{W}.rds in this directory
#' @param map.rds.dir for tiles of width W will look here for a file named map{W}.rds in this directory

fragCounter = function(bam, cov = NULL, midpoint = FALSE,window = 200, gc.rds.dir, map.rds.dir, minmapq = 1, paired = TRUE, outdir = NULL) {
  cov = PrepareCov(bam, cov = NULL, midpoint = FALSE, window = 200, minmapq = 1, paired = TRUE, outdir)
  cov = correctcov_stub(cov, gc.rds.dir = gc.rds.dir, map.rds.dir = map.rds.dir)
  cov$reads.corrected = multicoco(cov, numlevs = 1, base = max(10, 1e5/window), mc.cores = 1, fields = c('gc', 'map'), iterative = T, mono = T)$reads.corrected
  if (!is.null(outdir)) {
    out.rds = paste(outdir, '/cov.rds', sep = '')
    out.corr = paste(gsub('.rds$', '', out.rds), '.corrected.bw', sep = '')
    if (!is.null(tryCatch({library(rtracklayer); 'success'}, error = function(e) NULL))) {
      cov.corr.out = cov
      cov.corr.out$score = cov$reads.corrected
      cov.corr.out$score[is.na(cov.corr.out$score)] = -1
      cov.corr.out = cov.corr.out[width(cov.corr.out)==window] ## remove any funky widths at end of chromosome
      export(cov.corr.out[, 'score'], out.corr, 'bigWig', dataFormat = 'fixedStep')
    }
    saveRDS(cov, paste(gsub('.rds$', '', out.rds), '.rds', sep = ''))
  }
  return(cov)
  cat('done\n')
}



