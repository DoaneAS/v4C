############################################################################
## Ashley Stephen Doane
## Weill Cornell Medicine  asd2007@med.cornell.edu

## This program is free software: you can redistribute it and/or modify it
## under the terms of the GNU Lesser General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## This program is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
## GNU General Public License for more details.

## You should have received a copy of the GNU Lesser General Public License
## along with this program.  If not, see <http://www.gnu.org/licenses/>.
###############################################################################






#' @useDynLib v4C
#' @importFrom Rcpp sourceCpp
#' @importFrom InteractionSet InteractionSet
#' @importFrom data.table as.data.table
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name gi.straw
#' @title gi.straw
#' @description
#' queries .hic object via straw API https://github.com/theaidenlab/straw/tree/master/R
#' to return interactionSet
#'
#'
#'
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC", "VC_SQRT")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return interactionSet
#' @export
#' @author Ashley S Doane
gi.straw = function(hic, gr, norm = "KR", type = 'BP', res = 1e4, mc.cores = 1, ...)
{
  hic = normalizePath(hic)

  gr = reduce(gr.stripstrand(gr[, c()]))
  gr <- gr.nochr(gr)
  #gr <- diffloop::rmchr(gr)

  grs = as.data.table(gr)[, paste(seqnames, start, end, sep = ":")]
  n = length(gr)
  grs.pairs = chr2.pairs = chr1.pairs = NULL
  grs.singles = paste(grs, grs)
  chr1.singles = chr2.singles = as.character(seqnames(gr))
  if (n>1)
  {
    pairs = t(combn(n, 2)) ##pairs
    grs.pairs = paste(grs[pairs[,1]], grs[pairs[,2]])   ## combinations
    #' Tuesday, Nov 07, 2017 04:43:52 PM - Julie: adding as.character()
    chr1.pairs = as.character(seqnames(gr)[pairs[,1]])
    chr2.pairs = as.character(seqnames(gr)[pairs[,2]])
  }

  str = paste(norm, hic, c(grs.singles, grs.singles), type, as.integer(res))
  #str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  chr1 = as.character(c(chr1.singles, chr1.pairs))
  chr2 = as.character(c(chr2.singles, chr2.pairs))
  gi1 = straw_R(str)

  a1 <- data.frame(start = gi1$x, end= gi1$x + res -1, chr=chr1.singles)
  a2 <- data.frame(start = gi1$y, end= gi1$y + res -1, chr=chr1.singles)
  a1 <- makeGRangesFromDataFrame(a1)
  a2 <- makeGRangesFromDataFrame(a2)
  a1 <- diffloop::addchr(a1)
  a2 <- diffloop::addchr(a2)


  gi <- GInteractions(a1, a2)
  mcols(gi)$cts <- gi1$counts
  colnames(mcols(gi)) <- basename(hic)
  iset <- InteractionSet(assays = as.matrix(gi1$counts), interactions = gi)



  return(gi)
  # out = rbindlist(mcmapply(str = str, chr1 = chr1, chr2 = chr2,
  #                          FUN = function(str, chr1, chr2)
  #                          {
  #                            dt = as.data.table(.Call("_gTrack_straw_R", str))
  #                            dt[, ":="(chr1 = chr1, chr2 = chr2)]
  #                            return(dt)
  #                          }, SIMPLIFY = FALSE,  mc.cores = mc.cores))
  #

}









#' @useDynLib v4C
#' @importFrom Rcpp sourceCpp
#' @importFrom InteractionSet InteractionSet
#' @importFrom data.table as.data.table
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name getHic
#' @title getHic
#' @description
#' queries .hic object via straw API https://github.com/theaidenlab/straw/tree/master/R
#' to return sparse upper triangle matrix
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC", "VC_SQRT")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return interactionSet
#' @export
#' @author Ashley S Doane
getHic = function(hic, gr, norm = "KR", type = 'BP', res = 1e4, mc.cores = 1, ...)
{
  hic = normalizePath(hic)

  gr = GenomicRanges::reduce(gr.stripstrand(gr[, c()]))
  gr <- gr.nochr(gr)
  #gr <- diffloop::rmchr(gr)

  grs = as.data.table(gr)[, paste(seqnames, start, end, sep = ":")]
  n = length(gr)
  grs.pairs = chr2.pairs = chr1.pairs = NULL
  grs.singles = paste(grs, grs)
  chr1.singles = chr2.singles = as.character(seqnames(gr))
  #' if (n>1)
  #' {
  #'   pairs = t(combn(n, 2)) ##pairs
  #'   grs.pairs = paste(grs[pairs[,1]], grs[pairs[,2]])   ## combinations
  #'   #' Tuesday, Nov 07, 2017 04:43:52 PM - Julie: adding as.character()
  #'   chr1.pairs = as.character(seqnames(gr)[pairs[,1]])
  #'   chr2.pairs = as.character(seqnames(gr)[pairs[,2]])
  #' }

  str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  #str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  chr1 = as.character(c(chr1.singles, chr1.pairs))
  chr2 = as.character(c(chr2.singles, chr2.pairs))
  out = data.table::rbindlist(mcmapply(str = str, chr1 = chr1, chr2 = chr2,
                                       FUN = function(str, chr1, chr2)
                                       {
                                         dt = data.table::as.data.table(
                                           tryCatch(v4C:::straw_R(str),
                                                    error = function (e)
                                                      data.table()))
                                         return(dt)
                                       }, SIMPLIFY = FALSE,  mc.cores = mc.cores), fill = TRUE)
  out = out[!is.na(counts), ]
  if (!nrow(out))
    stop('Query resulted in no output, please check .hic file or input coordinates')
  out[, end1 := start1+res-1]
  out[, end2 := start2+res-1]
  #out[, str1 := paste0(chr1, ':', start1, '-', end1)]
  #out[, str2 := paste0(chr2, ':', start2, '-', end2)]
  out <- out[, c("chr1", "start1", "end1", "chr2", "start2", "end2", "counts")]
  out <- unique(out)
  a1 <- data.frame(start = out$start1, end= out$end1, chr=out$chr1)
  a2 <- data.frame(start = out$start2, end= out$end2, chr=out$chr2)
  a1 <- makeGRangesFromDataFrame(a1, seqinfo = gr.nochr(hg38s))
  a2 <- makeGRangesFromDataFrame(a2, seqinfo = gr.nochr(hg38s))

  gi <- GInteractions(a1, a2)

  gi$counts <- out$counts
  #gi1 = v4C:::straw_R(str)
  #dfout <- cbind.data.frame(chr= paste0("chr",chr1.singles), region1=gi1$start1, region2=gi1$start2, IF=gi1$counts)
  return(gi)
}




#' @useDynLib v4C
#' @importFrom Rcpp sourceCpp
#' @importFrom InteractionSet InteractionSet
#' @importFrom data.table as.data.table
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name getHicSparse
#' @title getHic
#' @description
#' queries .hic object via straw API https://github.com/theaidenlab/straw/tree/master/R
#' to return sparse upper triangle matrix
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC", "VC_SQRT")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return interactionSet
#' @export
#' @author Ashley S Doane
getHicSparse = function(hic, gr, norm = "KR", type = 'BP', res = 1e4, mc.cores = 1, ...)
{
  hic = normalizePath(hic)

  gr = GenomicRanges::reduce(gr.stripstrand(gr[, c()]))
  gr <- gr.nochr(gr)
  #gr <- diffloop::rmchr(gr)

  grs = as.data.table(gr)[, paste(seqnames, start, end, sep = ":")]
  n = length(gr)
  grs.pairs = chr2.pairs = chr1.pairs = NULL
  grs.singles = paste(grs, grs)
  chr1.singles = chr2.singles = as.character(seqnames(gr))
  #' if (n>1)
  #' {
  #'   pairs = t(combn(n, 2)) ##pairs
  #'   grs.pairs = paste(grs[pairs[,1]], grs[pairs[,2]])   ## combinations
  #'   #' Tuesday, Nov 07, 2017 04:43:52 PM - Julie: adding as.character()
  #'   chr1.pairs = as.character(seqnames(gr)[pairs[,1]])
  #'   chr2.pairs = as.character(seqnames(gr)[pairs[,2]])
  #' }

  str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  #str = paste(norm, hic, c(grs.singles, grs.pairs), type, as.integer(res))
  chr1 = as.character(c(chr1.singles, chr1.pairs))
  chr2 = as.character(c(chr2.singles, chr2.pairs))
  out = data.table::rbindlist(mcmapply(str = str, chr1 = chr1, chr2 = chr2,
                                       FUN = function(str, chr1, chr2)
                                       {
                                         dt = data.table::as.data.table(
                                           tryCatch(v4C:::straw_R(str),
                                                    error = function (e)
                                                      data.table()))
                                         return(dt)
                                       }, SIMPLIFY = FALSE,  mc.cores = mc.cores), fill = TRUE)
  out = out[!is.na(counts), ]
  if (!nrow(out))
    stop('Query resulted in no output, please check .hic file or input coordinates')

  out[, `:=`(end1, start1 + res - 1)]
  out[, `:=`(end2, start2 + res - 1)]
  out <- out[, c("chr1", "start1", "end1", "chr2", "start2",
                 "end2", "counts")]
  out <- unique(out)
  return(out)
}



gr.stripstrand = function(gr)
{
  GenomicRanges::strand(gr) = "*"
  return(gr)
}





#' @useDynLib v4C
#' @importFrom Rcpp sourceCpp
#' @importFrom InteractionSet InteractionSet
#' @importFrom data.table as.data.table
#' @importFrom data.table rbindlist
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name getViewPt
#' @title getViewPt
#' @description
#' queries .hic object via straw API https://github.com/theaidenlab/straw/tree/master/R
#' to return virtual 4C data
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query representing viewpoint
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC", "VC_SQRT")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param fwd number of bins downstream of viewpoint
#' @param rev number of bins upstream of viewpoint
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return interactionSet
#' @export
#' @author Ashley S Doane
getViewPt = function(hic, gr, norm = "KR", type = 'BP', res = 5e3, fwd=20, rev=20, mc.cores = 1, ...)
{
  hic = normalizePath(hic)

  binsize = res
  fwds = binsize * fwd
  revs = binsize* rev
  grf <- resize(gr, width = fwds, fix = "start")
  grr <- resize(gr, width = revs, fix = "end")
  sinfo = seqinfo(gr)
  gr <- resize(gr, width = res, fix="center")

  windows = unlist(c(tile(grr, width = binsize), tile(grf, width = binsize)))
  #windows <- GenomicRanges::reduce(windows)
  windows = gr.stripstrand(windows[,c()])
  windows <- diffloop::rmchr(windows)
  wgrs = as.data.table(windows)[, paste(seqnames, start, end, sep = ":")]
  n = length(windows)

  grx = GenomicRanges::reduce(gr.stripstrand(gr[, c()]))
  grx <- diffloop::rmchr(grx)

  grs = as.data.table(grx)[, paste(seqnames, start, end, sep = ":")]
  #n = length(grx)
  grs.pairs = chr2.pairs = chr1.pairs = NULL
  grs.singles = paste( grs,  wgrs)
  chr1.singles = chr2.singles = as.character(seqnames(windows))


  #str = paste(norm, hic, c(grs.singles, grs.singles), type, as.integer(res))
  str = paste(norm, hic, grs.singles, type, as.integer(res))
  chr1 = as.character(c(chr1.singles))
  chr2 = as.character(c(chr2.singles))
  out = rbindlist(mcmapply(str = str, chr1 = chr1, chr2 = chr2,
                           FUN = function(str, chr1, chr2)
                           {
                             dt = as.data.table(.Call("_v4C_straw_R", str))
                             dt[, ":="(chr1 = chr1, chr2 = chr2)]
                             return(dt)
                           }, SIMPLIFY = FALSE,  mc.cores = mc.cores))

  out[, x2 := x+res-1]
  out[, y2 := y+res-1]
  out[, str1 := paste0(chr1, ':', x, '-', x2)]
  out[, str2:= paste0(chr2, ':', y, '-', y2)]
  gr.out = unique(dt2gr(rbind(out[, .(seqnames = chr1, start = x, end = x2)],
                              out[, .(seqnames = chr2, start = y, end = y2)])))
  gr.out$grs = gr.string(gr.out)
  a1 = dt2gr(out[, .(seqnames = chr1, start = x, end = x2, seqinfo = sinfo)])
  a2 = dt2gr(out[, .(seqnames = chr2, start = y, end = y2, seqinfo=sinfo)])
  a1 = diffloop::addchr(a1)
  a2 = diffloop::addchr(a2)

  gi = InteractionSet::GInteractions(anchor1 = a1, anchor2 = a2)
  mcols(gi)$counts <- out$counts
  # out[, i1 := match(str1, gr.out$grs)]
  # out[, j1 := match(str2, gr.out$grs)]
  # out[, i := pmin(i1, j1)]
  # out[, j := pmax(i1, j1)]
  # out$chr  = paste0("chr", out$chr1)
  # a1 <- makeGRangesFromDataFrame(out,seqnames.field = "chr",start.field = "x", end.field = "y",ignore.strand = TRUE,
  #                                 starts.in.df.are.0based = TRUE)
  # a2 <- makeGRangesFromDataFrame(a2, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  #
  # gi = GInteractions(a1, a2 )
  # dfout <- cbind.data.frame(chr= paste0("chr",chr1.singles), region1=gi1$x, region2=gi1$y, IF=gi1$counts)
  return(gi)
}




#' @useDynLib v4C
#' @importFrom Rcpp sourceCpp
#' @importFrom InteractionSet InteractionSet
#' @importFrom data.table as.data.table
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name getmultiHic
#' @title getHic
#' @description
#' queries .hic object via straw API https://github.com/theaidenlab/straw/tree/master/R
#' to return sparse upper triangle matrix
#' @keywords straw
#' @param hic path to .hic file
#' @param gr granges to query
#' @param norm string specifying normalization to apply ("KR" (default), "NONE", VC", "VC_SQRT")
#' @param res resolution to query (default 10kb)
#' @param type "BP" vs "FRAG"
#' @param mc.cores parallelization to apply for query (useful when length(gr)>2)
#' @return interactionSet
#' @export
#' @author Ashley S Doane
getmultiHic = function(hic, grs, norm = "VC_SQRT", type = 'BP', res = 1e4, mc.cores = 1, ...)
{
  dfout = rbindlist(mcmapply(gr = grs,
                              FUN = function(gr=grs)
                              {
                                dt <- getHic(hic, gr = gr, norm=norm, res=res)

                                return(dt)
                              }, SIMPLIFY = FALSE,  mc.cores = mc.cores))

  return(dfout)
}


#' @importFrom InteractionSet InteractionSet
#' @importFrom InteractionSet GInteractions
#' @importFrom data.table as.data.table
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name readDeLoop
#' @title getHic
#' @description
#' reads cLoops differential loops
#' to return GInteractionSet
#' @keywords straw
#' @param deLoop path to deLoop file
#' @return GInteractions object
#' @export
#' @author Ashley S Doane
readDeLoop <- function(deLoop="/mnt/athena/projectshg38/datasets/hichip/cloop/H3K27acControl.deloop") {
  a <- read_tsv(deLoop)
  colnames(a) <- make.names(colnames(a))
  a$chr <- str_split_fixed( a$iva, ':', 2)[,1]
  str_split_fixed( a$iva, ':', 2)[,1]
  a$start <- str_split_fixed(
    str_split_fixed( a$iva, ':', 2)[,2], '-', 2)[,1]

  a$end <- str_split_fixed( a$iva, '-', 2)[,2]
  a1 <- a
  a1 <- a1 %>% dplyr::select(chr, start, end, loopId)

  a$chr <- str_split_fixed( a$ivb, ':', 2)[,1]
  str_split_fixed( a$ivb, ':', 2)[,1]
  a$start <- str_split_fixed(
    str_split_fixed( a$ivb, ':', 2)[,2], '-', 2)[,1]

  a$end <- str_split_fixed( a$ivb, '-', 2)[,2]
  a2 <- a
  a2 <- a %>% dplyr::select(chr, start, end, loopId)

  a1 <- makeGRangesFromDataFrame(a1, keep.extra.columns = TRUE)
  a2 <- makeGRangesFromDataFrame(a2, keep.extra.columns = TRUE)

  gi = GInteractions(a1, a2 )
  #ax <- a %>% dplyr::select(loopId:`poisson_p-value_corrected`)
  ax <- a %>% dplyr::select(- c(start, end, chr))

  mcols(gi) <- ax
  return(gi)
}




#' @importFrom InteractionSet InteractionSet
#' @importFrom InteractionSet GInteractions
#' @importFrom data.table as.data.table
#' @importFrom diffloop rmchr
#' @importFrom diffloop addchr
#' @name readDeLoop
#' @title getHic
#' @description
#' reads cLoops differential loops
#' to return GInteractionSet
#' @keywords straw
#' @param deLoop path to deLoop file
#' @return GInteractions object
#' @export
#' @author Ashley S Doane
readCLoop <- function(deLoop="/mnt/athena/projectshg38/datasets/hichip/cloop/H3K27acControl.deloop") {
  a <- read_tsv(deLoop)
  colnames(a) <- make.names(colnames(a))
  a <- a %>% dplyr::filter(a$significant==1)
  a$chr <- str_split_fixed( a$iva, ':', 2)[,1]
  str_split_fixed( a$iva, ':', 2)[,1]
  a$start <- str_split_fixed(
    str_split_fixed( a$iva, ':', 2)[,2], '-', 2)[,1]

  a$end <- str_split_fixed( a$iva, '-', 2)[,2]
  a1 <- a
  a1 <- a1 %>% dplyr::select(chr, start, end, loopId)

  a$chr <- str_split_fixed( a$ivb, ':', 2)[,1]
  str_split_fixed( a$ivb, ':', 2)[,1]
  a$start <- str_split_fixed(
    str_split_fixed( a$ivb, ':', 2)[,2], '-', 2)[,1]

  a$end <- str_split_fixed( a$ivb, '-', 2)[,2]
  a2 <- a
  a2 <- a %>% dplyr::select(chr, start, end, loopId)

  a1 <- makeGRangesFromDataFrame(a1, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)
  a2 <- makeGRangesFromDataFrame(a2, keep.extra.columns = TRUE, starts.in.df.are.0based = TRUE)

  gi = GInteractions(a1, a2 )
  #ax <- a %>% dplyr::select(loopId:`poisson_p-value_corrected`)
  ax <- a %>% dplyr::select(- c(start, end, chr))

  mcols(gi) <- ax
  return(gi)
}


