#!/usr/bin/env/Rscript 
## pipeline for pairwise alignment of two genomes using lastz, and UCSC's kent tools. 
## Assumes lastz is in system path, as well as kent binaries. 
##
## Input: Target and Query 2bit files
## Output: Pairwise genome alignment of known chromosomes

## lastz 
## kent binaries download instructions (http://genomewiki.ucsc.edu/index.php/DoBlastzChainNet.pl#Parasol_Job_Control_System)

## From CNEr package (https://bioconductor.riken.jp/packages/3.7/bioc/vignettes/CNEr/inst/doc/PairwiseWholeGenomeAlignment.html)

library(CNEr)
library(rtracklayer)
library(optparse)
library(readr)
library(stringr)

# command line options for multiple alignment pipeline
option_list = list(
  make_option(c("-t", "--target"),
              type="character",
              default=NULL,
              help="Target genome 2bit file (absolute path).",
              metavar = "character"),
  make_option(c("-q", "--query"),
              type="character",
              default = NULL,
              help="Query genome 2bit file (absolute path).",
              metavar = "character"),
  make_option(c("-T", "--targetSizes"),
              type="character",
              help="Target chrom.sizes file specifies chroms to include. Uses known chroms from twoBit file if not given.",
              default=NULL,
              metavar = "character"),
  make_option(c("-Q", "--querySizes"),
              type="character",
              default=NULL,
              help="Query chrom.sizes file specifies chroms to include. Uses known chroms from twoBit file if not given.",
              metavar = "character"),
  make_option(c("-d", "--outdir"),
              type="character",
              default = NULL,
              help = "Output directory (absolute path).",
              metavar = "character"),
  make_option(c("-p", "--threads"),
              type="numeric",
              default = NULL,
              help = "Number of threads to use for lastz",
              metavar = "numeric"))

# fetches known (cannonical chromosomes) 
# from twoBitFile
getKnownChroms <- function(path) {
  file <- file.path(path)
  chroms <- seqnames(seqinfo(TwoBitFile(file)))  
  known <- chroms[!grepl("_", chroms)]
  str_sort(known, numeric=T)
  return(known)
}

opt_parser = OptionParser(option_list = option_list)
opt = parse_args(opt_parser)

## check command line options
if (is.null(opt$query) | is.null(opt$target) | is.null(opt$outdir)){
  print_help(opt_parser)
  stop("Error, must supply all command line arguments.", call. = FALSE)
}

## create output directory if does note exits
axtDir <- opt$outdir
if (!dir.exists(axtDir)) {
  dir.create(axtDir, recursive = TRUE)
}
## target genome 
assemblyTarget <- file.path(opt$target)
if (!(file.exists(assemblyTarget))) {
  stop("Error, target genome .2bit file not found!")
}
print(paste0("Using target assembly: ",assemblyTarget))

## query genome
assemblyQuery <- file.path(opt$query)
if (!(file.exists(assemblyQuery))) {
  stop("Error, query genome .2bit file not found!")
}
print(paste0("Using query assembly: ", assemblyQuery))

## use known chroms if target sizes not specified
if (is.null(opt$targetSizes)) {
  ## extract known chroms target
  chromsTargetKnown <- getKnownChroms(opt$target)
  print(paste0("Using known target chroms:", paste(chromsTargetKnown, collapse = " ")))
} else {
  ## use opt$tSizes
  tSizes <- file.path(opt$targetSizes)
  if (!file.exists(tSizes)) {
    stop("tSizes file does not exist!")
  }
  chromsTargetKnown <- read_delim(tSizes, delim = "\t", col_names = F)$X1
  print(paste0("Using target chroms from sizes file: ", paste(chromsTargetKnown, collapse = " ")))
}

## use known chroms if query sizes not specified
if (is.null(opt$querySizes)) {
  ## extract known chroms query
  chromsQueryKnown <- getKnownChroms(opt$query)
  print(paste0("Using known query chroms: ", paste(chromsQueryKnown, collapse = " ")))
} else  {
 ## use opt qSizes 
 qSizes <- file.path(opt$querySizes)
 print(qSizes)
 if (!file.exists(qSizes)) {
   stop("qSizes file does not exist!")
 }
 chromsQueryKnown <- read_delim(qSizes, delim = "\t", col_names = F)$X1
 print(chromsQueryKnown)
 print(paste0("Using query chroms from sizes file: ", paste(chromsQueryKnown, collapse = " ")))
}

if (is.null(opt$threads)) {
  cores <- 2L
} else {
  cores <- opt$threads 
}

# run lastz
lavs <- lastz(assemblyTarget, assemblyQuery, 
              chrsTarget=chromsTargetKnown,
              chrsQuery=chromsQueryKnown,
              outputDir=axtDir,
              distance="near", mc.cores=cores, binary = 'lastz-1.04.00')

## lav files to psl files conversion
psls <- lavToPsl(lavs, removeLav=TRUE, binary="lavToPsl")


## Join close alignments
chains <- axtChain(psls, assemblyTarget=assemblyTarget,
                   assemblyQuery=assemblyQuery, distance="near",
                   removePsl=TRUE, binary="axtChain")


## Sort and combine
allChain <- chainMergeSort(chains, assemblyTarget, assemblyQuery,
                           allChain=file.path(axtDir,
                                              paste0(sub("\\.2bit$", "", basename(assemblyTarget),
                                                         ignore.case=TRUE), ".", 
                                                     sub("\\.2bit$", "", basename(assemblyQuery), 
                                                         ignore.case=TRUE), ".all.chain")),
                           removeChains=FALSE, binary="chainMergeSort")


## Filtering out chains
allPreChain <- chainPreNet(allChain, assemblyTarget, assemblyQuery,
                           allPreChain=file.path(axtDir,
                                                 paste0(sub("\\.2bit$", "", 
                                                            basename(assemblyTarget),
                                                            ignore.case = TRUE), ".", 
                                                        sub("\\.2bit$", "",
                                                            basename(assemblyQuery),
                                                            ignore.case = TRUE),
                                                        ".all.pre.chain")),
                           removeAllChain=FALSE, binary="chainPreNet")

## Keep the best chain and add synteny information
netSyntenicFile <- chainNetSyntenic(allPreChain, assemblyTarget, assemblyQuery,
                                    netSyntenicFile=file.path(axtDir,
                                                              paste0(sub("\\.2bit$", "",
                                                                         basename(assemblyTarget),
                                                                         ignore.case = TRUE), ".",
                                                                     sub("\\.2bit$", "",
                                                                         basename(assemblyQuery),
                                                                         ignore.case = TRUE),
                                                                     ".noClass.net")),
                                    binaryChainNet="chainNet", binaryNetSyntenic="netSyntenic")


## create net.axt file
netToAxt(netSyntenicFile, allPreChain, assemblyTarget, assemblyQuery,
         axtFile=file.path(axtDir,
                           paste0(sub("\\.2bit$", "",
                                      basename(assemblyTarget),
                                      ignore.case = TRUE), ".",
                                  sub("\\.2bit$", "",
                                      basename(assemblyQuery),
                                      ignore.case = TRUE),
                                  ".net.axt")),
         removeFiles=FALSE,
         binaryNetToAxt="netToAxt", binaryAxtSort="axtSort")
