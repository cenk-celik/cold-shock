rnaLocal <- getwd()

##----Install libraries----
libraries <- rlang::quos(Rsubread)
invisible(lapply(lapply(libraries, rlang::quo_name),
                 library,
                 character.only = TRUE
))

##----Set the directory for reference genome & annotation file
refLocal <- paste0(rnaLocal,"ref/") 
dir.create(refLocal, recursive = TRUE, showWarnings = FALSE)

##----Download genome and annotation files----
faFTP <- "ftp://ftp.ensembl.org/pub/release-104/fasta/caenorhabditis_elegans/dna/" # link for reference genome
gtfFTP <- "ftp://ftp.ensembl.org/pub/release-104/gtf/caenorhabditis_elegans/" # link for annotation file

faFile <- "Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.gz" # name of reference genome
gtfFile <- "Caenorhabditis_elegans.WBcel235.104.gtf.gz" # name of annotation file

download.file(paste0(faFTP, faFile), paste0(refLocal, faFile)) # download reference genome
download.file(paste0(gtfFTP, gtfFile), paste0(refLocal, gtfFile)) # download annotation

##----Build the index----
idxLocal <- paste0(rnaLocal,"index/") # create index using reference genome
dir.create(idxLocal, recursive= TRUE, showWarnings = FALSE) # directory for index

system(paste0("gunzip ",paste0(refLocal, faFile))) # unzip reference genome archive
refFile <- list.files(refLocal, pattern = "fa$", full = TRUE) # find reference genome file
buildindex(basename = paste0(idxLocal,"Ce104"), reference = refFile) # unzip with the name

##----Build Index----
mcIndex <- paste0(idxLocal,"ce104") # name of index
setwd(rnaLocal)

bamLocal <- paste0(rnaLocal,"bam/") # directory to save aligned samples
dir.create(bamLocal, recursive = TRUE, showWarnings = FALSE) # create directory for .BAM

mergedLocal <- paste0(rnaLocal, "merged")
dir.create(mergedLocal, recursive = TRUE, showWarnings = FALSE)

setwd(mergedLocal)

folderList <- list.dirs(path = dataLocal, full.names = TRUE, recursive = FALSE)
dataList <- list.dirs(path = dataLocal, full.names = FALSE, recursive = FALSE)

##----Align----
lapply(dataList, function(x){ # function to align all raw data # replace folderList with dataList
  setwd(x) # select each folder in merged/ one by one
  read1 <- list.files(getwd(), pattern = "*_R1_cat.fastq.gz*", full.names = TRUE) # get read one
  fileStub <- gsub("_lib.*", "", basename(x)) 
  align(index = mcIndex,
        readfile1 = read1, # 5' read
        output_file = paste0(bamLocal, fileStub,".bam")) # rename the .BAM files
  setwd(mergedLocal)
})

##----Prepare Annotation File----
system(paste0("gunzip ", paste0(refLocal, gtfFile)))
gtf <- list.files(refLocal, pattern=".*gtf$", full.names = TRUE)

##----Feature Counts----
bamFiles <- list.files(bamLocal, pattern = ".*bam$")
setwd(bamLocal) #temporarily set directory to .BAM folder

featurecounts <- featureCounts(files = bamFiles, GTF.featureType = "exon",
                               GTF.attrType = "gene_id", annot.ext = gtf,
                               isGTFAnnotationFile = TRUE, isPairedEnd = FALSE)

setwd(rnaLocal)

knitr::kable(featurecounts$stat)

##----Save featureCounts table----
dir.create(paste0(rnaLocal, "Rdata"), recursive = FALSE, showWarnings = FALSE)
save(featurecounts, file = paste0(rnaLocal, "Rdata/featurecounts2.RData"))

##----Session info----
sessionInfo()