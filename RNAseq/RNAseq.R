##########################
## Import Aligned Reads ##
##########################
library(GenomicRanges); library(Rsamtools); 
biocLite("leeBamViews")
library(leeBamViews) # Load required libraries.

testFile <- system.file("bam", "isowt5_13e.bam", package = "leeBamViews")

?readBamGappedAlignments
aligns <- readBamGappedAlignments(testFile) # Imports a BAM alignment file (here yeast example) and stores it as a GappedAlignments object.
class(aligns)

str(aligns)


# aligns2 <- readGAlignmentPairsFromBam(testFile)
# class(aligns2)
# str(aligns2)
# rm(aligns2)

rname(aligns) <- sub("^Sc", "", rname(aligns)); rname(aligns) <- sub("13", "XIII", rname(aligns)) # Required for data consistency.

alignscan <- scanBam(testFile); names(alignscan[[1]]) # Imports BAM data into a nested list object. 

######################################################
## Organize Annotation Data in a GRangesList Object ##
######################################################
library(GenomicFeatures)
txdb <- makeTranscriptDbFromUCSC(genome="sacCer2", tablename="sgdGene") # Creates a TranscriptDb object from transcript annotations available at the UCSC Genome Browser.
exonRanges <- exonsBy(txdb, "tx") # Stores exon data as a GRangesList object.

#######################################################
## Create Annotation Objects from Custom Data Frames ##
#######################################################
transcripts <- data.frame(tx_id=1:3, tx_chrom="chr1", tx_strand=c("-", "+", "+"), tx_start=c(1, 2001, 2001), tx_end=c(999, 2199, 2199))
splicings <- data.frame(tx_id=c(1L, 2L, 2L, 2L, 3L, 3L), exon_rank=c(1, 1, 2, 3, 1, 2), exon_start=c(1, 2001, 2101, 2131, 2001, 2131), exon_end=c(999, 2085, 2144, 2199, 2085, 2199), cds_start=c(1, 2022, 2101, 2131, NA, NA), cds_end=c(999, 2085, 2144, 2193, NA, NA))
myTxdb <- makeTranscriptDb(transcripts, splicings) # Example for creating TranscriptDb object from two data frames.
exonsBy(myTxdb, "tx") # Returns exon data as GRangesList object.

###################################################
## Counting Reads that Overlap Annotation Ranges ##
###################################################
counts <- countOverlaps(exonRanges, aligns) # Counts matching reads per transcript.
numBases <- sum(width(reduce(exonRanges))) # Length of exon union per gene. 
geneLengthsInKB <- (numBases/1000) # Conversion to kbp.
millionsMapped <- sum(counts)/1e+06 # Factor for converting to million of mapped reads. 
rpm <- counts/millionsMapped # RPK: reads per kilobase of exon model.
rpkm <- rpm/geneLengthsInKB # RPKM: reads per kilobase of exon model per million mapped reads. Note: the results are stored in a named vector that matches the index of the initial GRangesList object!!!

sortedRPKM <- sort(rpkm); highScoreGenes <- tail(sortedRPKM)
txs <- transcripts(txdb, vals = list(tx_id = names(highScoreGenes)), columns = c("tx_id", "gene_id"))
systNames <- as.vector(unlist(elementMetadata(txs)["gene_id"])) # Example for returning the six most highly expressed transcripts.

###################################################################
## Counting Reads for Non-Annotation Ranges (Intergenic Regions) ##
###################################################################
filtData <- subsetByOverlaps(aligns, exonRanges) # Returns all reads that align in non-exon ranges. This includes introns which are rare in this example.
filtData2 <- subsetByOverlaps(aligns, transcriptsBy(txdb, "gene")) # Returns all reads that align in intergenic regions (excludes introns).
cov <- coverage(filtData) # Calculates coverage for non-exon ranges.
cov <- cov[13] # Subsetting of chromosome 13 because data for the other chromosomes is missing in this example. 
islands <- slice(cov, lower = 1) # Filters for areas with a continuous read coverage of >=1.
transcribedRegions <- islands[width(islands) > 1000] # Returns only those transcribed regions that are at least 1000 bp long.

#########################################################################
## Avoid Double Counting and Ambiguous Mappings with summarizeOverlaps ##
#########################################################################
## Note: when counting single range features, a GRanges object is passed
## on to summarizeOverlaps and with multiple range features a GRangesList.
## It also accepts a ScanBamParam() function under the param argument to filter 
## the reads stored in a bam file, e.g. by "NH" tag to exclude reads with multiple mappings.  
library(GenomicRanges); library(Rsamtools)
## Create some sample data
rd <- GappedAlignments(letters[1:3], seqnames = Rle(rep("chr1",3)),
      pos = as.integer(c(100, 500, 900)),
      cigar = rep("300M", 3), strand = strand(rep("+", 3)))
gr1 <- GRanges("chr1", IRanges(start=100, width=101), strand="+", ID="gene1_exon1")
gr2 <- GRanges("chr1", IRanges(start=100, width=151), strand="+", ID="gene1_exon2")
gr3 <- GRanges("chr1", IRanges(start=900, width=101), strand="+", ID="gene1_exon3")
gr4 <- GRanges("chr1", IRanges(start=500, width=101), strand="+", ID="gene2_exon1")
gr <- c(gr1, gr2, gr3, gr4)
grl <- split(gr, gsub("_.*", "", elementMetadata(gr)$ID))

## Count reads overlapping with exonic regions of genes
assays(summarizeOverlaps(grl, rd, mode="Union", ignore.strand=TRUE))$counts
countOverlaps(grl, rd, ignore.strand=TRUE)

## Compare results with corresponding gene range counts
gr1 <- GRanges("chr1", IRanges(start=100, width=901), strand="+", ID="gene1")
gr2 <- GRanges("chr1", IRanges(start=500, width=101), strand="+", ID="gene2")
grgene <- c(gr1, gr2)
grgenel <- split(grgene, elementMetadata(grgene)$ID)
assays(summarizeOverlaps(grgenel, rd, mode="Union", ignore.strand=TRUE))$counts
countOverlaps(grgenel, rd, ignore.strand=TRUE)
