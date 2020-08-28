
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
}


library(dplyr)
library(data.table)
library(DEXSeq)

# Take a fcount file and convert it to dcounts for dexseq
message("Reading and adding Exon IDs for DEXSeq")

#--------------Cohort2--------#
#dcounts <- as.data.frame(fread(args[1]))
#dcounts = dcounts[,c(1, 7:27)]
##sampleData <- data.frame(row.names = c("CT1","CT2","AS1","AS2","AS3","AS4","AS5",
##                                       "AS6","AS7","AS8","AS9",
##                                       "CT3","CT4","CT5","CT6","CT7","CT8","CT9","CT10","CT11","CT12"), 
##                         condition = c("ctrl","ctrl","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4",
##                                       "ASCD4","ASCD4","ASCD4","ASCD4",
##                                       "ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl"))

# #--------------Cohort1+2--------#
# dcounts <- as.data.frame(fread(args[1]))
# dcounts = dcounts[,c(1,7:57)]
# sampleData <- data.frame(row.names = c("CT1","CT2",
#                                        "AS1","AS2","AS3","AS4","AS5","AS6","AS7","AS8","AS9","AS10",
#                                        "AS11","AS12","AS13","AS14","AS15","AS16","AS17","AS18",
#                                        "CT3","CT4","CT5","CT6","CT7","CT8","CT9","CT10","CT11","CT12",
#                                        "CT13","CT14","CT15","CT16","CT17","CT18","CT19","CT20","CT21","CT22",
#                                        "CT23","CT24","CT25","CT26","CT27","CT28","CT29","CT30","CT31","CT32",
#                                        "CT33"), 
#                          condition = c("ctrl","ctrl",
#                                        "ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4",
#                                        "ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4",
#                                        "ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
#                                        "ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
#                                        "ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
#                                        "ctrl"))
# colnames(dcounts) <- c("GeneID", rownames(sampleData) )

#--------------Cohort1   --------#
dcounts <- as.data.frame(fread(args[1]))
dcounts = dcounts[,c(1,9:17,27:43,54:57)]
sampleData <- data.frame(row.names = c("AS1","AS2","AS3","AS4","AS5","AS6","AS7","AS8","AS9",
                                       "CT1","CT2","CT3","CT4","CT5","CT6","CT7","CT8","CT9","CT10","CT11","CT12",
                                       "CT13","CT14","CT15","CT16","CT17","CT18","CT19","CT20","CT21"), 
                         condition = c("ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4","ASCD4",
                                       "ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
                                       "ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl","ctrl",
                                       "ctrl"))

colnames(dcounts) <- c("GeneID", rownames(sampleData) )


# # remove the non-expressed genes
# x= transform(dcounts, amount1=ave(UT1, GeneID, FUN=sum),amount2=ave(HD1, GeneID, FUN=sum),
#                     amount3=ave(as.numeric(as.character(UT2)), GeneID, FUN=sum),amount4=ave(HD2, GeneID, FUN=sum),
#                     length=ave(as.character(UT1), GeneID, FUN=length))
# x2 <- x[ rowSums(x[,6:9]) >=10, ]
# dcounts = dcounts[dcounts$GeneID %in% x2$GeneID,]

x = dcounts %>%
  group_by(GeneID) %>%
  summarise_all(funs(sum))
x2= transform(x, gene.exp=rowSums(x[,-1], na.rm=T))
x3 <- x2[x2$gene.exp > 100, ]
dcounts = dcounts[dcounts$GeneID %in% x3$GeneID,]


#####
id <- as.character(dcounts[,1])
n <- id
split(n,id) <- lapply(split(n ,id), seq_along )
rownames(dcounts) <- sprintf("%s%s%03.f",id,":E",as.numeric(n))
dcounts <- dcounts[,2:ncol(dcounts)]

dcounts <- dcounts[substr(rownames(dcounts), 1, 1) != "_", ] #remove _ from beginnning of gene name 

## get genes and exon names out
splitted <- strsplit(rownames(dcounts), ":")
exons <- sapply(splitted, "[[", 2)
genesrle <- sapply(splitted, "[[", 1)


## parse the flattened file
aggregates <- read.delim("/well/jknight/ping/gtfs/Homo_sapiens.GRCh37.87_DEXSeq.gtf", stringsAsFactors = FALSE, 
                         header = FALSE)
colnames(aggregates) <- c("chr", "source", "class", "start", 
                          "end", "ex", "strand", "ex2", "attr")
aggregates$strand <- gsub("\\.", "*", aggregates$strand)
aggregates <- aggregates[which(aggregates$class == "exonic_part"),] #########
aggregates$attr <- gsub("\"|=|;", "", aggregates$attr)
aggregates$gene_id <- sub(".*gene_id\\s(\\S+).*", "\\1", 
                          aggregates$attr)
# trim the gene_ids to 255 chars in order to match with featurecounts
longIDs <- sum(nchar(unique(aggregates$gene_id)) > 255)
#warning(paste0(longIDs, 
#             " aggregate geneIDs were found truncated in featureCounts output"), 
#       call. = FALSE)
aggregates$gene_id <- substr(aggregates$gene_id,1,255)

transcripts <- gsub(".*transcripts\\s(\\S+).*", "\\1", 
                    aggregates$attr)
transcripts <- strsplit(transcripts, "\\+")
exonids <- gsub(".*exonic_part_number\\s(\\S+).*", "\\1", # exonic_part_number
                aggregates$attr)
exoninfo <- GRanges(as.character(aggregates$chr), IRanges(start = aggregates$start, 
                                                          end = aggregates$end), strand = aggregates$strand)
names(exoninfo) <- paste(aggregates$gene_id, exonids, 
                         sep = ":E")

names(transcripts) <- names(exoninfo) 
#  if (!all(rownames(dcounts) %in% names(exoninfo))) {
#     stop("Count files do not correspond to the flattened annotation file")
#  }
matching <- match(rownames(dcounts), names(exoninfo))
stopifnot(all(names(exoninfo[matching]) == rownames(dcounts)))
stopifnot(all(names(transcripts[matching]) == rownames(dcounts)))
dxd <- DEXSeqDataSet(dcounts, sampleData, design = ~sample + exon + condition:exon, exons, 
                     genesrle, exoninfo[matching], transcripts[matching])


# ######  the wrapper function DEXseq(), which runs the analysis shown above through a single function call.
# dxr = DEXSeq(dxd)
#---------------------------------------------
# BPPARAM <- MulticoreParam(workers = 8)
#---------------------------------------------

dxd = estimateSizeFactors( dxd )
message("##Estimating dispersions...")
dxd = estimateDispersions( dxd )
message("##Running test...")
dxd = testForDEU( dxd )
message("##estimate relative exon usage fold changes...")
dxd = estimateExonFoldChanges( dxd,  fitExpToVar="condition")
message("##DEXSeqResults...")
res <- DEXSeqResults(dxd)
#---------------------------------------------
message("##significant DEU cases...")
table( res$padj < 0.05 )
#message("generate HTML report...")
#DEXSeqHTML( res, FDR=0.05, color=c("#FF000080", "#0000FF80") )
#---------------------------------------------#
message("##Save results")
results = as.data.frame(res)
write.csv(results, "DEXSeq.results.CD4_cohort1.csv")

write.table(res, 
            file = paste0("DEXSeq.results.CD4_cohort1", ".txt"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")
message("Summarizing results on gene level...")
pgq <- perGeneQValue(res, p = "pvalue")
save(dxd, res, pgq, file = paste0("DEXSeq.results.CD4_cohort1", ".Rdata"))
#tmp <- cbind(gene = names(pgq), "adjP" = pgq)
#write.table(tmp, 
#       file = paste0("AS.CD4.DEXSeq_gene.level", ".txt"), 
#    col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")