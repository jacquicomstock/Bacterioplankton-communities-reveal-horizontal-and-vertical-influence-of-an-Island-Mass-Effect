
path <- "/home/ec2-user"
fnFs <- sort(list.files(path, pattern="_R1_001.fastq", full.names = TRUE))
fnRs <- sort(list.files(path, pattern="_R2_001.fastq", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_L"), `[`, 1)
filt_path <- file.path(path, "filtered")
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq"))
filtRs <- file.path(filt_path, paste0(sample.names, "_R_filt.fastq"))
for(i in seq_along(fnFs)) {
  fastqPairedFilter(c(fnFs[i], fnRs[i]), c(filtFs[i], filtRs[i]),
                    truncLen=c(220,130), 
                    maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                    compress=TRUE, verbose=TRUE)
}
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
names(derepFs) <- sample.names
names(derepRs) <- sample.names
dadaFs.lrn <- dada(derepFs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errF <- dadaFs.lrn[[1]]$err_out
dadaRs.lrn <- dada(derepRs, err=NULL, selfConsist = TRUE, multithread=TRUE)
errR <- dadaRs.lrn[[1]]$err_out
saveRDS(errF, "/home/ec2-user/errF.rds")
saveRDS(errR, "/home/ec2-user/errR.rds")
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)
saveRDS(dadaFs, "/home/ec2-user/dadaFs_N.rds")
saveRDS(dadaRs, "/home/ec2-user/dadaRs_N.rds")
mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
saveRDS(mergers, "/home/ec2-user/dadaFs_N.rds")
seqtab <- makeSequenceTable(mergers[names(mergers) != "Mock"])
dim(seqtab)
table(nchar(getSequences(seqtab)))
seqtab.nochim <- removeBimeraDenovo(seqtab, verbose=TRUE)
sum(seqtab.nochim)/sum(seqtab)
dim(seqtab.nochim)
saveRDS(seqtab.nochim, "/home/ec2-user/dadaFs_Nseqtab.nochim.rds")
taxa <- assignTaxonomy(seqtab.nochim, "/home/ec2-user/silva_nr_v132_train_set.fa")
unname(head(taxa))
saveRDS(taxa, "/home/ec2-user/taxa.rds")
write.table(cbind(t(seqtab.nochim) , taxa), "/home/ec2-user/seqtab-nochimtaxa.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
write.table(taxa,"/home/ec2-user/taxa.txt", sep="\t", row.names=TRUE, col.names=NA, quote=FALSE)
