## ═══════════════════════════════════════════════════════════════════════════ ##
##  1. Input data
## ═══════════════════════════════════════════════════════════════════════════ ##
file.to.Jurkat.ENCODE.pool <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIVasRNAMAP", "Illumina.seq", "RRBS.seq", "RRBS3.pHCC5", "analyse", "RRBS.Jurkat.mapped.sorted.pools.sam")

file.path.to.814.pos <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "analyse", "RRBS1_bismark_bt2.mapped.sorted.sam") # pooled replicates
file.path.to.814.pos.rep1 <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "analyse", "RRBS1.rep1.mapped.sorted.sam")
file.path.to.814.pos.rep2 <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "analyse", "RRBS1.rep2.mapped.sorted.sam")

file.path.to.814.neg <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "analyse", "RRBS2_bismark_bt2.mapped.sorted.sam") # pooled replicates
file.path.to.814.neg.rep1 <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "analyse", "RRBS2.rep1.mapped.sorted.sam")
file.path.to.814.neg.rep2 <- file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "analyse", "RRBS2.rep2.mapped.sorted.sam")

###
f814.pos <- processBismarkAln(location = file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "RRBS1", "RRBS1_bismark_bt2.mapped.sorted.sam"), sample.id = "f814.pos", assembly = "hg38")
f814.neg <- processBismarkAln(location = file.path("","home", "labadmin", "Documents", "Arbeitplatz", "Projekten", "HIV.epigenetic.inheritance", "Illumina.seq", "RRBS.seq", "RRBS2", "RRBS2_bismark_bt2.mapped.sorted.sam"), sample.id = "f814.neg", assembly = "hg38")

## ═══════════════════════════════════════════════════════════════════════════ ##
##  2. File list
## ═══════════════════════════════════════════════════════════════════════════ ##
file.path.pool <- list(file.path.to.814.pos, file.path.to.814.neg)
f814 <- processBismarkAln(location = file.path.pool, sample.id = list("f814.pos", "f814.neg"), treatment = c(0, 1), assembly = "hg38")

file.path.814.reps <- list(file.path.to.814.pos.rep1, file.path.to.814.pos.rep2, file.path.to.814.neg.rep1, file.path.to.814.neg.rep2)
f814.reps <- processBismarkAln(location = file.path.814.reps, sample.id = list("f814.pos.rep1", "f814.pos.rep2", "f814.neg.rep1", "f814.neg.rep2"), treatment = c(0, 0, 1, 1), assembly = "hg38")

file.path.Jurkat.814.pos <- list(file.to.Jurkat.ENCODE.pool, file.path.to.814.pos)
f814.pos.Jurkat.norl <- processBismarkAln(location = file.path.Jurkat.814.pos, sample.id = list("Jurkat", "f814.pos"), treatment = c(0, 1), assembly = "hg38")

file.path.Jurkat.814.neg <- list(file.to.Jurkat.ENCODE.pool, file.path.to.814.neg)
f814.neg.Jurkat.norl <- processBismarkAln(location = file.path.Jurkat.814.neg, sample.id = list("Jurkat", "f814.neg"), treatment = c(0, 1), assembly = "hg38")

## ═══════════════════════════════════════════════════════════════════════════ ##
##  3. Filtering samples based on read coverage
## ═══════════════════════════════════════════════════════════════════════════ ##
f814.filter <- filterByCoverage(f814, lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
f814.reps.filter <- filterByCoverage(f814.reps, lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
f814.pos.Jurkat.norl.filter <- filterByCoverage(f814.pos.Jurkat.norl, lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
f814.neg.Jurkat.norl.filter <- filterByCoverage(f814.neg.Jurkat.norl, lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9) 

###
f814.pos.filter <- filterByCoverage(f814.pos, lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
f814.neg.filter <- filterByCoverage(f814.neg, lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)

## ═══════════════════════════════════════════════════════════════════════════ ##
##  4. Merging samples
## ═══════════════════════════════════════════════════════════════════════════ ##
f814.filter.merge <- methylKit::unite(f814.filter, destrand = F)
f814.reps.filter.merge <- methylKit::unite(f814.reps.filter, destrand = F)
f814.pos.Jurkat.norl.filter.merge <- methylKit::unite(f814.pos.Jurkat.norl.filter, destrand = F)
f814.neg.Jurkat.norl.filter.merge <- methylKit::unite(f814.neg.Jurkat.norl.filter, destrand = F) 

## ═══════════════════════════════════════════════════════════════════════════ ##
##  5. PCA
## ═══════════════════════════════════════════════════════════════════════════ ##
PCASamples(f814.reps.filter.merge)

## ═══════════════════════════════════════════════════════════════════════════ ##
##  6. Differentially methylated bases (DMBs)
## ═══════════════════════════════════════════════════════════════════════════ ##
f814.filter.merge.diff <- calculateDiffMeth(f814.filter.merge, overdispersion = "MN", adjust = "BH")
f814.pos.Jurkat.norl.filter.merge.diff <- calculateDiffMeth(f814.pos.Jurkat.norl.filter.merge, overdispersion = "MN", adjust = "BH")
f814.neg.Jurkat.norl.filter.merge.diff <- calculateDiffMeth(f814.neg.Jurkat.norl.filter.merge, overdispersion = "MN", adjust = "BH")

## ═══════════════════════════════════════════════════════════════════════════ ##
##  7. Differentially methylated regions (DMRs)
## ═══════════════════════════════════════════════════════════════════════════ ##
f814.filter.merge.tiles <- tileMethylCounts(f814.filter.merge, win.size=1000,step.size=1000)
f814.filter.merge.tiles.diff <- calculateDiffMeth(f814.filter.merge.tiles, mc.cores = 2)

f814.pos.Jurkat.norl.filter.merge.tiles <- tileMethylCounts(f814.pos.Jurkat.norl.filter.merge, win.size=1000,step.size=1000)
f814.pos.Jurkat.norl.filter.merge.tiles.diff <- calculateDiffMeth(f814.pos.Jurkat.norl.filter.merge.tiles, mc.cores = 2)

f814.neg.Jurkat.norl.filter.merge.tiles <- tileMethylCounts(f814.neg.Jurkat.norl.filter.merge, win.size=1000,step.size=1000)
f814.neg.Jurkat.norl.filter.merge.tiles.diff <- calculateDiffMeth(f814.neg.Jurkat.norl.filter.merge.tiles, mc.cores = 2)

## ═══════════════════════════════════════════════════════════════════════════ ##
##  8. Volcano plots
## ═══════════════════════════════════════════════════════════════════════════ ##

## Functions
S4.to.df.diff.sig <- function(S4) {
  df <- data.frame(chr = S4$chr, start = S4$start, end = S4$end, pvalue = S4$pvalue, qvalue = S4$qvalue, meth.diff = S4$meth.diff)
  
  df.sig <- df %>% dplyr::mutate(sig = case_when(pvalue < 0.05 & qvalue < 0.05 & meth.diff > 25 ~ "up", pvalue < 0.05 & qvalue < 0.05 & meth.diff < -25 ~ "down", pvalue > 0.05 | qvalue > 0.05 ~ "no.sig")) %>% na.omit()
  
  df.sig <- df.sig %>% dplyr::filter(chr == "1" | chr == "2" | chr == "3" | chr == "4" | chr == "5" | chr == "6" | chr == "7" | chr == "8" | chr == "9" | chr == "10" | chr == "11" | chr == "12" | chr == "13" | chr == "14" | chr == "15" | chr == "16" | chr == "17" | chr == "18" | chr == "19" | chr == "20" | chr == "21" | chr == "22" | chr == "X")

Richtung.chr.num <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")

df.sig$chr <- factor(df.sig$chr, levels = Richtung.chr.num)
  
  return(df.sig)
}

volcano.plotting <- function(df) {
  
  df.sig <- df %>% dplyr::filter(sig != "no.sig")
 
  df.no.sig.frac <- df %>% dplyr::filter(sig == "no.sig") %>% dplyr::sample_frac(0.3)
  
  df.chimeric <- bind_rows(df.sig, df.no.sig.frac)
  
  Kolory_volcano <- c("orangered","skyblue2","grey")
  names(Kolory_volcano) <- c("up", "down", "no.sig")

 df.chimeric$label <- NA
 df.chimeric$label[df.chimeric$sig != "no.sig"] <- df.chimeric$gene[df.chimeric$sig != "no.sig"]

 top.list <- df.chimeric[order(df.chimeric$pvalue), ] 
 top.list$top <- ""
 top.list$top[1:20] <- top.list$gene[1:20]
  
  ggplot(df.chimeric, aes(x = meth.diff, y = -log10(pvalue)))+geom_point(aes(col = sig), alpha = 0.5)+theme_bw()+scale_color_manual(values= Kolory_volcano)+geom_text_repel(data = top.list, aes(label = top), max.overlaps=Inf)+theme(axis.title.x=element_text(size=12), axis.text.x=element_text(size=12, colour = "black",angle=0,vjust=0.5), axis.title.y = element_text(size=12),axis.text.y = element_text(size = 12, colour = "black"))+xlab("Differential methylation (genomic region)")+ylab("P value (-log10)")
}












