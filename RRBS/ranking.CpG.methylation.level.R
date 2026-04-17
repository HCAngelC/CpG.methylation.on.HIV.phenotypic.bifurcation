S4.to.df <- function(S4) {
  df <- data.frame(chr = S4$chr, start = S4$start, end = S4$end, coverage = S4$coverage, numCs = S4$numCs, numTs = S4$numTs)
  
  df <- df %>% dplyr::filter(chr == "1" | chr == "2" | chr == "3" | chr == "4" | chr == "5" | chr == "6" | chr == "7" | chr == "8" | chr == "9" | chr == "10" | chr == "11" | chr == "12" | chr == "13" | chr == "14" | chr == "15" | chr == "16" | chr == "17" | chr == "18" | chr == "19" | chr == "20" | chr == "21" | chr == "22" | chr == "X")

Richtung.chr.num <- c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12", "13", "14", "15", "16", "17", "18", "19", "20", "21", "22", "X")

df$chr <- factor(df$chr, levels = Richtung.chr.num)
  
  return(df)
}

selection.top.coverage <- function(S4) {
  
  df <- S4.to.df(S4)  

  n.top1 <- 0.01*nrow(df)
  n.top5 <- 0.05*nrow(df)
  n.top10 <- 0.1*nrow(df)
  
  df <- df %>% dplyr::arrange(desc(coverage))
  
  df.high.top1 <- df[1:n.top1,] %>% dplyr::mutate(cat = "high", rank = "top1")
  df.high.top5 <- df[1:n.top5,] %>% dplyr::mutate(cat = "high", rank = "top5")
  df.high.top10 <- df[1:n.top10,] %>% dplyr::mutate(cat = "high", rank = "top10")
  
  df.low.top1 <- df[-c(1:n.top1),] %>% dplyr::mutate(cat = "low", rank = "top1")
  df.low.top5 <- df[-c(1:n.top5),] %>% dplyr::mutate(cat = "low", rank = "top5")
  df.low.top10 <- df[-c(1:n.top10),] %>% dplyr::mutate(cat = "low", rank = "top10")
  
  df.pool <- dplyr::bind_rows(df.high.top1, df.high.top5, df.high.top10, df.low.top1, df.low.top5, df.low.top10)
  
  df.pool <- df.pool %>% dplyr::mutate(meth.percent = (numCs/(numCs+numTs)))
  
  Vergleichung.rank <- c("top1", "top5", "top10")
  
  df.pool$rank <- factor(df.pool$rank, levels = Vergleichung.rank)
  
  return(df.pool)
}
