write.genepop <- function(genind_obj, outfilebase, genepop_header) {
  
  df1 <- genind2df(genind_obj, oneColPerAll = TRUE)
  #figure out which columns only have 2 digits (length of 2)
  genos <- df1 %>%
    dplyr::select(-pop)
  temp <- genos %>% 
    dplyr::mutate_all(as.numeric)
  
  twoers <- temp[temp < 100 & !(is.na(temp))]
  twoers_conv <- paste0("0",as.character(twoers))
  
  #replace
  genos[temp < 100 & !(is.na(temp))] <- twoers_conv
  
  loc1 <- seq(1, ncol(genos), 2)
  loc2 <- loc1+1
  
  newgenos <- as.data.frame(matrix(nrow=nrow(genos), ncol=nLoc(genind_obj)))
  
  for(i in 1:length(loc1)) {
   newgenos[,i] <- paste0(genos[,loc1[i]], genos[,loc2[i]])
   tempname <- names(genos)[loc1[i]]
   names(newgenos[,i]) <- gsub(pattern = "\\.1", "", tempname)
  }
  
  # replace NANA with NA
  newgenos[newgenos=="NANA"] <- NA
  
  #put sampleID row names and pop column back on
  df_dat <- cbind(pop=df1$pop, newgenos)
  rownames(df_dat) <- row.names(df1)
  
  allele.digits <- nchar(genind_obj@all.names[[1]][1])
  df_dat[is.na(df_dat)] <- paste(rep("0", allele.digits * 2), 
                                 collapse = "")
  numloci <- length(names(genind_obj@all.names))
  N <- dim(genind_obj@tab)[1]
  OS <- as.character(Sys.info()["sysname"])
  genout_file <- paste(outfilebase, "_genepop.gen", sep = "")
  write.table(genepop_header, genout_file, sep = "\t", 
              quote = F, col.names = F, row.names = F)
  for (locus in 1:numloci) {
    write.table(names(genind_obj@all.names)[locus], genout_file, 
                sep = "\t", quote = F, col.names = F, row.names = F, 
                append = T)
  }
  for (pop in unique(pop(genind_obj))) {
    write.table("POP", genout_file, sep = "\t", quote = F, 
                col.names = F, row.names = F, append = T)
    output <- df_dat[which(pop(genind_obj) == pop),2:ncol(df_dat)]
    output <- data.frame( paste0(rownames(output),",") ,output)
    write.table(output, genout_file, sep = "\t", quote = F, 
                col.names = F, row.names = F, append = T)
    #cbind(row.names(output),output)
  } 
}
