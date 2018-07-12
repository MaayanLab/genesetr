getLibInfo = function(lib){
  lib = genesetr::toList(lib)
  lib_lengths = unlist(lapply(lib,length))
  quart = as.numeric(quantile(lib_lengths))
  return(data.frame(
    n_gene_sets = length(lib),
    n_unique_genes = length(unique(unlist(lib))),
    mean_lib_length = mean(lib_lengths),
    sd_lib_length = sd(lib_lengths),
    quartile1_lib_length = quart[1],
    quartile2_lib_length = quart[2],
    quartile3_lib_length = quart[3],
    quartile4_lib_length = quart[4]
  ))
}

getLibGeneFreqs = function(lib){
  x = as.data.frame(table(unlist(lib)))
  x$Var1 = as.character(x$Var1)
  return(x)
}
