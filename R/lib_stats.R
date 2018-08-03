getLibInfo = function(lib){
  lib = genesetr::toList(lib)
  lib = genesetr::removeLibWeights(lib)
  set_lengths = unlist(lapply(lib,length))
  quart = as.numeric(quantile(lib_lengths))
  return(data.frame(
    n_gene_sets = length(lib),
    n_unique_genes = length(unique(unlist(lib))),
    mean_lib_length = mean(set_lengths),
    sd_set_length = sd(set_lengths),
    quartile1_set_length = quart[1],
    quartile2_set_length = quart[2],
    quartile3_set_length = quart[3],
    quartile4_set_length = quart[4]
  ))
}

getLibGeneFreqs = function(lib){
  x = as.data.frame(table(unlist(lib)))
  colnames(x)[1] = "gene"
  x$gene = as.character(x$gene)
  return(x)
}
