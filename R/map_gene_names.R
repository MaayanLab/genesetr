#map genes to HGNC-approved symbols
HGNCapproved = function(genes, untranslatable.rm = F){
  message("Using HGNC 2018 approved symbols.")
  idx = match(genes,genesetr::hgnc_dict$syn)
  genes[!is.na(idx)] = genesetr::hgnc_dict[na.omit(idx),"approved"]
  if(untranslatable.rm == T) genes[is.na(idx)] = NA
  return(genes)
}

#map all genes in a library to HGNC-approved symbols
#should change untranslatable.rm to untranslatable.na
lib2HGNC = function(lib, untranslatable.rm = F){
  message("Using HGNC 2018 approved symbols.")
  message("Note gene set names will not be mapped.")

  if(genesetr::isList(lib)){
    return(lapply(lib,function(set){
      return(genesetr::HGNCapproved(set,untranslatable.rm))
    }))

  }else if(genesetr::isMatrix(lib)){
    rownames(lib) = genesetr::HGNCapproved(rownames(lib),untranslatable.rm)
    return(lib)

  }else if(genesetr::isLongDF(lib)){
    lib$gene = HGNCapproved(lib$gene,untranslatable.rm)
    return(lib)

  }else{
    error("Gene set library is in unrecognized format.")
  }

}
