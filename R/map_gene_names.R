#map genes to HGNC-approved symbols

HGNCapproved = function(genes, untranslatable.na = F){
  message("Using HGNC 2018 approved symbols.")
  idx = match(genes,genesetr::hgnc_dict_2018$synonyms)
  genes[!is.na(idx)] = genesetr::hgnc_dict_2018[na.omit(idx),"approved"]
  if(untranslatable.na == T) genes[is.na(idx)] = NA
  return(genes)
}

#map all genes in a library to HGNC-approved symbols
#should change untranslatable.rm to untranslatable.na
lib2HGNC = function(lib, untranslatable.na = F){
  message("Using HGNC 2018 approved symbols.")
  message("Note gene symbols appearing in gene set names will not be mapped.")

  if(genesetr::isList(lib)){
    return(lapply(lib,function(set){
      return(genesetr::HGNCapproved(set,untranslatable.na))
    }))

  }else if(genesetr::isMatrix(lib)){
    rownames(lib) = genesetr::HGNCapproved(rownames(lib),untranslatable.na)
    return(lib)

  }else if(genesetr::isLongDF(lib)){
    lib$gene = HGNCapproved(lib$gene,untranslatable.na)
    return(lib)

  }else{
    error("Gene set library is in unrecognized format.")
  }

}

ensembl2HGNC = function(ids, untranslatable.na = F){
  #get everything to ensembl gene id
  p_match = match(ids, genesetr::ensembl_ids$ensembl_peptide_id)
  t_match = match(ids, genesetr::ensembl_ids$ensembl_transcript_id)

  ids[!is.na(t_match)] = genesetr::ensembl_ids[na.omit(t_match),"ensembl_gene_id"]
  ids[!is.na(p_match)] = genesetr::ensembl_ids[na.omit(p_match),"ensembl_gene_id"]

  #get everything to HGNC id
  hgnc_match = match(ids,genesetr::hgnc_2018$Ensembl.ID.supplied.by.Ensembl.)
  ids[!is.na(hgnc_match)] = genesetr::hgnc_2018[na.omit(hgnc_match), "Approved.Symbol"]
  if(untranslatable.na) ids[is.na(hgnc_match)] = NA
  return(ids)
}
