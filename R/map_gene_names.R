#' Determine whether there are weights associated with genes.
#'
#' @param genes a character vector of genes
#' @param untranslatable.na boolean, default is false. If a gene does not have an HGNC
#' mapping, then the function will return an NA for that gene
#' @param weighted boolean representing whether genes have weights
#' @return a character vector with genes mapped to approved HGNC symbols
#'
#' @examples
#' geneset = c("FOXO1", "MYC", "GATA3", "BP1")
#' mapped_geneset = HGNCapproved(geneset, untranslatable.na = T)
#'


HGNCapproved = function(genes, untranslatable.na = F, weighted = F){

  if(weighted){
    scores = genesetr::getStringElement(genes,",",2)
    genes = genesetr::getStringElement(genes,",",1)
  }
  idx = match(genes,names(genesetr::hgnc_dict))
  genes[!is.na(idx)] = genesetr::hgnc_dict[na.omit(idx)]
  if(untranslatable.na == T) genes[is.na(idx)] = NA
  if(weighted){
    genes = paste(genes,scores,sep = ",")
  }
  return(genes)
}

#' Map all genes in a genesetr library object to HGNC-approved symbols
#'
#' @param lib a genesetr library object
#' @param untranslatable.na boolean, default is false. If a gene does not have an HGNC
#' mapping, then the function will return an NA for that gene
#'
#' @return a genesetr library object with all genes mapped to approved HGNC symbols
#'
#' @examples
#' lib = genesetr::chea_example
#' mapped_lib = lib2HGNC(lib, untranslatable.na = T)
#'

lib2HGNC = function(lib, untranslatable.na = F){
  message("Using HGNC 2018 approved symbols.")
  message("Note gene symbols appearing in gene set names will not be mapped.")

  if(genesetr::isList(lib)){
    weighted = isWeightedLib(lib)
    return(lapply(lib,function(set){
      return(genesetr::HGNCapproved(set,
        untranslatable.na,
        weighted = weighted))
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

#' Maps Ensembl transcript and gene ids to HGNC-approved gene symbols
#'
#' @param ids a character vector of ensembl ids
#' @param untranslatable.na boolean, default is false. If an id does not have an HGNC
#' mapping, then the function will return an NA for that gene
#'
#' @return a character vector with all ids mapped to approved HGNC symbols
#'
#' @examples
#' ids = genesetr::ensembl_ids$ensembl_gene_id[1:10]
#' HGNC_symbols = ensembl2HGNC(ids, untranslatable.na = T)
#'

ensembl2HGNC = function(ids, untranslatable.na = F){
  #get everything to ensembl gene id
  p_match = match(ids, genesetr::ensembl_ids$ensembl_peptide_id)
  t_match = match(ids, genesetr::ensembl_ids$ensembl_transcript_id)

  ids[!is.na(t_match)] = genesetr::ensembl_ids[na.omit(t_match),"ensembl_gene_id"]
  ids[!is.na(p_match)] = genesetr::ensembl_ids[na.omit(p_match),"ensembl_gene_id"]

  #get everything to HGNC id
  hgnc_match = match(ids,genesetr::hugo_df$Ensembl.ID.supplied.by.Ensembl.)
  ids[!is.na(hgnc_match)] = genesetr::hugo_df[na.omit(hgnc_match), "Approved.Symbol"]
  if(untranslatable.na) ids[is.na(hgnc_match)] = NA
  return(ids)
}

