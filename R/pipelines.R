#' All Pairwise Comparisons Between Gene Sets of Two Libraries
#'
#' @param lib1 A library of gene sets in List, LongDF or Mat format.
#' @param lib2 A library of gene sets in List, LongDF or Mat format.
#' @param background The size of the "universe" of genes to use in building contingency tables from pairwise intersections of gene sets in the two libraries. Default is NULL-- background will be computed as the size of the union of unique genes in lib1 and lib2.
#' @param method Metric by which to assess gene set association. Options are "FET" and "OR".
#' @return A data frame with all pairs of genesets and association metrics
#' @examples
#' #load example gene set libraries
#' data("encode_example", "chea_example")
#' intersects = pairwiseSetOverlap(encode_example, chea_example)

pairwiseSetOverlap = function(lib1,lib2,background = NULL, method = "FET"){
  tab = genesetr::pairwiseContingTables(lib1,lib2, background)
  if("FET" %in% method){
    tab$FET.p.val = genesetr::fastFET(tab$a,tab$b,tab$c, tab$d)
  }
  if("OR" %in% method){
    tab$OR = genesetr::oddsRatio(tab$a,tab$b,tab$c,tab$d)
  }
  if("jaccard" %in% method){
    tab$jaccard = genesetr::jaccardInd(tab$a,tab$b,tab$c,tab$d)
  }
  return(tab)
}


