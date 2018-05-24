pairwiseSetOverlap = function(lib1,lib2,background = NULL, method = c("FET","OR")){
  tab = genesetr::pairwiseContingTables(lib1,lib2, background)
  if("FET" %in% method){
    tab$FET.p.val = genesetr::fastFET(tab$a,tab$b,tab$c, tab$d)
  }
  if("OR" %in% method){
    tab$OR = genesetr::oddsRatio(tab$a,tab$b,tab$c,tab$d)
  }
  return(tab)
}

shuffleGeneSetNames = function(lib){

}
