
libCompare = function(lib1,lib2, method = c("jaccard index")){
  results = genesetr::pairwiseContingTables(lib1, lib2, background = 20000)
  if(method == "jaccard index"){
    results$jaccard = genesetr::jaccardInd(results$a, results$b, results$c, results$d)
    return(results)
  }else{
    error("Indicate method argument.")
  }
}

