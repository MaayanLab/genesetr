
libCompare = function(lib1,lib2, method = "jaccard"){
  lib1 = genesetr::removeLibWeights(lib1)
  lib2 = genesetr::removeLibWeights(lib2)
  results = genesetr::pairwiseContingTables(lib1, lib2, background = 20000)
  if("jaccard" %in% method){
    results$jaccard = genesetr::jaccardInd(results$a, results$b, results$c, results$d)
    return(results)
  }else{
    error("Indicate method argument.")
  }
}

