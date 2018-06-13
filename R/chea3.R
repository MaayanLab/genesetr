queryChea3 = function(
  geneset,
  set_name = "untitled",
  n_results = 10,
  libnames = names(genesetr::libs),
  integrate_libs = T,
  method = c("FET"),
  background = NULL){

  query = list(geneset)
  names(query) = set_name


  if(!all(libnames %in% names(genesetr::libs))){
    error(paste("The following libraries are unavailable:",
      paste(setdiff(libnames,names(genesetr::libs)),collapse = " "),
      "These are the libraries available:",
      paste(names(genesetr::libs),collapse = " ")))}

  temp_libs = genesetr::libs[libnames]
  results = lapply(temp_libs,function(x){
    df = genesetr::pairwiseSetOverlap(query, x,
      background = background, method = method)
    df = df[order(df$FET.p.val),][1:n_results,]
    return(df)
  })
  return(jsonlite::toJSON(results))
}



