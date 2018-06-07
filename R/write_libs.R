writeGMT = function(lib,outfile){
  if(genesetr::isMat(lib)||genesetr::isLongDF(lib)){
    lib = genesetr::toList(lib)
  }
  if(isList(lib)){
    if (file.exists(outfile)) {
      warning("Replacing existing GMT file.")
      file.remove(outfile)
    }
    #put name of set as first term in set
    new_lib = as.list(names(lib))
    names(new_lib) = names(lib)
    new_lib = lapply(new_lib,function(sub){
      return(c(sub,lib[[sub]]))
    })
w
    #write lib to disk
    lapply(X = new_lib, FUN = function(x) {
      write(x, append = T, file = outfile, ncolumns = length(x), sep = "\t")
    })
  }else{
    error("Library not in correct format to write gmt.")
  }
}
