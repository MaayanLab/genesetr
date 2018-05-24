loadGMT = function(filename){
  lib = genesetr::read.ttsv(filename)
  lib = lapply(lib,function(y) return(y = as.character(y[!is.na(y)])))
  names(lib) = toupper(gsub(" ","",gsub("-","",gsub("\\.","",names(lib)))))
  genesetr::checkLib(lib)
  return(lib)
}

loadMat = function(filename){

}

loadLongDF = function(filename){

}

checkLib = function(x){
  if(genesetr::isLongDF(x)){
    if(any(colnames(x) != c("gene","set_name")) || ncol(x)!=2) genesetr::dfFormatWarn()
    if(any(duplicated(paste(x[,1],x[,2])))) genesetr::dupesWarn()

  }else if(genesetr::isList(x)){
    if(any(unlist(lapply(x,length))<1)) genesetr::emptyWarn()
    if(any(unlist(lapply(x,function(y){return(any(duplicated(y)))})))) genesetr::dupesWarn()
  }else if(genesetr::isMat(x)){
    if(any(colSums(x)==0)) genesetr::emptyWarn()
    if(any(!(mat == 1 || mat == 0))) warning("Matrix library format values should be [0,1]")
    if(any(mat > 1)) genesetr::dupesWarn()
    if(nrow(x)<ncol(x)) warning(paste("Num cols detected to be > than num rows. ",
      "Note that for matrix library format, genes are rows and sets are columns."))
  }else{
    error("Library is in unrecognized format.")
  }
}

removeDupes = function(x){
  genesetr::checkLib(x)
  if(genesetr::isLongDF(x)) return(x[!duplicated(paste(x[,1],x[,2]))])
  else if(genesetr::isList(x)) return(lapply(x, function(y)return(y[!duplicated(y)])))
  else if(genesetr::isMat(x)) return(x)
  else error("Unrecognized library format.")

}
removeEmptySets = function(x){
  genesetr::checkLib(x)
  if(genesetr::isLongDF(x)) return(x)
  else if(genesetr::isList(x)) return(x[lapply(x,length)>0])
  else if(genesetr::isMat(x)) return(x[,colSums(x)>0])
  else error("Unrecognized library format.")
}

dupesWarn = function(){
  warning("Set(s) in your library contain duplicated genes. Suggested: removeDupes().")
}

emptyWarn = function(){
  warning("Your library contains empty set(s). Suggested: removeEmptySets().")
}
dfFormatWarn = function(){
  warning("Long DF format: first column should contain",
    "/'gene/' and second column should contain /'set_name/'.")
}

