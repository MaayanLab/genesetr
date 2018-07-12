#' Load GMT from disk.
#'
#' @param filename The path to the file.
#' @return A named list of character vectors.
#'
#' A GMT is a tab-separated file where each line can have a variable number of terms where
#' each line corresponds to a set (e.g. a gene set). The first term in each line
#' corresponds to the set label.
#'
#' @examples
#' gmt = loadGMT("user/dir/my.gmt")

loadGMT = function(filename){
  lib = genesetr::read.ttsv(filename)
  lib = lapply(lib,function(y) return(y = as.character(y[!is.na(y)])))
  names(lib) = toupper(gsub(" ","",gsub("-","",gsub("\\.","",names(lib)))))
  genesetr::checkLib(lib)
  return(lib)
}

#' Load a MatLib from disk.
#'
#' @param filename The path to the file.
#' @return A sparse matrix representing the gene library.
#'
#' A MatLib is a (usually sparse) matrix representing a gene set library. The column labels are gene sets and the rows are genes. Membership of a gene to a gene set is designated by a non-zero entry in the matrix.
#'
#' @examples
#' gmt = loadMat("user/dir/my.gmt")
loadMat = function(filename){
return(read.table(filename, header = T, row.names = T,
  stringsAsFactors=F, quote="", comment.char="", sep="\t"))
}

#' Load a LongDF format gene set library from disk.
#'
#' @param filename The path to the file.
#' @return A two-column data frame representing a gene set library.
#'
#' The first column is the gene members and the second caolumn is the gene set names.
#'
#' @examples
#' gmt = loadLongDF("user/dir/my.gmt")
loadLongDF = function(filename){
return(read.table(filename, header = T, stringsAsFactors=F,
  quote="", comment.char="", sep="\t"))
}

#' Check library format compliance with genesetr formats.
#'
#' @param filename Library object.
#' @return TRUE if compliant, FALSE with accompanying warning if not compliant.
#'
#' The first column is the gene members and the second caolumn is the gene set names.
#'
#' @examples
#' gmt = loadLongDF(myLib)
#' print(checkLib(myLib))
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

#' Remove duplicate genes from gene set.
#'
#' @param filename Genesetr library object.
#' @return Returns library object with duplicate genes removed.
#'
#'
#'
#' @examples
#' gmt = loadLongDF('/user/libs/mylib.tsv')
#' myLib = removeDupes(myLib)
removeDupes = function(x){
  genesetr::checkLib(x)
  if(genesetr::isLongDF(x)) return(x[!duplicated(paste(x[,1],x[,2]))])
  else if(genesetr::isList(x)) return(lapply(x, function(y)return(y[!duplicated(y)])))
  else if(genesetr::isMat(x)) return(x)
  else error("Unrecognized library format.")

}

#' Remove empty gene sets from a library
#'
#' @param filename Genesetr library object.
#' @return Returns library object with empty sets removed.
#'
#' @examples
#' gmt = loadLongDF('/user/libs/mylib.tsv')
#' myLib = removeEmptySets(myLib)
removeEmptySets = function(x){
  genesetr::checkLib(x)
  if(genesetr::isLongDF(x)) return(x)
  else if(genesetr::isList(x)) return(x[lapply(x,length)>0])
  else if(genesetr::isMat(x)) return(x[,colSums(x)>0])
  else error("Unrecognized library format.")
}

dupesWarn = function(){
  # warning("Set(s) in your library contain duplicated genes. Suggested: removeDupes().")
}

emptyWarn = function(){
  # warning("Your library contains empty set(s). Suggested: removeEmptySets().")
}
dfFormatWarn = function(){
  warning("Long DF format: first column should contain",
    "/'gene/' and second column should contain /'set_name/'.")
}

