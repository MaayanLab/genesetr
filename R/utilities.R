read.ttsv = function(file, header=TRUE, sep="\t", ...) {

  n = max(count.fields(file, sep=sep), na.rm=TRUE)
  x = readLines(file)

  .splitvar = function(x, sep, n) {
    var = unlist(strsplit(x, split=sep))
    length(var) = n
    return(var)
  }

  x = do.call(cbind, lapply(x, .splitvar, sep=sep, n=n))
  x = apply(x, 1, paste, collapse=sep)
  out = read.csv(text=x, sep=sep, header=header, ...)
  return(out)
}

removeFactors = function(df){

}

getStringElement = function(string, delimiter, element){
  return(unlist(sapply(strsplit(string,delimiter),"[",element)))
}

dfcols.tochar = function(df){
  return(data.frame(lapply(df, as.character), stringsAsFactors=FALSE))
}
