toMat = function(x){
  if(isMat(x)){
    warning("Input x is already in library matrix format. Returning x.")
    return(x)
  }else if(isList(x)){
    x = toLongDF(x)
  }else if(!isLongDF(x)){
    error("Library is in unrecognized format.
      Must be library list or library long DF format. See documentation.")
  }
  x$val = 1
  x = reshape2::dcast(x,gene~set_name,fill = 0, drop = F,value.var = "val")
  rownames(x) = x[,1]
  x = x[,-1]
  return(as.matrix(x))
}

toLongDF = function(x){
  if(isMat(x)){
    x = toList(x)
  }
  if(isList(x)){
    x <- stack(x)
    x <- data.frame(lapply(x, as.character), stringsAsFactors=FALSE)
    colnames(x) <- c("gene","set_name")
    return(x)
  }
  else if(isLongDF(x)){
    warning("Input x is already in long df format. Returning x.")
    return(x)
  }
  else{
    error("Library is in unrecognized format.
      Must be library list or library matrix format. See documentation.")
  }
}

toList = function(x){
  if(genesetr::isMat(x)||genesetr::isLongDF(x)){
    if(genesetr::isMat(x)){
      x = as.data.frame(x)
      x$gene = rownames(x)
      x = reshape2::melt(x,id.vars = "gene")
      x = x[x$value>0,]
      x = data.frame(gene = x$gene, set_name = x$variable,stringsAsFactors = F)
    }
    x = plyr::dlply(x,plyr::.(set_name),function(sub){return(sub$gene)})
    return(x)
  }else if(genesetr::isList(x)){
    warning("Input x is already in list format. Returning x.")
    return(x)
  }else{
    error("Library is in unrecognized format.
      Must be library long df or library matrix format. See documentation.")
  }
}

isMat = function(x){
  return(is.matrix(x))
}

isList = function(x){
 return(inherits(x,"list"))
}

isLongDF = function(x){
  return(inherits(x,"data.frame"))
}

removeSetWeights = function(x){
  return(unlist(sapply(strsplit(x,","),"[",1)))
}

removeLibWeights = function(x){
  return(lapply(x,removeSetWeights))
}

isWeightedLib = function(x){
  if(genesetr::isList(x)){
    return(all(unlist(lapply(x,function(set){
      return(all(grepl(",",set)))
    }))))
  }
}
