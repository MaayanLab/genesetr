pairwiseContingTables = function(lib1, lib2,background = NULL){

  if(is.null(background)){
    background = length(unique(c(getGenes(lib1),getGenes(lib2))))
  }

  x = genesetr::pairwiseIntersects(lib1,lib2)
  x = reshape2::melt(x)
  colnames(x) = c("set1","set2","intersect")
  x$intersect = as.integer(x$intersect)
  x$set1 = as.character(x$set1)
  x$set2 = as.character(x$set2)

  n1 = plyr::ddply(toLongDF(lib1),.(set_name),nrow)
  n2 = plyr::ddply(toLongDF(lib2),.(set_name),nrow)

  x$len_set1 = as.integer(n1[match(x$set1,n1$set_name),2])
  x$len_set2 = as.integer(n2[match(x$set2,n2$set_name),2])
  x$background = as.integer(rep(background,nrow(x)))

  x$a = as.integer(x$intersect)
  x$b = x$len_set2 - x$intersect
  x$c = x$len_set1 - x$intersect
  x$d = x$background - x$len_set1 - x$len_set2 + x$intersect

  return(x)
}

pairwiseIntersects = function(lib1, lib2){
  lib1 = toLongDF(lib1)
  lib2 = toLongDF(lib2)
  genes = intersect(getGenes(lib1),getGenes(lib2))
  lib1 = lib1[lib1$gene %in% genes,]
  lib2 = lib2[lib2$gene %in% genes,]
  return(as.matrix(crossprod(table(lib1),table(lib2))))
}

getGenes = function(x){
  if(isMat(x)) return(rownames(x))
  else if(isLongDF(x)) return(unique(x[,1]))
  else if(isList(x)) return(unique(unlist(x)))
  else error("Unrecognized library format.")
}
