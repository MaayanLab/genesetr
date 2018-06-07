#' Contingency Tables
#'
#' 2x2 contingency tables
#'
#' @param lib1 A library of gene sets in List, LongDF or Mat format.
#' @param lib2 A library of gene sets in List, LongDF or Mat format.
#' @param background The size of the universe of genes under consideration. User may set to a constant. Default is NULL which results in background being set to the size of the union of the unique genes of lib1 and lib2.
#' @return A data frame with the following fields:
#' \describe{
#'   \item{a}{x}
#'   \item{b}{x}
#'   \item{c}{x}
#'   \item{d}{x}
#'   \item{background}{x}
#'   \item{len_set1}{x}
#'   \item{len_set2}{x}
#'   ...
#' }
#' @examples
#' #load example gene set libraries
#' data("encode_example","chea_example")
#' conting_tab = pairwiseContingTables(lib1, lib2,background= 30000)
#' conting_tab = pairwiseContingTables(lib1, lib2)
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

#' All Pairwise Intersections Between Two Libraries
#' Quickly computes the size of the intersection between all pairs of sets within a library.
#'
#' @param lib1 A library of gene sets in List, LongDF or Mat format.
#' @param lib2 A library of gene sets in List, LongDF or Mat format.
#' @return A matrix with lib1 gene set names along the rows and lib2 gene set names on the columns populated with the sizes of the intersections.
#' @examples
#' #load example gene set libraries
#' data("encode_example", "chea_example")
#' intersects = pairwiseIntersects(encode_example, chea_example)
pairwiseIntersects = function(lib1, lib2){
  lib1 = toLongDF(lib1)
  lib2 = toLongDF(lib2)
  genes = intersect(getGenes(lib1),getGenes(lib2))
  lib1 = lib1[lib1$gene %in% genes,]
  lib2 = lib2[lib2$gene %in% genes,]
  return(as.matrix(crossprod(table(lib1),table(lib2))))
}

#' Unique Genes in a Library
#'
#'Gets unique genes in a library
#'
#' @param x A library of sets in List, LongDF or Mat format.
#' @return A character vector of all unique gene names in a library.
#' @examples
#' #load example gene set libraries
#' data("encode_example")
#' unique_encode_genes = getGenes(encode_example)
getGenes = function(x){
  if(isMat(x)) return(rownames(x))
  else if(isLongDF(x)) return(unique(x[,1]))
  else if(isList(x)) return(unique(unlist(x)))
  else error("Unrecognized library format.")
}
