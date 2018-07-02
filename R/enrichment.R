#' Fast version of Fisher Exact Test
#'
#'
#' Quickly compute FET p-values for n 2x2 contingency tables T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}
#'
#' @param a A vector of values where a_{i} corresponds to t_{i}
#' @param b A vector of values where b_{i} corresponds to t_{i}
#' @param c A vector of values where c_{i} corresponds to t_{i}
#' @param d A vector of values where d_{i} corresponds to t_{i}
#' @return A vector of p-values P={p_{1}, p_{2},...,p_{n-1},p_{n}} where p_i corresponds to the FET result of t_i.
#' @examples
#' T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}, with each t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'
#' Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  5     6        setD  20     45
#'¬setB  123   500       ¬setC  10    100     ¬setD  60     1000
#'
#'
#'
#' a=c(3, 5, 20)
#' b=c(100, 6, 45)
#' c=c(123, 10, 60)
#' d=c(500, 100, 1000)
#' pvals = fastFET(a,b,c,d)
fastFET = function(a, b, c, d, alternative = "greater"){
  checkContingTableVals(a,b,c,d)

  l_a = length(a)
  #if a,b,c,d are length 1 use FET from standard R distribution
  if(l_a == 1) return(fisher.test(matrix(c(a,c,b,d),ncol = 2),
    alternative = alternative)$p.value)

  #populate factorial look-up vector
  f = numeric(length = max(rowSums(cbind(a,b,c,d))))
  f[1] = log(1)
  for(i in 2:length(f)){
    f[i] = f[i-1] + log(i)
  }

  pval = cbind(numeric(length = l_a),numeric(length = l_a))
  coeff = 1
  min = cbind(pmin(b,c),rep(NA,l_a))
  sides = 1

  if(alternative == "less"){
    coeff = -1
    min = cbind(pmin(a,d),rep(NA,l_a))
  }

  if(alternative == "two.sided"){
    coeff = c(-1,1)
    min = cbind(pmin(a,d),pmin(b,c))
    sides = 2
  }
for(j in 1:sides){
  for(i in 1:l_a){
    temp_a = a[i]:(a[i]+coeff[j]*min[i,j])
    temp_b = b[i]:(b[i]-coeff[j]*min[i,j])
    temp_c = c[i]:(c[i]-coeff[j]*min[i,j])
    temp_d = d[i]:(d[i]+coeff[j]*min[i,j])

    #replace 0s with 1 since R indexing starts at 1
    temp_a[temp_a==0] = 1
    temp_b[temp_b==0] = 1
    temp_c[temp_c==0] = 1
    temp_d[temp_d==0] = 1

    invariant = (f[a[i]+b[i]] + f[c[i]+d[i]] + f[a[i]+c[i]] + f[b[i]+d[i]])-
      f[a[i]+b[i]+c[i]+d[i]]#n
    pval[i,j] = sum(exp(invariant -
        (f[temp_a] + f[temp_b] + f[temp_c] + f[temp_d])))
  }}
  #computed as in fisher.exact in R stats package
  #unclear whether this is appropriate -Ali 5/21/2018
  if(alternative == "two.sided"){
    pval[,1] = pmin(pval[,1],pval[,2])
  }
  return(pval[,1])
}

#' Odds Ratio for Multiple Contigency Tables
#'
#'
#' Compute odds ratios for n 2x2 contingency tables T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}
#'
#' @param a A vector of values where a_{i} corresponds to t_{i}
#' @param b A vector of values where b_{i} corresponds to t_{i}
#' @param c A vector of values where c_{i} corresponds to t_{i}
#' @param d A vector of values where d_{i} corresponds to t_{i}
#' @return A vector of odds ratios R={r_{1}, r_{2},...,r_{n-1},r_{n}} where r_i corresponds to the odds ratio of t_i.
#' @examples
#' T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}, with each t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'
#' Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  5     6        setD  20     45
#'¬setB  123   500       ¬setC  10    100     ¬setD  60     1000
#'
#'
#'
#' a=c(3, 5, 20)
#' b=c(100, 6, 45)
#' c=c(123, 10, 60)
#' d=c(500, 100, 1000)
#' or = oddsRatio(a,b,c,d)
oddsRatio = function(a,b,c,d){
  checkContingTableVals(a,b,c,d)
  return((a/b)/(c/d))
}

#' Jaccard Index for Multiple Contigency Tables
#'
#'
#' Compute jaccard index for n 2x2 contingency tables T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}
#'
#' @param a A vector of values where a_{i} corresponds to t_{i}
#' @param b A vector of values where b_{i} corresponds to t_{i}
#' @param c A vector of values where c_{i} corresponds to t_{i}
#' @param d A vector of values where d_{i} corresponds to t_{i}
#' @return A vector of Jaccard indices j={j_{1}, j_{2},...,j_{n-1},j_{n}} where j_i corresponds to the jaccard index of t_i.
#' @examples
#' T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}, with each t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'
#' Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  5     6        setD  20     45
#'¬setB  123   500       ¬setC  10    100     ¬setD  60     1000
#'
#'
#'
#' a=c(3, 5, 20)
#' b=c(100, 6, 45)
#' c=c(123, 10, 60)
#' d=c(500, 100, 1000)
#' ji = jaccardInd(a,b,c,d)
jaccardInd = function(a,b,c,d){
  checkContingTableVals(a,b,c,d)
  j = a/(a+b+c)
}

#' Jaccard Distance for Multiple Contigency Tables
#'
#'
#' Compute jaccard distance for n 2x2 contingency tables T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}
#'
#' @param a A vector of values where a_{i} corresponds to t_{i}
#' @param b A vector of values where b_{i} corresponds to t_{i}
#' @param c A vector of values where c_{i} corresponds to t_{i}
#' @param d A vector of values where d_{i} corresponds to t_{i}
#' @return A vector of Jaccard distances j={j_{1}, j_{2},...,j_{n-1},j_{n}} where j_i corresponds to the jaccard distance of t_i.
#' @examples
#' T = \{t_{1},t_{2},...,t_{n-1},t_{n}\}, with each t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'
#' Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  5     6        setD  20     45
#'¬setB  123   500       ¬setC  10    100     ¬setD  60     1000
#'
#'
#'
#' a=c(3, 5, 20)
#' b=c(100, 6, 45)
#' c=c(123, 10, 60)
#' d=c(500, 100, 1000)
#' jd = jaccardDist(a,b,c,d)
jaccardDist = function(a,b,c,d){
  return(1-jaccardInd(a,b,c,d))
}

checkContingTableVals = function(a,b,c,d){
  l_a = length(a)
  l_b = length(b)
  l_c = length(c)
  l_d = length(d)

  if(any(!rep(l_a,3)==c(l_b,l_c,l_d))){
    stop("Vectors a,b,c and d must all be of equal length.")
  }
  if(!is.numeric(c(a, b, c, d))){
    stop("All elements of a,b,c, and d must be numeric.")
  }
  if(any(is.na(c(a,b,c,d)))){
    warning("NA elements will produce NA result.")
  }
}
