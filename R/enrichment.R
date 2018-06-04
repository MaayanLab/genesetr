#' Fast version of Fisher Exact Test
#'
#'
#' Quickly compute FET p-values for n 2x2 contingency tables T = \{t_1,t_2,...,t_n-1,t_n\}, with each t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'
#' @param a A vector of values where a_i corresponds
#' @param b A vector of values where b_i corresponds to
#' @param c A vector of values where c_i corresponds to
#' @param d A vector of values where d_i corresponds to
#' @return A vector of p-values P={p_1, p_2,...,p_n-1,p_n} where p_i corresponds to the FET result of t_i.
#' @examples
#' Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  a     b        setD  a     b
#'¬setB  123   500       ¬setC  c     d       ¬setD  c     d
#'
#'
#'
#' a=c()
#' b=c()
#' c=c()
#' d=c()
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

#' Odds Ratio
#'
#'
#' For a set of n 2x2 contingency tables T = \{t_1,t_2,...,t_n-1,t_n\}, with each contingency table t_i of the form:
#'            setAn ¬setA
#'        setB  a     b | a+b
#'       ¬setB  c     d | c+d
#'          ------------|-----
#'             a+c   b+d| a+b+c+d
#'Computes an odds ratio ri for each ci.
#'
#' @param a A vector of numbers a = {a1, a2,...,an-1,an} where a_i corresponds to continency table ti
#' @param b A vector of numbers b = {b1, b2,...,bn-1,bn} where b_i corresponds to continency table ti
#' @param c A vector of numbers c = {c1, c2,...,cn-1,cn}where c_i corresponds to continency table ti
#' @param d A vector of numbers d = {d1, d2,...,dn-1,dn}where d_i corresponds to continency table ti
#' @return A vector of odds ratios R={r_1, r_2,...,r_n-1,r_n} where r_i corresponds to the FET result of t_i.
#' @examples
#' Given three contingency tables:
#'      setA  ¬setA          setA  ¬setA         setB  ¬setB
#' setB  3     100        setC  a     b        setD  a     b
#'¬setB  123   500       ¬setC  c     d       ¬setD  c     d
#'
#'
#'
#' a=c()
#' b=c()
#' c=c()
#' d=c()
#' pvals = fastFET(a,b,c,d)
oddsRatio = function(a,b,c,d){
  checkContingTableVals(a,b,c,d)
  return(b*c/a*d)
}

jaccardInd = function(a,b,c,d){
  checkContingTableVals(a,b,c,d)
  j = a/(a+b+c)
}

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
