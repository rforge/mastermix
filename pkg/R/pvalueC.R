#' @title pvalueC
#' @author Oyvind Bleka <Oyvind.Bleka.at.fhi.no>
#' @description Calculates Tippet-p-value in C.
#' @export
#' @details This function calculates restricted sum for all values in LR.table greater than LR.observed and summing corresponding values in P.value.
#' @param LR.suspect Observed LR to restrict sum on.
#' @param LR.table matrix of all possible LR.
#' @param P.table corresponding values to sum assigned in LR.table.
#' @return Restriced sum from calculation.


pvalueC <- function( LR.suspect, LR.table, P.table ){
 prepare.input <- function( LR.table, P.table ) {
 #
 # Converting the tables of likelihood ratios and probabilities to lists,
 # keeping only the unique LR-scores for each marker, and sorting everything
 # in proper order.
 #
 # Input:
 # LR.table  - Matrix of likelihood ratio scores for all genotypes of all markers. 
 #             Each row correspond to a marker, and it may contain 0's (MxG matrix)
 # P.table   - Matrix of population probabilities for all genotypes of all markers.
 #             Each row correspond to a marker, and it may contain 0's (MxG matrix)
 #
 # The function returns a list containing two lists, LR.list and P.list, containing the unique 
 # non-zero values in LR.table and P.table. For each marker the elements are sorted in descending
 # order by LR-score.
 #
  M <- nrow( LR.table ) # The number of markers
  LR.list <- P.list <- vector( "list", M )
  LR.max <- apply( LR.table, 1, max )
  
  # Smart sorting of markers
  ix <- order( LR.max, decreasing=T )
  LR.table <- LR.table[ix,]
  P.table <- P.table[ix,]
  LR.max <- LR.max[ix]
  
  for( marker in 1:M ){
    # Finding unique likelihood ratios and summing probabilities accordingly
    lr <- unique( LR.table[marker,] )
    pr <- tapply( P.table[marker,], match( LR.table[marker,], lr ), sum )
    # Discarding 0's
    ix <- which( lr > 0 )
    lr <- lr[ix]
    pr <- pr[ix]
    # Sorting in descending order by LR
    ix <- order( lr, decreasing=T )
    LR.list[[marker]] <- lr[ix]
    P.list[[marker]] <- pr[ix]
  }
  Inflation.factor <- cumprod( LR.max[M:1] )[M:1]
  return( list( LR.list=LR.list, P.list=P.list, Inflation.factor=Inflation.factor ) )
 }
 #
 # Computes the p-value for LR.suspect
 #
 # Input:
 # LR.suspect  - The likelihood ratio for the suspect genotype profile (1x1 positive value)
 # LR.table    - The pre-computed likelihood ratios for every genotype of every marker (MxG matrix)
 # P.table     - The population probabilities for every genotype of every marker (MxG matrix)
  
  # The number of markers
  M <- nrow( LR.table )
  
  # Preparing input
  prepped <- prepare.input( LR.table, P.table )
  LR.list <- prepped$LR.list
  Plist  <- as.numeric(unlist(prepped$P.list))
  inflationfactor <- as.numeric(prepped$Inflation.factor)
  
  # The number of genotypes for each marker
  G <- as.integer(sapply( LR.list, length ))

  # Initial values
  pvalue <- as.numeric(1)
  LRsuspect  = as.numeric(LR.suspect)
  LRlist <- as.numeric(unlist(LR.list))

  p.value <- .C("prodrecfunction",pvalue, LRsuspect, LRlist,Plist,M,G,inflationfactor,PACKAGE="mastermix")[[1]]
  return( p.value )
}


