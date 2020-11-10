#  -----------------------------------------------------------------------------
#  ParseBlastData
#  2020/11/09
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#  Sort the blast output data to return a vector with the allele ids in the 
#  order of the key 
#  -----------------------------------------------------------------------------

ParseBlastData = function( mlstGenes, blastData )
{
  mlstResults = sapply( mlstGenes, function(x)
  {
    isMatch = blastData$geneId == x
    if ( !TRUE %in% isMatch ) return( NA )
    return( blastData$alleleId[ isMatch ] )
  })
  
  return( mlstResults )
}

#  -----------------------------------------------------------------------------