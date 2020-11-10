#  -----------------------------------------------------------------------------
#  ParseBlastData
#  2020/11/09
#  Ryan D. Crawford
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

# ------------------------------------------------------------------------------