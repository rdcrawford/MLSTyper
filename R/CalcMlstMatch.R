# ------------------------------------------------------------------------------
# GetBestMlstMatches
# 2018/07/20
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# Iterate over the MLST key until an exact match is found. If no exact match
# is found return the best match followed by an astrix
# ------------------------------------------------------------------------------

CalcMlstMatch = function(
  mlstKey,    # The reference key with the mlst types
  mlstProfile # List with the blast data
  )
{
  matchScore = 0
  bestMatch  = 0
  idx        = 1
  geneRange  = 2:ncol( mlstKey )
  
  while ( idx < nrow( mlstKey ) )
  {
    curSrcore = sum( sapply( geneRange, function(j) 
    {
      if ( is.na( mlstProfile[ j - 1 ] ) ) return( FALSE )
      return( mlstKey[ idx, j ] == mlstProfile[ j - 1 ] )
    }))
   
    if ( curSrcore > matchScore )
    {
      matchScore = curSrcore
      bestMatch  = mlstKey[ idx, 1 ]
      if ( matchScore == length( mlstProfile ) ) break
    }
    
    idx = idx + 1 
  }
  
  if ( matchScore == length( mlstProfile ) ) return( bestMatch )
  return(  paste0( bestMatch, '*' ) )
}

# ------------------------------------------------------------------------------
