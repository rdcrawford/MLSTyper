# ------------------------------------------------------------------------------
# GetBestMlstMatches
# 2018/07/20
# Ryan D. Crawford
# ------------------------------------------------------------------------------
# For each genome, find the best match for the MLST profile. Return the 
# results as a character matrix. The last column with the mutations in the
# mismatch genes. 
# ------------------------------------------------------------------------------

FindMlst = function(
  mlstKey, # The reference key with the mlst types
  blList   # List with the blast data
  )
{
  # Get the gene ids 
  mlstGenes = colnames( mlstKey )[ 2:ncol( mlstKey ) ]
  
  # Make a matrix with the mlst results
  mlstResults = t( future.apply::future_sapply( 1:length( blList ),
    function(i) ParseBlastData( mlstGenes, blList[[i]] )
    ))
  
  # Find the best match for the mlst profiles
  mlstMatches = future.apply::future_sapply( 1:nrow( mlstResults ),
    function(i) CalcMlstMatch( mlstKey, mlstResults[ i, ] )
    )
  
  # Add the ST matches to the output data
  mlstResults = cbind( mlstMatches, mlstResults )
  
  # Find if there are any mismatches in the blast data
  isPerfectMatch = future.apply::future_sapply( 1:length( blList ),
    function(i) sum( blList[[i]]$nMismatches ) == 0
    )
  
  # If there are any instances where there is not a perfect match,
  # add the genes where there is a mismatch to the output data
  if ( FALSE %in% isPerfectMatch )
  {
    mismatchAlleles = vector( "character", length( blList ) )
    
    for ( i in which( !isPerfectMatch ) )
    {
      isMismatch = blList[[i]]$nMismatches != 0
      mismatchAlleles[i] = paste( paste( 
        blList[[i]][ isMismatch, 1], 
        blList[[i]][ isMismatch, 2],
        blList[[i]][ isMismatch, 3],
        sep = ',' 
        ), collapse = ";" )
      
      # If there is not already an astrix on the resul, add one
      if ( !grepl( '*', mlstResults[ i, 1 ], fixed = TRUE ) )
         mlstResults[ i, 1 ] = paste0( '*',   mlstResults[ i, 1 ] )
    }
    mlstResults = cbind( mlstResults, mismatchAlleles )
  }
  
  colnames( mlstResults )[1] = "ST"
  row.names( mlstResults ) = NULL
  return( mlstResults )
}
# ------------------------------------------------------------------------------