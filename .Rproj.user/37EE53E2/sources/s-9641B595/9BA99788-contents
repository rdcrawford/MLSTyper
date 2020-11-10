#  -----------------------------------------------------------------------------
#' ReadInKey
#  2020/11/09
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#' @description Look up the mlst type for an organisms of interest
#' @param mlstPath Path containing the mlst data
#' @return A character matirx with the mlst types
#  -----------------------------------------------------------------------------

ReadInKey = function( mlstPath )
{
  # Read in the tsv
  mlstData = 
    strsplit( readLines( paste0( mlstPath, "key.tsv" ) ), split = '\t' )
  
  # Transform the data to a character matirx
  nTypes = length( mlstData )
  mlstKey = matrix(
    unlist( mlstData[2:nTypes] ),
    nrow = nTypes - 1,
    byrow = TRUE
    )
  
  # Set the columns names to the gene names
  colnames( mlstKey ) = sapply( mlstData[[1]], function(x)
  {
    comps = strsplit(x, split = '_')[[1]]
    return( comps[ length( comps ) ] )
  }, USE.NAMES = FALSE )
  
  return(  mlstKey )
}

# ------------------------------------------------------------------------------