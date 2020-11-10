#  -----------------------------------------------------------------------------
#' PrintMlstSpecies
#  2020/11/07
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#' @description  Download and format a MLST database from PubMLST.
#' @param species If looking for a specific species, prints all of the
#'   availible databases that start with the same letter as the input
#' @return NULL
#' @export
#  -----------------------------------------------------------------------------

PrintMlstSpecies = function( mlstDbs, species )
{
  # Find the species of interest
  speciesNames = sapply( 1:length( mlstDbs ), function(i) mlstDbs[[i]][[1]] )
  
  if ( missing( mlstDbs ) )
  {
    # Read in the xml file with information on the mlst databases
    mlstDbs = xml2::as_list(
      xml2::read_xml( "http://pubmlst.org/data/dbases.xml" )
      )[[1]]
  }
  
  
  if ( missing( species ) )
  {
    cat( "Here are the availible species:\n" )
  } else {
    
    firstChar = substr( species, 1, 1 )
    
    isCharMatch = 
      grepl( paste0( '^', firstChar ), speciesNames, ignore.case = TRUE )
    speciesNames = speciesNames[ isCharMatch ]
    cat( 
      "Here are the availible species that start with ",
      firstChar, ":\n", sep = ''
      ) 
  }
  
  cat( paste( "  -- ", speciesNames ), sep = '' )

}

# ------------------------------------------------------------------------------
