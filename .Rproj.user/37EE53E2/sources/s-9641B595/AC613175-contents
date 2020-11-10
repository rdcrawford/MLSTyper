
#  -----------------------------------------------------------------------------
#' GetMlstPath
#  2020/11/07
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#' @description Get the path to the MLST database if the species has not 
#'   been downloaded, download the species.
#' @param scheme Organisms to look up
#' @return String with the path to the mlst data
#' @export
#  -----------------------------------------------------------------------------

GetMlstPath = function( species )
{
  mlstDirs = system( 
    paste0("ls -d ", .libPaths()[1], "/MLSTyper/data/*" ),
    intern = TRUE 
    )
  idx = grep( gsub( ' ', '_', species ), mlstDirs, ignore.case = TRUE )
  if ( length( idx ) == 1 ) return( paste0( mlstDirs[ idx ], '/' ) )
  return( DownloadMlstData( species ) )
}

#  -----------------------------------------------------------------------------