#  -----------------------------------------------------------------------------
#' GetMlst
#  2020/11/07
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
#' @description Look up the mlst type for an organisms of interest
#' @param fastaDir Directory containing the fasta files. By default, this uses
#'   all of the files in the directory. If there are additional files in this 
#'   directory "fastaExt" can also be supplied to select only the appropriate 
#'   files. This argument is incompatible with "fastaFiles," which can be used
#'   to supply the paths to the fasta files as a character vector.
#' @param fastaExt Optional file extension for the fasta files.
#' @param fastaFiles A character vector with the paths to the the fasta files 
#'   for each genome can be input. Incompatible with "fasta dir."
#' @param outDir Directory to write the alignments, other output files, 
#'   and temp files created by cd-hit and mafft. Defaults to the current 
#'   working directory.
#' @return A data.frame with the mlst types of the organisms of interests
#' @export
#  -----------------------------------------------------------------------------

GetMlst = function(
  species,    # species
  fastaDir,   # Directory containing the fasta files
  fastaFiles, # Fasta files for the input genomes
  fastaExt,   # Extension used on the fasta files
  genomeIds,  # Vector of genomes to
  outFile     # Optional. Directory to write the output files
  )
{
  startTime = Sys.time() # Start the timer
  # ---- Parse the input arguments ---------------------------------------------

  # Make sure that there aren't incompatible arguments with respect to 
  # fastaFiles and fastaDir
  if ( ( !missing(fastaFiles) && !missing(fastaDir) ) ||  
    ( missing(fastaFiles) && missing(fastaDir) ) 
    )
  {
    errorMsg = paste0(
      "Input arguments must contain one of \"fastaFiles\" (",
      "character vector with paths to fasta files)",
      " or \"fastaDir\"(directory containing fasta files)"
      )
    stop( errorMsg )
    
  } else if ( missing( fastaFiles )  ) {

    # Get the paths to the fasta files 
    fastaFiles = GetFilePaths( fastaDir, fastaExt )
  }
    
  if ( missing( species ) ) stop( "Must input a species!" )

  # If no genome ids were created, make some using the fasta files
  if ( missing( genomeIds ) )
  {
    genomeIds = 
      sapply( fastaFiles, ExtractGenomeNameFromPath, USE.NAMES = FALSE )
  }
  
  # ---- Find the target genes using the gene ids ------------------------------

  cat(
    "\n\nGetting MLST:\n",
    "  -- ", length( fastaFiles ), " genomes were input\n",
    sep = ''
    )
  
  # Blast the fasta files against the mlst data-base
  cat("\nStep 1: getting path to the MLST data-base\n")
  mlstPath = GetMlstPath( species )
  blDb     = paste0( mlstPath, "mlst_genes.fasta" )
  mlstKey  = ReadInKey( mlstPath )
  stepTime = GetSplit( startTime )
  
  # Find the best match for each allele in the mlst genes
  cat("\nStep 2: blasting the fasta files\n")
  blList = future.apply::future_sapply( fastaFiles,
    function(x) list( RunBlast( x, blDb ) )
    )
   stepTime = GetSplit( stepTime )
  
  # Identify orthologous genes with cd-hit
  cat("\nStep 3: finding mlst matches\n")
  mlstData = cbind( genomeIds, FindMlst( mlstKey, blList ) )
  stepTime = GetSplit( stepTime )
  
  # If an output file was requested, write the table to the output
  if ( !missing( outFile ) )
  {  
    # Identify orthologous genes with cd-hit
    cat("\nStep 4: writing output file\n")
    
    # Write the data on the alignment genes
    write.table( 
      mlstData, 
      file      = outFile, 
      append    = FALSE, 
      sep       = "\t",
      row.names = FALSE,
      col.names = TRUE
      )  
    stepTime = GetSplit( stepTime )
  }
  
  cat( "\nRun complete\n" )
  endTime = GetSplit( startTime )
  cat( "\n\n" )
  
  return( mlstData )
}

# ------------------------------------------------------------------------------