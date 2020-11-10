#  -----------------------------------------------------------------------------
#' DownloadMlstData
#  2020/11/07
#  Ryan D. Crawford
#  -----------------------------------------------------------------------------
# species = "Klebsiella"
# speciesspecies = "Klebsiella pneumoniae"
#  -----------------------------------------------------------------------------

DownloadMlstData = function( species )
{
  # Read in the xml file with information on the mlst databases
  mlstDbs = xml2::as_list(
    xml2::read_xml( "http://pubmlst.org/data/dbases.xml" )
    )[[1]]
  
  # Find the species of interest
  isSpecies = sapply( 1:length( mlstDbs ),
    function(i) grepl( species, mlstDbs[[i]][[1]] )
    )
  
  # Check that the species was found
  numHits = sum( isSpecies )
  if ( numHits == 0 )
  {
    stop( paste( species, "was not found...") )
  } else if ( numHits > 1 ) {
    hits = sapply( which( isSpecies ),
      function(i) mlstDbs[[i]][[1]]
      )
    errStr = paste(
      "There was more than one match. Select one of:\n",
      paste(  paste( "  -- ", hits), collapse = '' ),
      sep = ''
      )
    stop( errStr )
  }

  # Make a directory to store this species data
  dataDir = 
    paste0( .libPaths()[1], "/MLSTyper/data/", gsub( ' ', '_', species ), '/' )
  if ( !file.exists( dataDir ) ) system( paste( "mkdir", dataDir ) )
    
  # Download the key
  speciesData = mlstDbs[[ which( isSpecies ) ]][[2]]
  profileUrl  = speciesData[[1]][[2]][[1]][[1]]
  cmd = paste0( "wget --output-document=", dataDir, "key.tsv ", profileUrl )
  system( cmd )

  # Download the fasta files
  numGenes = length( speciesData[[1]][[3]] )
  geneFas = sapply ( seq( numGenes ), function(i)
  {
    geneId  = trimws( speciesData[[1]][[3]][[i]][[1]][[1]] )
    geneFa  = paste0( dataDir, geneId, ".fasta " )
    geneUrl = speciesData[[1]][[3]][[i]][[2]][[1]]
    cmd     = paste( "wget --output-document", geneFa, geneUrl )
    system( cmd )
    return( geneFa )
  })
  
  # Make the blast database
  genesFile = paste0( dataDir, "mlst_genes.fasta" )
  system( paste0( "cat ",  paste( geneFas, collapse = ' ' ), '>', genesFile ) )
  system( paste( "makeblastdb -in", genesFile, "-dbtype nucl") )
  return( dataDir )
}

# ------------------------------------------------------------------------------
