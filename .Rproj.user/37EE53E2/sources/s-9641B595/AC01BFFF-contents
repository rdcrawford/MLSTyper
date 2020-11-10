#include <Rcpp.h>
#include <Rcpp.h>
#include <sstream>
#include "BlastData.h"
using namespace Rcpp;

//[[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------
// Parse Blast Data 
// Ryan D. Crawford
// 2019/10/10
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

BlastData::BlastData(
  const std::vector<std::string> &fastaPaths,
  const std::string &blastDB 
  )
{
  // blast each genome against the 
  for (std::string fa : fastaPaths) 
  {
    std::string blastCmd = createBlastCmd( fa, blastDB );
    runBlast( blastCmd );
  }
}

// Create the blast command to run 
std::string BlastData::createBlastCmd( std::string faFile, std::string dbPath )
{
  std::string blastCmd;
  blastCmd =  "blastn ";
  blastCmd += "-outfmt \"6 sseqid nident slen\"";
  blastCmd += " -query ";
  blastCmd +=  faFile;
  blastCmd += " -db ";
  blastCmd += dbPath;
  return blastCmd;
}

void BlastData::parseMlstHeader(
  std::string header, std::string &gene, int &allele
  )
{
  // Find the positions of the last and next to last underscores
  auto lastUs       = header.find_last_of( '_' ) + 1;
  auto nextToLastUs = header.rfind( '_', lastUs - 2 ) + 1;
  
  // Get the gene between the last underscore and the end of the gene
  gene   = header.substr( nextToLastUs, lastUs - nextToLastUs - 1 );
  
  // Get the integer corresponding to the unique allele between the last
  // underscore and the end of the string
  allele = stoi( header.substr( lastUs, header.length() - lastUs + 1 ) );
}

DataFrame BlastData::createDataFrame(  )
{
  DataFrame blastData = DataFrame::create( 
    Named("geneId")  = geneIds,
    _["alleleId"]    = alleleIds,
    _["nMismatches"] = nMismatches
    );
  return blastData;
}

bool BlastData::parseBlastEntry(
  std::stringstream &ss, std::string &geneId, int &alleleId, int &nMismatch
  )
{
  // std::string s = ss.str();
  // Rcout << s << std::endl;
  // The variables in the blast results
  std::string sAlgnId;
  int nident;
  int slen;
  
  // Find the gene and allele ID for the blast hit
  ss >> sAlgnId;
  parseMlstHeader( sAlgnId, geneId, alleleId );
  
  // Calculate the number of mismatches
  ss >> nident ;
  ss >> slen;
  nMismatch = slen - nident;
  
  // If the gene is present, return false to indicate that this
  // blast entry should not be added to the output data
  if ( FindGene( geneId, alleleId, nMismatch ) ) return false;
  return true;
}

bool BlastData::FindGene( 
  const std::string &geneId, const int &alleleId, const int &nMismatch
  )
{
  // Find element in a vector 
  auto it = std::find( geneIds.begin(), geneIds.end(), geneId ); 
  if ( it == geneIds.end() ) return false;
  
  // Check that this is the minimium number of mismatches
  auto idx = std::distance( geneIds.begin(), it );
  if ( nMismatch < nMismatches[ idx ] )
  {
    // Reassign the number of mis matches and allele id
    nMismatches[ idx ] = nMismatch;
    alleleIds[ idx ]   = alleleId;
  }
  return true;
}

void BlastData::runBlast( const std::string &blastCmd )
{
  // Create a buffer to store the blast results
  std::array<char, 128> buffer;
  std::unique_ptr< FILE, decltype( &pclose ) > pipe(
    popen( blastCmd.c_str(), "r" ), pclose
    );
  
  // If running the command failed throw an errow 
  if ( !pipe ) Rcpp::stop( "blastn command failed!" );
  
  // Append the results of each alignment to the vectors composing the 
  // blast data
  while ( fgets( buffer.data(), buffer.size(), pipe.get() ) != nullptr )
  {
    // For each alignment, assign the data to a string string stream
    std::stringstream ss;
    ss << buffer.data();

    std::string geneId;
    int         alleleId;
    int         nMismatch;
    
    // If this is a new allele inthe blast data, add it to the 
    // vectors to output 
    if ( parseBlastEntry( ss, geneId, alleleId, nMismatch ) )
    {
      geneIds.push_back( geneId );
      alleleIds.push_back( alleleId );
      nMismatches.push_back( nMismatch );
    }
    
    ss.clear();
  }
}

// -----------------------------------------------------------------------------
