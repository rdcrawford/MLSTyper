#include <Rcpp.h>
using namespace Rcpp;

// -----------------------------------------------------------------------------
// Parse Blast Data 
// Ryan D. Crawford
// 2019/10/10
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

//[[Rcpp::plugins(cpp11)]]

// ---- Class definition -------------------------------------------------------

class BlastData
{
public:
  // Ctor: input a vector of 
  BlastData(const std::vector<std::string> &fastaPaths, 
    const std::string &blastDB );
  
  // Dtor
  ~BlastData(){ ; }
  
  // Create a dataframe for this blast data
  DataFrame createDataFrame();
  
  // Find Best blast aligns
  
private:
  
  // Columns of the output data-frame
  std::vector< std::string > geneIds;
  std::vector< int >         alleleIds;
  std::vector< int >         nMismatches;
  
  // Parse the information on a blast alignment 
  bool parseBlastEntry( std::stringstream &ss, std::string &geneId, 
    int &alleleId, int &nMismatch );
  
  // Parse the header to get the 
  void parseMlstHeader( std::string header, std::string &gene, int &allele );
  
  // Run the blast command and append the results to the data frame
  void runBlast( const std::string &blastCmd );
  
  // Create the blast command to run 
  std::string createBlastCmd( std::string faFile, std::string dbPath);
  
  // Determine if the information on this gene has already been recorded
  bool FindGene( const std::string &geneId, const int &alleleId, 
    const int &nMismatch );
};

// -----------------------------------------------------------------------------
