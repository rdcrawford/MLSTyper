#include <Rcpp.h>
#include "BlastData.h"
using namespace Rcpp;

// -----------------------------------------------------------------------------
// RunBlast
// Ryan D. Crawford
// 2019/11/06
// -----------------------------------------------------------------------------
//
// -----------------------------------------------------------------------------

//[[Rcpp::plugins(cpp11)]]

// -----------------------------------------------------------------------------

// [[Rcpp::export]]
DataFrame RunBlast(
  const std::vector<std::string> &fastaPaths, 
  const std::string &blastDB 
  )
{
  // Create the blast data class object 
  BlastData blastData(fastaPaths, blastDB);
  
  // Return the dataframe 
  return blastData.createDataFrame();
}

// -----------------------------------------------------------------------------
