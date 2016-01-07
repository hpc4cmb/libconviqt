#include "conviqt.hpp"

namespace conviqt {

void hpsort_DL( levels::arr<double> &ra, levels::arr<long> &brr ) {

  // Sort ra and brr based on ra.
  // This version requires lots of work space but is fast.
  
  if ( ra.size() != brr.size() ) {
    throw std::runtime_error( "Incompatible dimensions" );
  }

  size_t i, n=ra.size();

  // Generate a vector of pairs

  std::vector< std::pair<double, size_t> > pairs(n);
  for (i=0; i<n; ++i) pairs[i] = std::make_pair(ra[i], i);

  // Sort the pairs

  std::sort( pairs.begin(), pairs.end(), [](std::pair<double,size_t> p1, std::pair<double,size_t> p2) -> bool {return p1.first < p2.first;} );

  // Pull the sorted ar-values

  for ( i=0; i<n; ++i ) {
    ra[i] = pairs[i].first;
  }

  // Sort brr

  std::vector<long> temp(n);

  for ( i=0; i<n; ++i ) {
    size_t j = pairs[i].second;
    temp[i] = brr[j];
  }

  memcpy( &(brr[0]), &(temp[0]), sizeof(long)*n );
  
}
  
  
void hpsort_arr( levels::arr<double> &ra, int sortcol ) {
  
  if ( sortcol < 0 || sortcol > 4 ) {
    throw std::runtime_error( "Sort column must be between 0 and 4." );
  }
  
  // Sort 5-column array, ra based on the column sortcol
  
  size_t i, n=ra.size()/5;
  
  // Generate a vector of pairs

  std::vector< std::pair<double, size_t> > pairs(n);
  for (i=0; i<n; ++i) pairs[i] = std::make_pair(ra[5*i+sortcol], i);

  // Sort the pairs

  std::sort( pairs.begin(), pairs.end(), [](std::pair<double,size_t> p1, std::pair<double,size_t> p2) -> bool {return p1.first < p2.first;} );

  // Order each column based on the sorted indices

  std::vector< double > temp(n);
  
  for ( int col=0; col<5; ++col ) {
    if ( col == sortcol ) {
      // These elements are already sorted
      for ( i=0; i<n; ++i ) {
	ra[5*i+col] = pairs[i].first;
      }
    } else {
      // Use temporary workspace to sort this column
      for ( i=0; i<n; ++i ) {
	size_t j = pairs[i].second;
	temp[i] = ra[5*j+col];
      }
      for ( i=0; i<n; ++i ) {
	ra[5*i+col] = temp[i];
      }
    }
  }

}

  
void hpsort_arrTheta( levels::arr<double> &ra ) {

  hpsort_arr( ra, 1 );
  
}

void hpsort_arrTOD( levels::arr<double> &ra ) {

  hpsort_arr( ra, 3 );
  
}

  
void hpsort_arrTime( levels::arr<double> &ra ) {

  hpsort_arr( ra, 4 );

}

void hpsort_DDcm( levels::arr<double> &ra, levels::arr<double> &brr ) {

  if ( ra.size() != brr.size() ) {
    throw std::runtime_error( "Incompatible dimensions" );
  }
  
  // Sort ra and brr based on ra.
  // This version requires lots of work space but is fast.
  
  size_t i, n=ra.size();

  // Generate a vector of pairs

  std::vector< std::pair<double, size_t> > pairs(n);
  for (i=0; i<n; ++i) pairs[i] = std::make_pair(ra[i], i);

  // Sort the pairs

  std::sort( pairs.begin(), pairs.end(), [](std::pair<double,size_t> p1, std::pair<double,size_t> p2) -> bool {return p1.first < p2.first;} );

  // Sort brr using ra as workspace

  for ( i=0; i<n; ++i ) {
    size_t j = pairs[i].second;
    ra[i] = brr[j];
  }

  memcpy( &(brr[0]), &(ra[0]), sizeof(double)*n );
  
  // Pull the sorted ar-values

  for ( i=0; i<n; ++i ) {
    ra[i] = pairs[i].first;
  }
  
}

} // namespace conviqt
