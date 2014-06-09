#ifndef __INJECTOR_HPP
#define __INJECTOR_HPP

#include <iostream>
#include "arr.h"
#include <string>
#include <regex>
#include <stdexcept>
#include <cfloat>

#include <boost/filesystem.hpp>
#include <boost/serialization/string.hpp>
#include <boost/mpi.hpp>

#include "fitsio.h"


using namespace std;



class fitsinfo {

public:

  fitsinfo() {}
  fitsinfo( string path, double start_time, double stop_time, long nrow, int hdu );

  string path;
  double start_time;
  double stop_time;
  long nrow;
  int hdu;

private:

  friend class boost::serialization::access;
  template<class Archive> void serialize(Archive & ar, const unsigned int version);

};



class injector {

public :
  
  injector ( string rootdir, string file_pattern, string det_pattern, double scale, bool add, int ntask, int id, int info=0 );

  int inject( arr < double > pntarr );
  
private :

  int find_file( double start_time );
  
  int seek_position( double start_time );
  
  long write_data( arr<double> pntarr, long offset );
  
  int get_startstop( string fname, double &start, double &stop, long &nrow, int &hdu );
  
  void test_fits_status( int status );

  void set_tscale( double t );
  
  const int ncol_ = 5; // Number of columns in the pntarr inputs
  const double tol_ = 1e-3;
  
  string rootdir_;
  string file_pattern_;
  string det_pattern_;
  bool add_;
  string fname_;
  int id_;
  int info_;
  double tscale_;
  double scale_;
  bool do_scale_;
  vector< fitsinfo > files;

  // This values pertain to the file presently open

  fitsinfo info;
  vector<double> time;
  vector<double> tod;
  long pos;

};


#endif
