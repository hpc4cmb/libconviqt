#include "injector.hpp"


float ntohf( const float inFloat ) {
  float retVal;
  char *floatToConvert = ( char* ) & inFloat;
  char *returnFloat = ( char* ) & retVal;
  
  // swap the bytes into a temporary buffer
  for (int i=0; i<4; ++i) returnFloat[i] = floatToConvert[3-i];

  return retVal;
}


double ntohd( const double inDouble ) {
  double retVal;
  char *doubleToConvert = ( char* ) & inDouble;
  char *returnDouble = ( char* ) & retVal;
  
  // swap the bytes into a temporary buffer
  for (int i=0; i<8; ++i) returnDouble[i] = doubleToConvert[7-i];

  return retVal;
}


injector::injector ( string rootdir, string file_pattern, string det_pattern, double scale, bool add, int ntask, int id, int info ) {
  rootdir_ = rootdir;
  file_pattern_ = file_pattern;
  det_pattern_ = det_pattern;
  if ( fabs(scale-1) > FLT_EPSILON ) {
    do_scale_ = true;
    scale_ = scale;
  } else {
    do_scale_ = false;
    scale_ = 1;
  }
  add_ = add;
  id_ = id;
  info_ = info;
  
  // Construct a list of matching files in the directory

  regex ftest( file_pattern_ );  

  cout.precision(16);

  vector< string > all_files;

  if ( id_ == 0 ) {

    if (info_ > 0) cout << id_ << " : Searching for files recursively " << endl;

    boost::filesystem::recursive_directory_iterator it = boost::filesystem::recursive_directory_iterator(rootdir_);
    boost::filesystem::recursive_directory_iterator end; // The default iterator can be used for end checking

    int i = -1;
    for ( ; it != end; ++it ) {
      
      //if ( !boost::filesystem::is_directory(*it) ) { // This test throws an exception on Edison

      string fname = it->path().string();
    
      if ( regex_match( fname, ftest ) ) all_files.push_back( fname );

      //}
      
    }

    if (info_ > 0) cout << id_ << " : Found " << all_files.size() << " files " << endl;
  }

  boost::mpi::communicator world;

  boost::mpi::broadcast( world, all_files, 0 );
  
  // Each process examines disjoint files and the results are collected in the end

  vector< fitsinfo > myfiles;

  for ( int i=0; i<all_files.size(); ++i ) {

    if ( i % ntask != id_ ) continue;

    double start, stop;
    long nrow;
    int hdu;
    bool double_precision;
    string fname = all_files[i];
    if ( !get_startstop( fname, start, stop, nrow, hdu, double_precision ) ) {
      if ( info_ > 0 ) cout << id_ << " : Found file : " << fname << " : " << start << " - " << stop << endl;
      myfiles.push_back( fitsinfo(fname, start, stop, nrow, hdu, double_precision) );
    }

  }

  // Broadcast the lists

  files.clear();
  vector< fitsinfo > recvfiles;  
  for ( int idsend=0; idsend<ntask; ++idsend ) {
    if ( id == idsend ) recvfiles = myfiles;
    boost::mpi::broadcast( world, recvfiles, idsend );
    files.insert( files.end(), recvfiles.begin(), recvfiles.end() );
  }

  if ( id_ == 0 ) {
    if ( files.size() == 0 )
      throw runtime_error( "No files match the search criteria: " + file_pattern + " and " + det_pattern );
    else
      cout << "Found " << files.size() << " exchange files matching " << file_pattern_ << endl;
  }

}


void injector::set_tscale( double t ) {
  // Uses the example time stamp to determine how to scale the time to OBT seconds.
  if ( t > 1e9 && t < 2e9 ) {
    tscale_ = 1.0; // Already right units
  } else if ( t > 1e13 && t < 1e15 ) {
    tscale_ = pow(2., -16.); // OBT ticks
  } else if ( t > 1e15 && t < 2e15 ) {
    tscale_ = 1e-6; // microseconds
  } else if ( t > 1e18 && t < 2e18 ) {
    tscale_ = 1e-9; // nanoseconds
  } else {
    throw runtime_error("Do not know how to scale the time stamps to seconds");
  }
}


void injector::test_fits_status( int status ) {
  if (status) {
    fits_report_error(stdout, status);
    throw runtime_error( "ERROR accessing fits file" );
  }
}


int injector::get_startstop( string fname, double &start, double &stop, long &nrow, int &hdu, bool &double_precision ) {

  // Examine the Planck Exchange file, fname, and retrieve start and stop times,
  // number or rows and which hdu matches det_pattern_.

  fitsfile *ffile;
  int status=0, hdutype, colnum, nhdu;
  char value[FLEN_VALUE], comment[FLEN_COMMENT];
  
  fits_open_file( &ffile, fname.c_str(), READONLY, &status );
  test_fits_status( status );

  // See if the file contains an HDU matching the detector name

  fits_get_num_hdus( ffile, &nhdu, &status );
  test_fits_status( status );

  regex dettest( det_pattern_ );

  hdu = -1;
  for ( int ihdu = 3; ihdu <= nhdu; ihdu++ ) {
    fits_movabs_hdu( ffile, ihdu, &hdutype, &status );
    test_fits_status( status );
    fits_read_key( ffile, TSTRING, "EXTNAME", value, comment, &status );
    test_fits_status( status );

    if ( regex_match( value, dettest ) ) {
      int typecode;
      long repeat, width;
      fits_get_coltype( ffile, 1, &typecode, &repeat, &width, &status );
      if ( typecode == TFLOAT )
        double_precision = false;
      else if ( typecode == TDOUBLE )
        double_precision = true;
      else
        throw runtime_error( string("Unsupported fits column type: ") + to_string(typecode) );
      test_fits_status( status );
      hdu = ihdu;
      break;
    }
  }

  if ( hdu < 1 ) return -1; // No match
  
  // Match found, get start and stop time

  fits_movabs_hdu( ffile, 2, &hdutype, &status );
  test_fits_status( status );

  fits_get_num_rows( ffile, &nrow, &status );
  test_fits_status( status );

  char key[8] = "OBT*";
  fits_get_colnum( ffile, CASEINSEN, key, &colnum, &status );
  test_fits_status(status);

  vector<long> obt(nrow);
  fits_read_col( ffile, TLONG, colnum, (long)1, (long)1, nrow, NULL, obt.data(), NULL, &status );
  test_fits_status(status);

  set_tscale(obt[0]);

  start = ((double)obt[0]) * tscale_;
  stop = ((double)obt[nrow-1]) * tscale_;

  //cout << " fits details : " << hdutype << " " <<  nrow << " " << colnum << " " << obt[0] << " " << obt[nrow-1] << " "<< start << " " << stop << tscale_ << endl;

  fits_close_file( ffile, &status );
  test_fits_status( status );

  return 0;
}


int injector::inject( arr < double > pntarr ) {

  // Loop over the TOD contained in the 4th column of pntarr
  // and the time stamps in the 5th column and
  // inject as many samples as possible from pntarr+offset to the
  // fits table starting at pos.

  long nwrite = pntarr.size() / ncol_;

  if ( nwrite == 0 ) return 0;

  long offset=0;

  double tstart;

  while ( nwrite > 0 ) {

    tstart = pntarr[offset + 4];

    find_file( tstart );

    seek_position( tstart );

    if ( info_ > 1 ) cout << id_ << " : Write position is in " << info.path << " hdu " << info.hdu << " row " << pos << endl;

    long n = write_data( pntarr, offset );
    
    offset += n * ncol_;
    nwrite -= n;

    if ( info_ > 0 ) cout << "Wrote " << n << " elements to " << info.path << " starting at " << pos << endl;

  }
  
  return 0;
}


int injector::find_file( double start_time ) {

  // Browses through the list of files to find the one with with matching time stamps
  // Loads the on-board time from the file so that it is efficient to seek inject positions
  // later. The file itself is only open during the read to allow concurrent access

  for ( vector<fitsinfo>::iterator it=files.begin(); it!=files.end(); ++it ) {
    if ( it->start_time < start_time+tol_ && start_time-tol_ < it->stop_time ) {

      if ( it->path == info.path ) return 0; // We already have this file open

      info = *it;

      time.resize( info.nrow );
      
      int status=0, hdutype;
      fitsfile *ffile;

      fits_open_file( &ffile, info.path.c_str(), READONLY, &status );
      test_fits_status( status );

      fits_movabs_hdu( ffile, 2, &hdutype, &status );
      test_fits_status( status );

      vector<long> obt(info.nrow);
      fits_read_col( ffile, TLONG, 1, (long)1, (long)1, info.nrow, NULL, obt.data(), NULL, &status );
      test_fits_status(status);

      fits_close_file( ffile, &status );
      test_fits_status( status );

      set_tscale(obt[0]);

      for ( long row=0; row<info.nrow; ++row ) time[row] = obt[row] * tscale_;

      return 0;
    }
  }

  // If we made it this far, no file matched the criteria
  
  throw runtime_error( "Unable to find an exchange file containing the start time: " + to_string(start_time) );

  return -1;
}


int injector::seek_position( double start_time ) {

  // Finds the first Exchange file time stamp that is >= start_time

  vector<double>::iterator low = lower_bound( time.begin(), time.end(), start_time-tol_ );

  if ( low == time.end() ) throw runtime_error( "Failed to find position in file = " + info.path + " start_time = " + to_string(start_time) );

  pos = low - time.begin();

  cout.precision(16);
  if ( info_ > 1 ) cout << "position: " << info.path << " " << pos << " / " << info.nrow << " " << start_time << " " << *low << endl;

  return 0;
}


long injector::write_data( arr<double> pntarr, long offset ) {

  // Write the maximal number of consecutive samples to the Planck Exchange File
  // The file is only open during the write operation to allow concurrent access

  long nwrite = 0;
  while ( abs( time[pos+nwrite] - pntarr[offset+nwrite*ncol_+4] ) < tol_ ) {
    ++nwrite;
    if ( pos+nwrite == time.size() ) break;
    if ( offset+nwrite*ncol_ == pntarr.size() ) break;
  }

  if ( nwrite == 0 ) {
    // Perhaps TOAST added missing samples that are not present in the TOD file. Try advancing the position
    while ( time[pos+nwrite] > pntarr[offset+nwrite*ncol_+4] - tol_ ) {
      if ( offset+nwrite*ncol_+4 > pntarr.size()-1 ) break;
      nwrite++;
    }
    if (nwrite > 0) return nwrite; // Virtual samples are not written
  }

  if ( nwrite == 0 ) {
    cerr.precision(16);
    cerr << id_ << " : No matching timestamps in " << info.path << " starting at " << pos << endl;
    for ( long i=0; i<10; ++i ) {
      if ( pos+i > time.size()-1 ) { cerr << "No more cached time stamps " << endl; break; }
      if ( offset+i*ncol_+4 > pntarr.size()-1 ) { cerr << "No more convolved data " << endl; break; }
      cerr << id_ << " : " << time[pos+i] << " " << pntarr[offset+i*ncol_+4] << " " << time[pos+i] - pntarr[offset+i*ncol_+4] << endl;
    }
    throw runtime_error(" Unable to find a single matching sample to write");
  }

  int status=0, hdutype;
  fitsfile *ffile;

  fits_open_file( &ffile, info.path.c_str(), READWRITE, &status );
  test_fits_status( status );

  fits_movabs_hdu( ffile, info.hdu, &hdutype, &status );
  test_fits_status( status );

  if ( false ) {

    // Old code to write using high level CFITSIO routines (VERY inefficient)

    vector<double> buffer(nwrite);
    
    if ( add_ ) {
      fits_read_col( ffile, TDOUBLE, 1, pos+1, (long)1, nwrite, NULL, buffer.data(), NULL, &status );
      test_fits_status(status);
      if ( do_scale_ ) 
        for ( long i=0; i<nwrite; ++i ) buffer[i] += pntarr[offset+i*ncol_+3] * scale_;
      else
        for ( long i=0; i<nwrite; ++i ) buffer[i] += pntarr[offset+i*ncol_+3];
    } else {
      if ( do_scale_ ) 
        for ( long i=0; i<nwrite; ++i ) buffer[i] = pntarr[offset+i*ncol_+3] * scale_;
      else
        for ( long i=0; i<nwrite; ++i ) buffer[i] = pntarr[offset+i*ncol_+3];
    }

    fits_write_col( ffile, TDOUBLE, 1, pos+1, (long)1, nwrite, buffer.data(), &status );
  } else {
    
    // New code: use raw I/O and bypass buffering and data translation

    long nbuffer;

    if ( info.double_precision )
      nbuffer = (long long) nwrite * (sizeof(double) + sizeof(char));
    else
      nbuffer = (long long) nwrite * (sizeof(float) + sizeof(char));

    vector<unsigned char> buffer( nbuffer );

    fits_read_tblbytes( ffile, pos+1, (long long)1, nbuffer, buffer.data(), &status );

    test_fits_status( status );

    int stride;
    if ( info.double_precision )
      stride = sizeof(double)/sizeof(unsigned char) + 1;
    else
      stride = sizeof(float)/sizeof(unsigned char) + 1;

    unsigned char * pbuffer = buffer.data();
    for ( long i=0; i<nwrite; ++i ) {
      if ( info.double_precision ) {
        double dtemp = pntarr[offset+i*ncol_+3];
        if ( do_scale_ ) dtemp *= scale_;
        if ( add_ ) dtemp += ntohd( *(double*)pbuffer );
        if ( htonl(47) != 47 ) dtemp = ntohd(dtemp); // Convert byte order
        memcpy(pbuffer, &dtemp, sizeof(double));
      } else {
        float ftemp = pntarr[offset+i*ncol_+3]; // Convert double to float
        if ( do_scale_ ) ftemp *= scale_;
        if ( add_ ) ftemp += ntohf( *(float*)pbuffer );
        if ( htonl(47) != 47 ) ftemp = ntohf( ftemp ); // Convert byte order
        memcpy(pbuffer, &ftemp, sizeof(float));
      }
      pbuffer += stride;
    }

    fits_write_tblbytes( ffile, pos+1, (long long)1, nbuffer, buffer.data(), &status);

    test_fits_status( status );
    
  }

  fits_close_file( ffile, &status );
  test_fits_status(status);

  return nwrite;
}


fitsinfo::fitsinfo( string path, double start_time, double stop_time, long nrow, int hdu, bool double_precision ) {
  this->path = path;
  this->start_time = start_time;
  this->stop_time = stop_time;
  this->nrow = nrow;
  this->hdu = hdu;
  this->double_precision = double_precision;
}


template<class Archive> void fitsinfo::serialize(Archive & ar, const unsigned int version) {
  ar & path;
  ar & start_time;
  ar & stop_time;
  ar & nrow;
  ar & hdu;
  ar & double_precision;
}
