#ifndef __CONVIQT_HPP__

// Include the autotools-provided configuration macros
#ifdef HAVE_CONFIG_H
#include"config.h"
#endif

#include <alm.h>

// The declarations go here

#define __CONVIQT_HPP__

class beam : public Alm<double> {
public :
  int read_beam( std::string fname );
private :
  ;
};

class sky : public Alm<double> {
public :
  int read_sky( std::string fname );
private :
  ;
};

class pointing {
public :
  ;
private :
  std::vector<double> arr;
};

class convolver {  
public :
  convolver( sky s, beam b );
  int convolve( pointing pnt );
  int set_sky( sky s );
  int set_beam( beam b );
private :
  beam b;
  sky s;
};

#endif
