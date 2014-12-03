//
//Module created by Gary Prezeau at JPL to calculate TOD exactly using using convolutions in spherical harmonic space performed by conviqt and LevelS macros without every writing anything to disk.
//
#include <complex>
#include <algorithm>
#include <iostream>
#include <iomanip>
#include <string>

//#include <toast_mpi.hpp>
#include "mpi.h"
//#include <toast.hpp>

#include "planck_rng.h"
#include "lsconstants.h"
#include "arr.h"
#include "cxxutils.h"
#include "focalplane_db.h"
#include "paramfile.h"
#include "pointing.h"
#include "vec3.h"
#include "rotmatrix.h"
#include "detector_pointing.h"
#include "sat_info.h"
#include "mpi_support.h"
#include "writers.h"
#include "iohandle_current.h"
#include "iohandle.h"
#include "repointing.h"
#include "levels_facilities.h"
#include "fftpack_support.h"
#include "wigner.h"
#include "alm.h"
#include "alm_dmcio.h"
#include "alm_powspec_tools.h"
#include "xcomplex.h"
#include "interpolator.h"
#include "trafos.h"
#include "walltimer.h"
#include <cmath>
#include <cstdlib>
#include <time.h>
#include <ctime>
#include "satquat.h"
#include "oofnoise.h"

using namespace std;

#define conv_acc 1.e-30

#define CMULT_VERBOSITY 0

/*
class toast_samples {

public :
  
  toast_samples ( ) { total_ = 0; nchunk_ = 0; }
  
  ~toast_samples ( ) { }
  
  long total_;
  long nchunk_;
  
  void operator () ( toast::chunk_p chnk, std::vector < toast::channel_p > chans ) {
    total_ += chnk->last() - chnk->first() + 1;
    nchunk_++;

    return;
  }

};


// This functor provides a callback method which is used to process every
// distributed data chunk and add data to the array.

class toast_chunk_process {

  public :
    
    toast_chunk_process ( bool galactic, string const & det_id, arr < double > * pntarr ) 
  {
      det_ = "/planck/" + det_id;
      gtype_ = toast::GEOM_THETAPHI;

      if ( galactic ) 
	{
	  coord_ = toast::PNTG_GALACTIC;
	} 
      else 
	{
	  coord_ = toast::PNTG_ECLIPTIC;
	}

      out_ = pntarr;
      good_ = 0;
  }
    
    ~toast_chunk_process ( ) 
  {

  }

    string det_;
    toast::geom_type gtype_;
    toast::pntg_coord coord_;
    arr < double > * out_;
    long good_;

    void operator () ( toast::chunk_p chnk, std::vector < toast::channel_p > chans ) {

      int64_t len = chnk->last() - chnk->first() + 1;

      // get streamset

      toast::node_p strmsetshr ( chnk->observation()->parent_ );
      toast::streamset_p strmset = strmsetshr->shared_ref < toast::streamset > ();
    
      // get telescope
  
      toast::node_p teleshr ( strmset->parent_ );
      toast::telescope_p tele = teleshr->shared_ref < toast::telescope > ();

      // Ideally, we would read the pointing for several channels at once,
      // but this would require the createTODPointingToast function to return
      // multiple channels in one call...
      //
      // Instead, we just search through the channels, select the one we want,
      // and get the pointing and data.

      std::vector < boost::numeric::ublas::vector < double > > chanpnt(1);
      chanpnt[0].resize( 3 * len );
      boost::numeric::ublas::vector < uint8_t > pntflags ( len );

      boost::numeric::ublas::vector < double > chandata ( len );
      boost::numeric::ublas::vector < uint8_t > flags ( len );
      boost::numeric::ublas::vector < double > timestamps ( len );


      std::vector < toast::channel_p > :: iterator itchan;
      long countc=0;
      for ( itchan = chans.begin(); itchan != chans.end(); ++itchan ) {
	//        if ( (*itchan)->name().find(det_) != string::npos ) {
        if ( (*itchan)->name() != det_ ) {
          // this is the channel we want!
          std::vector < toast::channel_p > chan(1);
          chan[0] = (*itchan);
          
          // read channel pointing
          
          tele->channel_geom_pointing ( chnk->observation(), chnk->first(), chnk->last(), chan, gtype_, coord_, chanpnt, pntflags );
          
          // read channel data and flags          
          
          (*itchan)->read ( chnk->observation(), "DEFAULT", chnk->first(), chnk->last(), chandata, flags );
          
          // read timestamps
          
          chnk->observation()->global_times ( chnk->first(), len, timestamps );
          
          int64_t chunkoff;
          long outoff;
          
          for ( int64_t i = 0; i < len; ++i ) {
	    if ( ( flags[i] | pntflags[i] ) == 0 ) {
              // no flags, data is good
              chunkoff = 3 * i;
	      //              outoff = 4 * good_;
              outoff = 5 * good_;
              
              (*out_)[ outoff ] = (chanpnt[0])( chunkoff + 1 );
              (*out_)[ outoff + 1 ] = (chanpnt[0])( chunkoff );
              (*out_)[ outoff + 2 ] = (chanpnt[0])( chunkoff + 2 );
              (*out_)[ outoff + 3 ] = chandata(i);
              (*out_)[ outoff + 4 ] = timestamps(i);
              ++good_;
            }
          }
          
        }
      }
      if ( CMULT_VERBOSITY > 1 ) cout << "good_ = " << good_ << "  len = " << len << endl; 
    }

};
*/
namespace {

double xpow (int expo, double val)
  { return (expo&1) ? -val : val; }

double wignerCalc00_halfpi (tsize n)
  {
  if (n&1) return 0.;
  double rec=1.;
  for (tsize k=1; k<n; k+=2)
    rec *= -(k/(k+1.));
  return rec;
  }

double wignerCalc00 (double theta, tsize n)
  {
  if (n==0) return 1.;
  double rec0=1., cth=cos(theta), rec1=cth;
  for (tsize k=2; k<=n; ++k)
    {
    double tmp = ((2*k-1.)*cth*rec1 - (k-1.)*rec0) / k;
    rec0 = rec1;
    rec1 = tmp;
    }
  return rec1;
  }

void wignerCalcHalfpi (tsize n, tsize mmax, arr2<double> &d)
{

  try {
    d.fast_alloc(2*n+1,2*mmax+1);
  } catch ( bad_alloc & e ) {
    cerr << " Out of memory allocating " << (2*n+1)*(2*mmax+1)*8./1024/1024 << "MB for wignerCalcHalfpi" << endl;
    throw;
  }
  
  if (n==0)
    {
    d[0][0]=1;
    return;
    }

  double d00;
  if ((n&1)==0)
    {
    double rec = d00 = d[n][mmax] = wignerCalc00_halfpi(n);
    if (1<=mmax) d[n-1][mmax-1] = rec;
    for (tsize k=2; k<=n; k+=2)
      {
      rec *= -sqrt(((n+k-1)*(n-k+2))/double((n-k+1)*(n+k)));
      d[n-k][mmax] = rec;
      if (k<=mmax) d[n][mmax-k] = rec;
      }
    for (tsize k=1; k<=n; k+=2)
      {
      d[n-k][mmax] = 0;
      if (k<=mmax) d[n][mmax-k] = 0;
      }
    }
  else
    {
    double rec = wignerCalc00_halfpi(n-1)*sqrt(n/(n+1.));
    d[n-1][mmax] = -rec;
    if (1<=mmax) d[n][mmax-1] = rec;
    for (tsize k=3; k<=n; k+=2)
      {
      rec *= -sqrt((n+k-1)*(n-k+2)/double((n-k+1)*(n+k)));
      d[n-k][mmax] = -rec;
      if (k<=mmax) d[n][mmax-k] = rec;
      }
    for (tsize k=0; k<=n; k+=2)
      {
      d[n-k][mmax] = 0;
      if (k<=mmax) d[n][mmax-k] = 0;
      }
    d00 = d[n-1][mmax];
    }

  arr<double> prod(n+1), sqtab(n);
  for (tsize i=0; i<n; ++i)
    sqtab[i] = 0.5*sqrt((2.*n-i)*(i+1.));

  for (tsize i=1; i<=mmax; ++i)
    {
    double di = double(i);
    double rec = prod[n-i] = -sqrt(.5*n)/di;

    for (tsize k=1; k<=n-i; ++k)
      prod[n-i-k] = rec = -sqtab[k] / (di + sqtab[k-1]*rec);

    tsize lowerLimit = i;
    if ((i==1) && ((n&1)==0)) lowerLimit = 2;
    double ratio=d00;
    for (tsize l=lowerLimit; l<=n; ++l)
      {
      ratio *= prod[l-i];
      d[n-l][mmax-i] = xpow(l-i,ratio);
      if (l<=mmax) d[n-i][mmax-l] = ratio;
      }
    if (i<mmax)
      d00 = -d[n-i][mmax-i-1];
    }
// now apply symmetry conditions to fill the rest of the array
  for (tsize j=1; j<=mmax; ++j)
    d[n][mmax+j] = xpow(n,d[n][mmax-j]);
  for (tsize i=1; i<=n; ++i)
    d[n+i][mmax] = xpow(i,d[n-i][mmax]);
  for (tsize i=1; i<=n; ++i)
    for (tsize j=1; j<=mmax; ++j)
      {
      d[n-i][mmax+j] = xpow(n+i,d[n-i][mmax-j]);
      d[n+i][mmax-j] = xpow(n+j,d[n-i][mmax-j]);
      d[n+i][mmax+j] = xpow(i+j,d[n-i][mmax-j]);
      }
  }

void wignerCalcGeneral (tsize n, tsize mmax, double theta, arr2<double> &d)
{
  
  try {
    d.fast_alloc(2*n+1,2*mmax+1);
  } catch ( bad_alloc & e ) {
    cerr << " Out of memory allocating " << (2*n+1)*(2*mmax+1)*8./1024/1024 << "MB for wignerCalcGeneral" << endl;
    throw;
  }

  double cth=cos(theta), sth=sin(theta);
  double xsth=1./sth;
  double fact1 = sqrt(2.*n)*cos(theta/2.)*sin(theta/2.);
  arr<double> prod1(n+1), prod2(n+1);

  arr<double> sqtab(n);
  for (tsize i=0; i<n; ++i)
    sqtab[i] = 0.5*sqrt((2.*n-i)*(i+1.));

  d[n][mmax] = wignerCalc00(theta,n);
  double d00 = d[n][mmax];
  double d00_2 = xpow(n,d00);
  double rec1 = prod1[n] = fact1/(n*cth);
  for (tsize k=1; k<n; ++k)
    prod1[n-k] = rec1 = -sqtab[k]/( -((n-k)*cth)*xsth + sqtab[k-1]*rec1 );
  double ratio1 = d00;
  for (tsize l=1; l<=n; ++l)
    {
    ratio1 *= prod1[l];
    d[n-l][mmax] = xpow(l,ratio1);
    if (l<=mmax)
      {
      d[n][mmax+l] = xpow(l,ratio1);
      d[n][mmax-l] = ratio1;
      }
    }
  if (n>0)
    {
    d00 = d[n][mmax+1];
    d00_2 = d[n][mmax-1];
    }
  for (tsize i=1; i<=mmax; ++i)
    {
    //I STILL NEED TO WORRY ABOUT THE PROD1 ARRAY HAVING VANISHING DENOMINATORS
    bool flag1=false, flag2=false;
    if (abs(d00) < 1e-14)
      {
      flag1=true;
      d[n-i][mmax-i] = -sqtab[n-(i-1)]/sqtab[n-i]*d[n-i+2][mmax-i];
      }
    if (abs(d00_2) < 1e-14)
      {
      flag2=true;
      d[n-i][mmax+i] = -sqtab[n-(i-1)]/sqtab[n-i]*d[n-i+2][mmax+i];
      }
    rec1 = prod1[n-i] = fact1/(n*cth-i);
    double rec2 = prod2[n-i] = -fact1/(n*cth+i);
    for (tsize k=1; k<=n-i; ++k)
      {
      if ((flag1) && (k==n-i)) continue;
      if (abs(prod1[n-i-k+1]) >= 1e14)
        prod1[n-i-k+1] = rec1 = 1e16;
      prod1[n-i-k] = rec1 = -sqtab[k]/( (i-(n-k)*cth)*xsth + sqtab[k-1]*rec1 );
      }
    for (tsize k=1; k<=n-i; ++k)
      {
      if ((flag2) && (k==n-i)) continue;
      if (abs(prod2[n-i-k+1]) >= 1e14)
        prod2[n-i-k+1] = rec2 = 1e16;
      prod2[n-i-k] = rec2 = -sqtab[k]/( (i+(n-k)*cth)*xsth + sqtab[k-1]*rec2 );
      }
    ratio1 = d00;
    double ratio2 = d00_2;
    if (flag1) ratio1 = d[n-i][mmax-i];
    if (flag2) ratio2 = d[n-i][mmax+i];
    for (tsize l=i; l<=n; ++l)
      {
      if ((flag1) && (l==i)) continue;
      ratio1 *= prod1[l-i];
      double p_li = xpow(l+i,1.);
      d[n-l][mmax-i] = p_li*ratio1;
      if (l<=mmax)
        d[n-i][mmax-l] = ratio1;
      }
    for (tsize l=i; l<=n; ++l)
      {
      if ((flag2) && (l==i)) continue;
      ratio2 *= -prod2[l-i];
      double p_li = xpow(l+i,1.);
      d[n-l][mmax+i] = p_li*ratio2;
      if (l<=mmax)
        d[n-i][mmax+l] = p_li*ratio2;
      }
    if (i<mmax)
      {
      d00 = -d[n-i][mmax-i-1];
      d00_2 = -d[n-i][mmax+i+1];
      }
    }
// now apply symmetry conditions to fill the rest of the array
  for (tsize i=1; i<=n; ++i)
    for (tsize j=0; j<=2*mmax; ++j)
      d[n+i][j] = xpow(j+i+mmax,d[n-i][2*mmax-j]);
  }

} // unnamed namespace

void wignerCalc (tsize n, tsize mmax, double theta, arr2<double> &d)
{
  try {
    d.fast_alloc (2*n+1,2*mmax+1);
  } catch ( bad_alloc & e ) {
    cerr << " Out of memory allocating " << (2*n+1)*(2*mmax+1)*8./1024/1024 << "MB for wignerCalc" << endl;
    throw;
  }

  tdiff mm=mmax;
  if (n == 0)
    {
      d[0][0]=1.;
      return;
    }
  if (abs(theta) >= twopi)
    {
      double thetasign = 1.;
      if (theta < 0) thetasign = -1.;
      long anglefac(0);
      double anglediff=2.*twopi;
      while ( abs(anglediff) > twopi)
	{
	  anglefac++;
	  anglediff = theta - thetasign*anglefac*twopi;
	}
      theta = theta - thetasign*anglefac*twopi;
    }
  if ((abs_approx(theta, 0., 1e-14)) || (abs_approx(abs(theta), twopi, conv_acc)))
    {
    d.fill(0.);
    for (tdiff k=-mm; k<=mm; ++k)
      d[n+k][mmax+k] = 1.;
    return;
    }
  if (abs_approx(theta, pi, 1e-14))
    {
    d.fill(0.);
    for (tdiff k=-mm; k<=mm; ++k)
      d[n+k][mmax+k] = xpow(n+k,1.);
    return;
    }
  if (abs_approx(theta, halfpi, 1e-14))
    {
    wignerCalcHalfpi (n,mmax,d);
    return;
    }
  wignerCalcGeneral (n,mmax,theta,d);
  }

void readinblm(int corenum, long newcorenum, long coreupper, Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &blmC, long lmax, long beammmax, bool pol, paramfile &par_file)
{
  string infile_beam = par_file.find<string>("infile_beam");
  long blmsize=(beammmax+1)*(beammmax+2) + 2*(beammmax+1)*(lmax-beammmax);
  if (corenum==newcorenum)
    { 
      pol ? read_Alm_from_dmc (infile_beam,blmT, blmG, blmC, lmax, beammmax)
	: read_Alm_from_dmc (infile_beam,blmT, lmax, beammmax);
      if ( CMULT_VERBOSITY > 1 ) cout << "Done reading beam alms in core = " << corenum << endl;
    }
  if(corenum==newcorenum) for (long ii=newcorenum+1; ii<coreupper; ii++) mpiMgr.sendRaw(&blmT(0,0).re, blmsize,ii);
  if (corenum != newcorenum) mpiMgr.recvRaw (&blmT(0,0).re, blmsize,newcorenum);
  if (pol)
    {
      if(corenum==newcorenum) for (long ii=newcorenum+1; ii<coreupper; ii++) mpiMgr.sendRaw(&blmG(0,0).re, blmsize,ii);
      if (corenum != newcorenum) mpiMgr.recvRaw (&blmG(0,0).re, blmsize,newcorenum);
      
      if(corenum==newcorenum) for (long ii=newcorenum+1; ii<coreupper; ii++) mpiMgr.sendRaw(&blmC(0,0).re, blmsize,ii);
      if (corenum != newcorenum) mpiMgr.recvRaw (&blmC(0,0).re, blmsize,newcorenum);
    }
}

void readinslm(int corenum, long newcorenum, long coreupper, Alm<xcomplex<float> > &slmInputT, Alm<xcomplex<float> > &slmInputG, Alm<xcomplex<float> > &slmInputC, long lmax, bool pol, paramfile &par_file)
{
  string infile_sky = par_file.find<string>("infile_sky");
  long slmsize=(lmax+1)*(lmax+2);
  if (corenum==newcorenum)
    { 
      pol ? read_Alm_from_dmc (infile_sky,slmInputT, slmInputG, slmInputC, lmax, lmax)
	: read_Alm_from_dmc (infile_sky,slmInputT, lmax, lmax);
      if ( CMULT_VERBOSITY > 1 ) cout << "Done reading sky alms in core = " << corenum << endl;
    }
  if(corenum==newcorenum) for (long ii=newcorenum+1; ii<coreupper; ii++) mpiMgr.sendRaw(&slmInputT(0,0).re, slmsize,ii);
  if (corenum != newcorenum) mpiMgr.recvRaw (&slmInputT(0,0).re, slmsize,newcorenum);
  if (pol)
    {
      if(corenum==newcorenum) for (long ii=newcorenum+1; ii<coreupper; ii++) mpiMgr.sendRaw(&slmInputG(0,0).re, slmsize,ii);
      if (corenum != newcorenum) mpiMgr.recvRaw (&slmInputG(0,0).re, slmsize,newcorenum);
  
      if(corenum==newcorenum) for (long ii=newcorenum+1; ii<coreupper; ii++) mpiMgr.sendRaw(&slmInputC(0,0).re, slmsize,ii);
      if (corenum != newcorenum) mpiMgr.recvRaw (&slmInputC(0,0).re, slmsize,newcorenum);
    }
  double fwhm = arcmin2rad*par_file.find<double>("fwhm_deconv",0.);
  pol ? smoothWithGauss (slmInputT, slmInputG, slmInputC, -fwhm)
      : smoothWithGauss (slmInputT, -fwhm);
}

void sift_down_DL(arr<double> &ra, arr<long> &brr, const int l, const int r) //FROM P-340 OF NR
{
  int j, jold;
  double a;
  long b;
  a=ra[l];
  b=brr[l];
  jold=l;
  j=2*l+1;
  while(j<=r)
    {
      if (j<r && ra[j] < ra[j+1]) j++;
      if (a>= ra[j]) break;
      ra[jold]=ra[j];
      brr[jold]=brr[j];
      jold=j;
      j=2*j+1;
    }
  ra[jold]=a;
  brr[jold]=b;
}

  
void hpsort_DL(arr<double> &ra, arr<long> &brr) //FROM P-340 OF NR
{
  int i;
  int n=ra.size();
  for(i=n/2-1; i>=0; i--)
    sift_down_DL(ra,brr,i,n-1);
  for (i=n-1; i>0; i--)
    {
      swap(ra[0], ra[i]);
      swap(brr[0], brr[i]);
      sift_down_DL(ra,brr,0,i-1);
    }
}

void sift_down_arrTime(arr<double> &ra, const int l, const int r) //FROM P-340 OF NR
{
  int j, jold;
  double a,c,d,f;
  double b;
  a=ra[5*l+4];
  b=ra[5*l+0];
  f=ra[5*l+1];
  c=ra[5*l+2];
  d=ra[5*l+3];
  jold=l;
  j=2*l+1;
  while(j<=r)
    {
      if (j<r && ra[5*j+4] < ra[5*(j+1)+4]) j++;
      if (a>= ra[5*j+4]) break;
      ra[5*jold+4]=ra[5*j+4];
      ra[5*jold+2]=ra[5*j+2];
      ra[5*jold+0]=ra[5*j+0];
      ra[5*jold+1]=ra[5*j+1];
      ra[5*jold+3]=ra[5*j+3];
      jold=j;
      j=2*j+1;
    }
  ra[5*jold+4]=a;
  ra[5*jold+2]=c;
  ra[5*jold+0]=b;
  ra[5*jold+1]=f;
  ra[5*jold+3]=d;
}

  
void hpsort_arrTime(arr<double> &ra) //FROM P-340 OF NR
{
  int i;
  int n=ra.size()/5;
  for(i=n/2-1; i>=0; i--)
    sift_down_arrTime(ra,i,n-1);
  for (i=n-1; i>0; i--)
    {
      swap(ra[4], ra[5*i+4]);
      swap(ra[0], ra[5*i+0]);
      swap(ra[1], ra[5*i+1]);
      swap(ra[2], ra[5*i+2]);
      swap(ra[3], ra[5*i+3]);
      sift_down_arrTime(ra,0,i-1);
    }
}

void sift_down_arrTheta(arr<double> &ra, const int l, const int r) //FROM P-340 OF NR
{
  int j, jold;
  double a,c,d,f;
  double b;
  a=ra[5*l+1];
  b=ra[5*l+0];
  f=ra[5*l+4];
  c=ra[5*l+2];
  d=ra[5*l+3];
  jold=l;
  j=2*l+1;
  while(j<=r)
    {
      if (j<r && ra[5*j+1] < ra[5*(j+1)+1]) j++;
      if (a>= ra[5*j+1]) break;
      ra[5*jold+1]=ra[5*j+1];
      ra[5*jold+2]=ra[5*j+2];
      ra[5*jold+0]=ra[5*j+0];
      ra[5*jold+4]=ra[5*j+4];
      ra[5*jold+3]=ra[5*j+3];
      jold=j;
      j=2*j+1;
    }
  ra[5*jold+1]=a;
  ra[5*jold+2]=c;
  ra[5*jold+0]=b;
  ra[5*jold+4]=f;
  ra[5*jold+3]=d;
}

  
void hpsort_arrTheta(arr<double> &ra) //FROM P-340 OF NR
{
  int i;
  int n=ra.size()/5;
  for(i=n/2-1; i>=0; i--)
    sift_down_arrTheta(ra,i,n-1);
  for (i=n-1; i>0; i--)
    {
      swap(ra[1], ra[5*i+1]);
      swap(ra[0], ra[5*i+0]);
      swap(ra[4], ra[5*i+4]);
      swap(ra[2], ra[5*i+2]);
      swap(ra[3], ra[5*i+3]);
      sift_down_arrTheta(ra,0,i-1);
    }
}

void sift_down_arrTOD(arr<double> &ra, const int l, const int r) //FROM P-340 OF NR
{
  int j, jold;
  double a,c,d,f;
  double b;
  a=ra[5*l+3];
  b=ra[5*l+0];
  f=ra[5*l+4];
  c=ra[5*l+2];
  d=ra[5*l+1];
  jold=l;
  j=2*l+1;
  while(j<=r)
    {
      if (j<r && ra[5*j+3] < ra[5*(j+1)+3]) j++;
      if (a>= ra[5*j+3]) break;
      ra[5*jold+3]=ra[5*j+3];
      ra[5*jold+2]=ra[5*j+2];
      ra[5*jold+0]=ra[5*j+0];
      ra[5*jold+4]=ra[5*j+4];
      ra[5*jold+1]=ra[5*j+1];
      jold=j;
      j=2*j+1;
    }
  ra[5*jold+3]=a;
  ra[5*jold+2]=c;
  ra[5*jold+0]=b;
  ra[5*jold+4]=f;
  ra[5*jold+1]=d;
}


void hpsort_arrTOD(arr<double> &ra) //FROM P-340 OF NR
{
  int i;
  int n=ra.size()/5;
  for(i=n/2-1; i>=0; i--)
    sift_down_arrTOD(ra,i,n-1);
  for (i=n-1; i>0; i--)
    {
      swap(ra[3], ra[5*i+3]);
      swap(ra[0], ra[5*i+0]);
      swap(ra[4], ra[5*i+4]);
      swap(ra[2], ra[5*i+2]);
      swap(ra[1], ra[5*i+1]);
      sift_down_arrTOD(ra,0,i-1);
    }
}

void sift_down_DDcm(arr<double> &ra, arr<double> &brr, const int l, const int r) //FROM P-340 OF NR                                                                   
{
  int j, jold;
  double a;
  double b;
  a=ra[l];
  b=brr[l];
  jold=l;
  j=2*l+1;
  while(j<=r)
    {
      if (j<r && ra[j] < ra[j+1]) j++;
      if (a>= ra[j]) break;
      ra[jold]=ra[j];
      brr[jold]=brr[j];
      jold=j;
      j=2*j+1;
    }
  ra[jold]=a;
  brr[jold]=b;
}


void hpsort_DDcm(arr<double> &ra, arr<double> &brr) //FROM P-340 OF NR                                                                                                
{
  int i;
  int n=ra.size();
  for(i=n/2-1; i>=0; i--)
    sift_down_DDcm(ra,brr,i,n-1);
  for (i=n-1; i>0; i--)
    {
      swap(ra[0], ra[i]);
      swap(brr[0], brr[i]);
      sift_down_DDcm(ra,brr,0,i-1);
    }
}

void ratiobetacalcgreatercm(long &ratiobeta, double theta, int cores, arr<double> &corethetaarr)
{
  for (long ii=0; ii<cores; ii++)
    {
      if (ratiobeta+1>=cores) break;
      if (theta<corethetaarr[ratiobeta+1]) break;
      ratiobeta++;
    }
}

void ratiobetacalcsmallercm(long &ratiobeta, double theta, int cores, arr<double> &corethetaarr)
{
  for (long ii=0; ii<100000; ii++)
    {
      if (ratiobeta<0) break;
      if (theta>corethetaarr[ratiobeta]) break;
      ratiobeta--;
    }
}
void thetaDeltaThetacm(double thetaini, int cores, int corenum, double &theta, double &deltatheta)
{
  if (corenum==0)
    {
      theta=0.;
      deltatheta=thetaini;
    }
  else if (corenum==cores-1)
    {
      deltatheta=1./cores;
      theta=thetaini-deltatheta;
    }
  else
    {
      double a=-cos(thetaini)/2.;
      double b=-sin(thetaini)+thetaini*cos(thetaini);
      double c=thetaini*sin(thetaini) - cos(thetaini)*thetaini*thetaini/2.-1./cores;
      theta=(-b-sqrt(b*b-4.*a*c))/(2.*a);
      deltatheta=thetaini-theta;
    }
}

void todRedistribution5cm (arr<double> pntarr, arr<int> inBetaSeg, arr<int> outBetaSeg, arr<int> &inBetaSegAcc, arr<int> &outBetaSegAcc, long &outBetaSegSize, int cores, arr<double> &pntarr2, long totsize, arr<int> &inOffset, arr<int> &outOffset, double ratiodeltas, arr<double> &corethetaarr)
{
  if ( CMULT_VERBOSITY > 1 ) cout << "Entered todRedistribution5cm" << endl;
  inBetaSegAcc.alloc(cores);
  inBetaSegAcc.fill(0);
  outBetaSegAcc.alloc(cores);
  outBetaSegAcc.fill(0);
  for (long ii=0; ii<cores; ii++)
    for (long jj=0; jj<=ii; jj++)
      {
      inBetaSegAcc[ii]+=inBetaSeg[jj];
      outBetaSegAcc[ii]+=outBetaSeg[jj];
      }
  outBetaSegSize=outBetaSegAcc[cores-1];
  inOffset.alloc(cores);
  inOffset.fill(0);
  outOffset.alloc(cores);
  outOffset.fill(0);
  for (long ii=1; ii<cores; ii++) {outOffset[ii]=outBetaSegAcc[ii-1]; inOffset[ii]=inBetaSegAcc[ii-1];}
  if (totsize>0)
    {
      try {
        pntarr2.alloc(5*totsize);
      } catch ( bad_alloc & e ) {
        cerr << " todRedistribution5cm : Out of memory allocating " << (5*totsize)*8./1024/1024 << "MB for pntarr2" << endl;
        throw;
      }
      pntarr2.fill(0);
      arr<int> countbeta(cores,0);
      for (long ii=0; ii<totsize; ii++)
	{
	  long nbeta, ratiobeta;
	  double theta=pntarr[5*ii+1];
	  if (theta >= halfpi) theta=pi-theta;
	  ratiobeta=theta/ratiodeltas;
	  if (theta>=corethetaarr[ratiobeta])
	    ratiobetacalcgreatercm(ratiobeta, theta, cores, corethetaarr);
	  else
	    ratiobetacalcsmallercm(ratiobeta, theta, cores, corethetaarr);
	  pntarr2[inOffset[ratiobeta] + countbeta[ratiobeta]]=pntarr[5*ii];
	  pntarr2[inOffset[ratiobeta] + countbeta[ratiobeta]+1]=pntarr[5*ii+1];
	  pntarr2[inOffset[ratiobeta] + countbeta[ratiobeta]+2]=pntarr[5*ii+2];
	  pntarr2[inOffset[ratiobeta] + countbeta[ratiobeta]+3]=pntarr[5*ii+3];
	  pntarr2[inOffset[ratiobeta] + countbeta[ratiobeta]+4]=pntarr[5*ii+4];
	  countbeta[ratiobeta]+=5;
	}
    }
  if ( CMULT_VERBOSITY > 1 ) cout << "Leaving todRedistribution5cm" << endl;
}


void weight_ncm (double x, arr<double> &wgt, long order)
  {
    int npoints=order+1;
    arr<double> base_wgt(npoints);
    base_wgt.fill(1);
    for (int m=0; m<npoints; ++m)
      for (int n=0; n<npoints; ++n)
	if (m!=n) base_wgt[m] *= m-n;
    for (int m=0; m<npoints; ++m) base_wgt[m] = 1./base_wgt[m];
    for (int m=0; m<npoints; ++m)
      wgt[m] = base_wgt[m];
    double mul1=x, mul2=x-npoints+1;
    for (int m=1; m<npoints; ++m)
      {
	wgt[m]*=mul1;
	wgt[npoints-m-1]*=mul2;
	mul1*=x-m;
	mul2*=x-npoints+m+1;
      }
  }

void todAnnulus_v3(arr3<xcomplex<double> > &tod1, arr3<xcomplex<double> > &Cmm, long lmax, long beammmax, long lmaxOut, arr<double> &cs, arr<double> &sn, arr<double> &cs0, arr<double> &sn0, long NThetaIndex1)
{
  long binElements=2*lmaxOut+1;
  for (long msky = -lmax; msky <= lmax; msky++)
    for (long beamIndex=0; beamIndex<=beammmax; beamIndex++)
      for (long lat=0; lat<NThetaIndex1; lat++)
	{
	  double dPhi=-halfpi;
	  double tmpR=Cmm(msky+lmaxOut, beamIndex, lat).real(), tmpI=Cmm(msky+lmaxOut, beamIndex, lat).imag();
	  Cmm(msky+lmaxOut, beamIndex, lat).real() =(cs[msky+lmax]*tmpR + sn[msky+lmax]*tmpI);
	  Cmm(msky+lmaxOut, beamIndex, lat).imag() =(cs[msky+lmax]*tmpI - sn[msky+lmax]*tmpR);
	}
  for (long beamIndex=0; beamIndex<=beammmax; beamIndex++)
    for (long lat=0; lat<NThetaIndex1; lat++)
      {
	arr<xcomplex<double> > Cmsky(binElements, 0);
	for (long msky=-lmax; msky<=lmax; msky++) Cmsky[msky+lmaxOut]=Cmm(msky+lmaxOut, beamIndex, lat);
	cfft p1(binElements);
	p1.backward(Cmsky);
	for (long msky=0; msky<binElements; msky++) tod1(msky, beamIndex, lat)=Cmsky[msky];
      }

  for (long msky = -lmaxOut; msky <= lmaxOut; msky++)
    for (long beamIndex=0; beamIndex<=beammmax; beamIndex++)
      for (long lat=0; lat<NThetaIndex1; lat++)
	{
	  double tmpR=tod1(msky+lmaxOut, beamIndex, lat).real(), tmpI=tod1(msky+lmaxOut, beamIndex, lat).imag();
	  tod1(msky+lmaxOut, beamIndex, lat).real() =(cs0[lmaxOut+msky]*tmpR + sn0[lmaxOut+msky]*tmpI);
	  tod1(msky+lmaxOut, beamIndex, lat).imag() =(cs0[lmaxOut+msky]*tmpI - sn0[lmaxOut+msky]*tmpR);
	  //Rotate in psi space
	  tmpR=tod1(msky+lmaxOut, beamIndex, lat).real(); tmpI=tod1(msky+lmaxOut, beamIndex, lat).imag();
	  tod1(msky+lmaxOut, beamIndex, lat).real() =(cs[beamIndex+lmax]*tmpR + sn[beamIndex+lmax]*tmpI);
	  tod1(msky+lmaxOut, beamIndex, lat).imag() =(cs[beamIndex+lmax]*tmpI - sn[beamIndex+lmax]*tmpR);
	}
  //Note that now, -lmax <= msky <= lmax corresponds to the angle 0 <= phi <= 2.*pi/2./(2.*lmax+1.) and -pi <= phi <= -2.pi/(2.*lmax+1.)
  //Finished with convolution and FFT over msky.
}

void todAnnulus(arr2<xcomplex<double> > &tod1, arr2<xcomplex<double> > &Cmm, long lmax, long beammmax, long lmaxOut, arr<double> &cs, arr<double> &sn, arr<double> &cs0, arr<double> &sn0)
{
  long binElements=2*lmaxOut+1;
  for (long beamIndex=0; beamIndex<=beammmax; beamIndex++)
    for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	double tmpR=Cmm[msky+lmaxOut][beamIndex].real(), tmpI=Cmm[msky+lmaxOut][beamIndex].imag();
	Cmm[msky+lmaxOut][beamIndex].real() =(cs[msky+lmax]*tmpR + sn[msky+lmax]*tmpI);
	Cmm[msky+lmaxOut][beamIndex].imag() =(cs[msky+lmax]*tmpI - sn[msky+lmax]*tmpR);
      }
  for (long beamIndex=0; beamIndex<=beammmax; beamIndex++)
    {
      arr<xcomplex<double> > Cmsky(binElements);
      for (long msky=0; msky<binElements; msky++)
	{
	  Cmsky[msky]=Cmm[msky][beamIndex];
	}
      cfft p1(binElements);
      p1.backward(Cmsky);
      for (long msky=0; msky<binElements; msky++)
	{
	  tod1[msky][beamIndex]=Cmsky[msky];
	}
      for (long msky = -lmaxOut; msky <= lmaxOut; msky++)
	{
	  double tmpR=tod1[msky+lmaxOut][beamIndex].real(), tmpI=tod1[msky+lmaxOut][beamIndex].imag();
	  tod1[msky+lmaxOut][beamIndex].real() =(cs0[lmaxOut+msky]*tmpR + sn0[lmaxOut+msky]*tmpI);
	  tod1[msky+lmaxOut][beamIndex].imag() =(cs0[lmaxOut+msky]*tmpI - sn0[lmaxOut+msky]*tmpR);
	  //Rotate in psi space
	  tmpR=tod1[msky+lmaxOut][beamIndex].real(); tmpI=tod1[msky+lmaxOut][beamIndex].imag();
	  tod1[msky+lmaxOut][beamIndex].real() =(cs[beamIndex+lmax]*tmpR + sn[beamIndex+lmax]*tmpI);
	  tod1[msky+lmaxOut][beamIndex].imag() =(cs[beamIndex+lmax]*tmpI - sn[beamIndex+lmax]*tmpR);
	}
    }
  //Note that now, -lmax <= msky <= lmax corresponds to the angle 0 <= phi <= 2.*pi/2./(2.*lmax+1.) and -pi <= phi <= -2.pi/(2.*lmax+1.)
  //Finished with convolution and FFT over msky.
}

void conviqt_hemiscm_pol_fast(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long beammmax, long lmax, long lmaxOut, arr2<xcomplex<double> > &tod1, arr2<xcomplex<double> > &tod2, double beta, long &countdlm, arr<long> &effM, arr<double> &dp1, int corenum)
{
  long dp1size=dp1.size();
  long binElements=2*lmaxOut+1;
  long psiElements=beammmax+1;
  arr2<xcomplex<double> > Cmm(binElements, psiElements,0.);
  arr2<xcomplex<double> > Cmm2(binElements, psiElements,0.);

  for (long ii=0; ii<=lmax; ii++)
    {
      double dii=double(ii);
      double dsignl=pow(-1.,dii);

      long beamupper(beammmax);
      if (ii<beammmax) beamupper=ii;
      long mskyUpper=(ii<effM[ii]) ? ii : effM[ii];
      for(long beamIndex = 0; beamIndex <= beamupper; beamIndex++)
	{
	  double dbeamIndex = double(beamIndex);
	  double dsignb=pow(-1.,dbeamIndex), dlb=dsignb*dsignl;
	  for (long msky = 0; msky <= mskyUpper; msky++)
	    {
	      double dmsky=double(msky);
	      double dsign=pow(-1.,dmsky), dsb=dsign*dsignb;
	      double dsignt=dsb*dsignl;
	      //	      if (countdlm>=dp1size) {cout << "dp1size = " << dp1size << "  countdlm = " << countdlm << "  corenum = " << corenum << "   ii = " << ii << "  beamIndex = " << beamIndex << "  msky = " << msky << "   effM[lmax]= " << effM[lmax] << endl; planck_fail("countdlm != dp1size");}
	      double dMatrixElementmskypos = dp1[countdlm];//Note that msky in dlm is located to the left of mb and they have to be interchanged in the convolution
	      countdlm++;
	      //	      if (countdlm>=dp1size) {cout << "dp1size = " << dp1size << "  countdlm = " << countdlm << "  corenum = " << corenum << "   ii = " << ii << "  beamIndex = " << beamIndex << "  msky = " << msky << "   effM[lmax]= " << effM[lmax] << endl; planck_fail("countdlm != dp1size");}
	      double dMatrixElementmskyneg = dp1[countdlm];//dp1[ii-msky][beamIndex+beamupper];
	      countdlm++;
	      double dMatrixElementmskypos2 = dlb*dMatrixElementmskyneg;
	      double dMatrixElementmskyneg2 = dlb*dMatrixElementmskypos;
	      double prod1=slmT(ii,msky).re*blmT(ii,beamIndex).re + slmG(ii,msky).re*blmG(ii,beamIndex).re + slmC(ii,msky).re*blmC(ii,beamIndex).re;
	      double prod3=slmT(ii,msky).im*blmT(ii,beamIndex).re + slmG(ii,msky).im*blmG(ii,beamIndex).re + slmC(ii,msky).im*blmC(ii,beamIndex).re;
	      double prod2=slmT(ii,msky).im*blmT(ii,beamIndex).im + slmG(ii,msky).im*blmG(ii,beamIndex).im + slmC(ii,msky).im*blmC(ii,beamIndex).im;
	      double prod4=slmT(ii,msky).re*blmT(ii,beamIndex).im + slmG(ii,msky).re*blmG(ii,beamIndex).im + slmC(ii,msky).re*blmC(ii,beamIndex).im;
	      double tmp_1=prod1 + prod2;
	      double tmp_2=prod3 - prod4;
	      double tmp_3=prod1 - prod2;
	      double tmp_4=-prod3 - prod4;
	      double xtmp_1=tmp_1 * dMatrixElementmskypos;
	      double xtmp_2=tmp_2 * dMatrixElementmskypos;
	      double xtmp_3=dsign*tmp_3 * dMatrixElementmskyneg;
	      double xtmp_4=dsign*tmp_4 * dMatrixElementmskyneg;
	      double xtmp_5=tmp_1 * dMatrixElementmskypos2;
	      double xtmp_6=tmp_2 * dMatrixElementmskypos2;
	      double xtmp_7=dsign*tmp_3 * dMatrixElementmskyneg2;
	      double xtmp_8=dsign*tmp_4 * dMatrixElementmskyneg2;
	      Cmm[msky+lmaxOut][beamIndex].real() += xtmp_1;
	      Cmm[msky+lmaxOut][beamIndex].imag() += xtmp_2;
	      Cmm2[msky+lmaxOut][beamIndex].real() += xtmp_5;
	      Cmm2[msky+lmaxOut][beamIndex].imag() += xtmp_6;
	      if (msky!=0)
		{
		  Cmm[-msky+lmaxOut][beamIndex].real() += xtmp_3;
		  Cmm[-msky+lmaxOut][beamIndex].imag() += xtmp_4;
		  Cmm2[-msky+lmaxOut][beamIndex].real() += xtmp_7;
		  Cmm2[-msky+lmaxOut][beamIndex].imag() += xtmp_8;
		}
	    }
	}
    }
  arr<double> cs(2*lmax+1), sn(2*lmax+1), cs0(binElements), sn0(binElements);
  for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	cs[lmax+msky] =cos(msky*dPhi);
	sn[lmax+msky] =sin(msky*dPhi);
      }
  for (long msky=-lmaxOut; msky<=lmaxOut; msky++)
    {
      cs0[lmaxOut+msky]=cos(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
      sn0[lmaxOut+msky]=sin(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
    }
  todAnnulus(tod1, Cmm, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0);
  todAnnulus(tod2, Cmm2, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0);
}
/*
inline bool exists_test3 (const std::string& name) 
{
  struct stat buffer;   
  return (stat (name.c_str(), &buffer) == 0); 
}
*/
void conviqt_hemiscm_pol_fill(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long beammmax, long lmax, long lmaxOut, arr2<xcomplex<double> > &tod1, arr2<xcomplex<double> > &tod2, double beta, arr<long> &effM, long &countdlm, int corenum, arr<double> &dlmelements)
{
  long binElements=2*lmaxOut+1;
  long psiElements=beammmax+1;
  arr2<xcomplex<double> > Cmm(binElements, psiElements,0.);
  arr2<xcomplex<double> > Cmm2(binElements, psiElements,0.);
  arr<double> signarr(lmax+1, 1);
  for (long ii=0; ii<=lmax; ii++)
    if (ii%2==1) signarr[ii]=-1.;

  for (long ii=0; ii<=lmax; ii++)
    {
      double dsignl=signarr[ii];
      long beamupper(beammmax);
      if (ii<beammmax) beamupper=ii;
      arr2<double> dp1;
      wignerCalc(ii, beamupper, beta, dp1);

      long mskyUpper=(ii<effM[ii]) ? ii : effM[ii];
      for(long beamIndex = 0; beamIndex <= beamupper; beamIndex++)
	{
	  double dsignb=signarr[beamIndex], dlb=dsignb*dsignl;
	  for (long msky = 0; msky <= mskyUpper; msky++)
	    {
	      double dmsky=double(msky);
	      double dsign=signarr[msky], dsb=dsign*dsignb;
	      double dsignt=dsb*dsignl;
	      double dMatrixElementmskypos = dsb*dp1[ii+msky][beamIndex+beamupper];//Note that msky in dlm is located to the left of mb and they have to be interchanged in the convolution
	      dlmelements[countdlm]=dMatrixElementmskypos;
	      countdlm++;
	      double dMatrixElementmskyneg = dsb*dp1[ii-msky][beamIndex+beamupper];
	      dlmelements[countdlm]=dMatrixElementmskyneg;
	      countdlm++;
	      double dMatrixElementmskypos2 = dlb*dMatrixElementmskyneg;
	      double dMatrixElementmskyneg2 = dlb*dMatrixElementmskypos;
	      double prod1=slmT(ii,msky).re*blmT(ii,beamIndex).re + slmG(ii,msky).re*blmG(ii,beamIndex).re + slmC(ii,msky).re*blmC(ii,beamIndex).re;
	      double prod3=slmT(ii,msky).im*blmT(ii,beamIndex).re + slmG(ii,msky).im*blmG(ii,beamIndex).re + slmC(ii,msky).im*blmC(ii,beamIndex).re;
	      double prod2=slmT(ii,msky).im*blmT(ii,beamIndex).im + slmG(ii,msky).im*blmG(ii,beamIndex).im + slmC(ii,msky).im*blmC(ii,beamIndex).im;
	      double prod4=slmT(ii,msky).re*blmT(ii,beamIndex).im + slmG(ii,msky).re*blmG(ii,beamIndex).im + slmC(ii,msky).re*blmC(ii,beamIndex).im;
	      double tmp_1=prod1 + prod2;
	      double tmp_2=prod3 - prod4;
	      double tmp_3=prod1 - prod2;
	      double tmp_4=-prod3 - prod4;
	      double xtmp_1=tmp_1 * dMatrixElementmskypos;
	      double xtmp_2=tmp_2 * dMatrixElementmskypos;
	      double xtmp_3=dsign*tmp_3 * dMatrixElementmskyneg;
	      double xtmp_4=dsign*tmp_4 * dMatrixElementmskyneg;
	      double xtmp_5=tmp_1 * dMatrixElementmskypos2;
	      double xtmp_6=tmp_2 * dMatrixElementmskypos2;
	      double xtmp_7=dsign*tmp_3 * dMatrixElementmskyneg2;
	      double xtmp_8=dsign*tmp_4 * dMatrixElementmskyneg2;
	      Cmm[msky+lmaxOut][beamIndex].real() += xtmp_1;
	      Cmm[msky+lmaxOut][beamIndex].imag() += xtmp_2;
	      Cmm2[msky+lmaxOut][beamIndex].real() += xtmp_5;
	      Cmm2[msky+lmaxOut][beamIndex].imag() += xtmp_6;
	      if (msky!=0)
		{
		  Cmm[-msky+lmaxOut][beamIndex].real() += xtmp_3;
		  Cmm[-msky+lmaxOut][beamIndex].imag() += xtmp_4;
		  Cmm2[-msky+lmaxOut][beamIndex].real() += xtmp_7;
		  Cmm2[-msky+lmaxOut][beamIndex].imag() += xtmp_8;
		}
	    }
	}
      dp1.dealloc();
    }
  arr<double> cs(2*lmax+1), sn(2*lmax+1), cs0(binElements), sn0(binElements);
  for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	cs[lmax+msky] =cos(msky*dPhi);
	sn[lmax+msky] =sin(msky*dPhi);
      }
  for (long msky=-lmaxOut; msky<=lmaxOut; msky++)
    {
      cs0[lmaxOut+msky]=cos(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
      sn0[lmaxOut+msky]=sin(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
    }
  todAnnulus(tod1, Cmm, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0);
  todAnnulus(tod2, Cmm2, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0);
}

void conviqt_hemiscm_v4(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, long beammmax, long lmax, long lmaxOut, arr3<xcomplex<double> > &tod1, arr3<xcomplex<double> > &tod2, long NThetaIndex1, arr<double> &rthetas)
{
  wignergen wgen(lmax,rthetas,conv_acc), wgen_neg(lmax,rthetas,conv_acc);
  wigner_estimator estimator(lmax,100);

  long binElements=2*lmaxOut+1;
  long psiElements=beammmax+1;
  arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
  Cmm.fill(0.);
  arr3<xcomplex<double> > Cmm2(binElements, psiElements, NThetaIndex1);
  Cmm2.fill(0.);
  for(long beamIndex = 0; beamIndex <= beammmax; beamIndex++)
    {
      double dbeamIndex = double(beamIndex);
      double dsignb=pow(-1.,dbeamIndex);
      for (long msky = 0; msky <= lmax; msky++)
	{
	  double dmsky=double(msky);
	  double dsign=pow(-1.,dmsky), dsb=dsign*dsignb;
	  estimator.prepare_m (beamIndex,msky);
	  if (estimator.canSkip(rthetas[NThetaIndex1-1])) continue; // negligible dmm
	  wgen.prepare(beamIndex,msky);
	  wgen_neg.prepare(beamIndex,-msky);
	  for (long lat=0; lat<NThetaIndex1; lat++)
	    {
	      int firstl1, firstl2;
	      const arr<double> &dmm=wgen.calc(lat,firstl1);
	      const arr<double> &dmmneg=wgen_neg.calc(lat,firstl2);
	      int firstl=(firstl1>firstl2) ? firstl1: firstl2;
	      for (long ii=firstl; ii<=lmax; ii++)
		{
		  double dii=double(ii);
		  double dsignl=pow(-1.,dii);
		  double dlb=dsignb*dsignl;
		  double dMatrixElementmskypos = dsb*dmm[ii];//Note that msky in dlm is located to the left of mb and they have to be interchanged in the convolution
		  double dMatrixElementmskyneg = dsb*dmmneg[ii];
		  double dMatrixElementmskypos2 = dlb*dMatrixElementmskyneg;
		  double dMatrixElementmskyneg2 = dlb*dMatrixElementmskypos;
		  double prod1=slmT(ii,msky).re*blmT(ii,beamIndex).re;
		  double prod3=slmT(ii,msky).im*blmT(ii,beamIndex).re;
		  double prod2=slmT(ii,msky).im*blmT(ii,beamIndex).im;
		  double prod4=slmT(ii,msky).re*blmT(ii,beamIndex).im;
		  double tmp_1=prod1 + prod2;
		  double tmp_2=prod3 - prod4;
		  double tmp_3=prod1 - prod2;
		  double tmp_4=-prod3 - prod4;
		  double xtmp_1=tmp_1 * dMatrixElementmskypos;
		  double xtmp_2=tmp_2 * dMatrixElementmskypos;
		  double xtmp_3=dsign*tmp_3 * dMatrixElementmskyneg;
		  double xtmp_4=dsign*tmp_4 * dMatrixElementmskyneg;
		  double xtmp_5=tmp_1 * dMatrixElementmskypos2;
		  double xtmp_6=tmp_2 * dMatrixElementmskypos2;
		  double xtmp_7=dsign*tmp_3 * dMatrixElementmskyneg2;
		  double xtmp_8=dsign*tmp_4 * dMatrixElementmskyneg2;
		  Cmm(msky+lmaxOut, beamIndex, lat).real() += xtmp_1;
		  Cmm(msky+lmaxOut, beamIndex, lat).imag() += xtmp_2;
		  Cmm2(msky+lmaxOut, beamIndex, lat).real() += xtmp_5;
		  Cmm2(msky+lmaxOut, beamIndex, lat).imag() += xtmp_6;
		  if (msky!=0)
		    {
		      Cmm(-msky+lmaxOut, beamIndex, lat).real() += xtmp_3;
		      Cmm(-msky+lmaxOut, beamIndex, lat).imag() += xtmp_4;
		      Cmm2(-msky+lmaxOut, beamIndex, lat).real() += xtmp_7;
		      Cmm2(-msky+lmaxOut, beamIndex, lat).imag() += xtmp_8;
		    }
		}
	    }
	}
    }
  arr<double> cs(2*lmax+1), sn(2*lmax+1), cs0(binElements), sn0(binElements);
  for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	cs[lmax+msky] =cos(msky*dPhi);
	sn[lmax+msky] =sin(msky*dPhi);
      }
  for (long msky=-lmaxOut; msky<=lmaxOut; msky++)
    {
      cs0[lmaxOut+msky]=cos(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
      sn0[lmaxOut+msky]=sin(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
    }
  todAnnulus_v3(tod1, Cmm, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0, NThetaIndex1);
  todAnnulus_v3(tod2, Cmm2, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0, NThetaIndex1);
}

void conviqt_hemiscm_pol_v4(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long beammmax, long lmax, long lmaxOut, arr3<xcomplex<double> > &tod1, arr3<xcomplex<double> > &tod2, long NThetaIndex1, arr<double> &rthetas)
{
  wignergen wgen(lmax,rthetas,conv_acc), wgen_neg(lmax,rthetas,conv_acc);
  wigner_estimator estimator(lmax,100);

  long binElements=2*lmaxOut+1;
  long psiElements=beammmax+1;
  arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
  Cmm.fill(0.);
  arr3<xcomplex<double> > Cmm2(binElements, psiElements, NThetaIndex1);
  Cmm2.fill(0.);
  for(long beamIndex = 0; beamIndex <= beammmax; beamIndex++)
    {
      double dbeamIndex = double(beamIndex);
      double dsignb=xpow(beamIndex,1);//pow(-1.,dbeamIndex);
      for (long msky = 0; msky <= lmax; msky++)
	{
	  double dmsky=double(msky);
	  double dsign=xpow(msky,1);//pow(-1.,dmsky), 
	  double dsb=dsign*dsignb;
	  estimator.prepare_m (beamIndex,msky);
	  if (estimator.canSkip(rthetas[NThetaIndex1-1])) continue; // negligible dmm
	  wgen.prepare(beamIndex,msky);
	  wgen_neg.prepare(beamIndex,-msky);
	  for (long lat=0; lat<NThetaIndex1; lat++)
	    {
	      int firstl1, firstl2;
	      const arr<double> &dmm=wgen.calc(lat,firstl1);
	      const arr<double> &dmmneg=wgen_neg.calc(lat,firstl2);
	      int firstl=(firstl1>firstl2) ? firstl1: firstl2;
	      for (long ii=firstl; ii<=lmax; ii++)
		{
		  double dii=double(ii);
		  double dsignl=xpow(ii,1);//pow(-1.,dii);
		  double dlb=dsignb*dsignl;
		  double dMatrixElementmskypos = dsb*dmm[ii];//Note that msky in dlm is located to the left of mb and they have to be interchanged in the convolution
		  double dMatrixElementmskyneg = dsb*dmmneg[ii];
		  double dMatrixElementmskypos2 = dlb*dMatrixElementmskyneg;
		  double dMatrixElementmskyneg2 = dlb*dMatrixElementmskypos;
		  double prod1=slmT(ii,msky).re*blmT(ii,beamIndex).re + slmG(ii,msky).re*blmG(ii,beamIndex).re + slmC(ii,msky).re*blmC(ii,beamIndex).re;
		  double prod3=slmT(ii,msky).im*blmT(ii,beamIndex).re + slmG(ii,msky).im*blmG(ii,beamIndex).re + slmC(ii,msky).im*blmC(ii,beamIndex).re;
		  double prod2=slmT(ii,msky).im*blmT(ii,beamIndex).im + slmG(ii,msky).im*blmG(ii,beamIndex).im + slmC(ii,msky).im*blmC(ii,beamIndex).im;
		  double prod4=slmT(ii,msky).re*blmT(ii,beamIndex).im + slmG(ii,msky).re*blmG(ii,beamIndex).im + slmC(ii,msky).re*blmC(ii,beamIndex).im;
		  double tmp_1=prod1 + prod2;
		  double tmp_2=prod3 - prod4;
		  double tmp_3=prod1 - prod2;
		  double tmp_4=-prod3 - prod4;
		  double xtmp_1=tmp_1 * dMatrixElementmskypos;
		  double xtmp_2=tmp_2 * dMatrixElementmskypos;
		  double xtmp_3=dsign*tmp_3 * dMatrixElementmskyneg;
		  double xtmp_4=dsign*tmp_4 * dMatrixElementmskyneg;
		  double xtmp_5=tmp_1 * dMatrixElementmskypos2;
		  double xtmp_6=tmp_2 * dMatrixElementmskypos2;
		  double xtmp_7=dsign*tmp_3 * dMatrixElementmskyneg2;
		  double xtmp_8=dsign*tmp_4 * dMatrixElementmskyneg2;
		  Cmm(msky+lmaxOut, beamIndex, lat).real() += xtmp_1;
		  Cmm(msky+lmaxOut, beamIndex, lat).imag() += xtmp_2;
		  Cmm2(msky+lmaxOut, beamIndex, lat).real() += xtmp_5;
		  Cmm2(msky+lmaxOut, beamIndex, lat).imag() += xtmp_6;
		  if (msky!=0)
		    {
		      Cmm(-msky+lmaxOut, beamIndex, lat).real() += xtmp_3;
		      Cmm(-msky+lmaxOut, beamIndex, lat).imag() += xtmp_4;
		      Cmm2(-msky+lmaxOut, beamIndex, lat).real() += xtmp_7;
		      Cmm2(-msky+lmaxOut, beamIndex, lat).imag() += xtmp_8;
		    }
		}
	    }
	}
    }
  arr<double> cs(2*lmax+1), sn(2*lmax+1), cs0(binElements), sn0(binElements);
  for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	cs[lmax+msky] =cos(msky*dPhi);
	sn[lmax+msky] =sin(msky*dPhi);
      }
  for (long msky=-lmaxOut; msky<=lmaxOut; msky++)
    {
      cs0[lmaxOut+msky]=cos(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
      sn0[lmaxOut+msky]=sin(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
    }
  todAnnulus_v3(tod1, Cmm, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0, NThetaIndex1);
  todAnnulus_v3(tod2, Cmm2, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0, NThetaIndex1);
}

void conviqt_hemiscm_single(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, long beammmax, long lmax, long lmaxOut, arr3<xcomplex<double> > &tod1, long NThetaIndex1, arr<double> &rthetas)
{
  wignergen wgen(lmax,rthetas,conv_acc), wgen_neg(lmax,rthetas,conv_acc);
  wigner_estimator estimator(lmax,100);

  long binElements=2*lmaxOut+1;
  long psiElements=beammmax+1;
  arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
  Cmm.fill(0.);
  for(long beamIndex = 0; beamIndex <= beammmax; beamIndex++)
    {
      double dbeamIndex = double(beamIndex);
      double dsignb=pow(-1.,dbeamIndex);
      for (long msky = 0; msky <= lmax; msky++)
	{
	  double dmsky=double(msky);
	  double dsign=pow(-1.,dmsky), dsb=dsign*dsignb;
	  estimator.prepare_m (beamIndex,msky);
	  if (estimator.canSkip(rthetas[NThetaIndex1-1])) continue; // negligible dmm
	  wgen.prepare(beamIndex,msky);
	  wgen_neg.prepare(beamIndex,-msky);
	  for (long lat=0; lat<NThetaIndex1; lat++)
	    {
	      int firstl1, firstl2;
	      const arr<double> &dmm=wgen.calc(lat,firstl1);
	      const arr<double> &dmmneg=wgen_neg.calc(lat,firstl2);
	      int firstl=(firstl1>firstl2) ? firstl1: firstl2;
	      for (long ii=firstl; ii<=lmax; ii++)
		{
		  double dii=double(ii);
		  double dsignl=pow(-1.,dii);
		  double dlb=dsignb*dsignl;
		  double dMatrixElementmskypos = dsb*dmm[ii];//Note that msky in dlm is located to the left of mb and they have to be interchanged in the convolution
		  double dMatrixElementmskyneg = dsb*dmmneg[ii];
		  double dMatrixElementmskypos2 = dlb*dMatrixElementmskyneg;
		  double dMatrixElementmskyneg2 = dlb*dMatrixElementmskypos;
		  double prod1=slmT(ii,msky).re*blmT(ii,beamIndex).re;
		  double prod3=slmT(ii,msky).im*blmT(ii,beamIndex).re;
		  double prod2=slmT(ii,msky).im*blmT(ii,beamIndex).im;
		  double prod4=slmT(ii,msky).re*blmT(ii,beamIndex).im;
		  double tmp_1=prod1 + prod2;
		  double tmp_2=prod3 - prod4;
		  double tmp_3=prod1 - prod2;
		  double tmp_4=-prod3 - prod4;
		  double xtmp_1=tmp_1 * dMatrixElementmskypos;
		  double xtmp_2=tmp_2 * dMatrixElementmskypos;
		  double xtmp_3=dsign*tmp_3 * dMatrixElementmskyneg;
		  double xtmp_4=dsign*tmp_4 * dMatrixElementmskyneg;
		  Cmm(msky+lmaxOut, beamIndex, lat).real() += xtmp_1;
		  Cmm(msky+lmaxOut, beamIndex, lat).imag() += xtmp_2;
		  if (msky!=0)
		    {
		      Cmm(-msky+lmaxOut, beamIndex, lat).real() += xtmp_3;
		      Cmm(-msky+lmaxOut, beamIndex, lat).imag() += xtmp_4;
		    }
		}
	    }
	}
    }
  arr<double> cs(2*lmax+1), sn(2*lmax+1), cs0(binElements), sn0(binElements);
  for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	cs[lmax+msky] =cos(msky*dPhi);
	sn[lmax+msky] =sin(msky*dPhi);
      }
  for (long msky=-lmaxOut; msky<=lmaxOut; msky++)
    {
      cs0[lmaxOut+msky]=cos(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
      sn0[lmaxOut+msky]=sin(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
    }
  todAnnulus_v3(tod1, Cmm, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0, NThetaIndex1);
}

void conviqt_hemiscm_pol_single(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long beammmax, long lmax, long lmaxOut, arr3<xcomplex<double> > &tod1, long NThetaIndex1, arr<double> &rthetas)
{
  wignergen wgen(lmax,rthetas,conv_acc), wgen_neg(lmax,rthetas,conv_acc);
  wigner_estimator estimator(lmax,100);

  long binElements=2*lmaxOut+1;
  long psiElements=beammmax+1;
  arr3<xcomplex<double> > Cmm(binElements, psiElements, NThetaIndex1);
  Cmm.fill(0.);
  for(long beamIndex = 0; beamIndex <= beammmax; beamIndex++)
    {
      double dbeamIndex = double(beamIndex);
      double dsignb=pow(-1.,dbeamIndex);
      for (long msky = 0; msky <= lmax; msky++)
	{
	  double dmsky=double(msky);
	  double dsign=pow(-1.,dmsky), dsb=dsign*dsignb;
	  estimator.prepare_m (beamIndex,msky);
	  if (estimator.canSkip(rthetas[NThetaIndex1-1])) continue; // negligible dmm
	  wgen.prepare(beamIndex,msky);
	  wgen_neg.prepare(beamIndex,-msky);
	  for (long lat=0; lat<NThetaIndex1; lat++)
	    {
	      int firstl1, firstl2;
	      const arr<double> &dmm=wgen.calc(lat,firstl1);
	      const arr<double> &dmmneg=wgen_neg.calc(lat,firstl2);
	      int firstl=(firstl1>firstl2) ? firstl1: firstl2;
	      for (long ii=firstl; ii<=lmax; ii++)
		{
		  double dii=double(ii);
		  double dsignl=pow(-1.,dii);
		  double dlb=dsignb*dsignl;
		  double dMatrixElementmskypos = dsb*dmm[ii];//Note that msky in dlm is located to the left of mb and they have to be interchanged in the convolution
		  double dMatrixElementmskyneg = dsb*dmmneg[ii];
		  double dMatrixElementmskypos2 = dlb*dMatrixElementmskyneg;
		  double dMatrixElementmskyneg2 = dlb*dMatrixElementmskypos;
		  double prod1=slmT(ii,msky).re*blmT(ii,beamIndex).re + slmG(ii,msky).re*blmG(ii,beamIndex).re + slmC(ii,msky).re*blmC(ii,beamIndex).re;
		  double prod3=slmT(ii,msky).im*blmT(ii,beamIndex).re + slmG(ii,msky).im*blmG(ii,beamIndex).re + slmC(ii,msky).im*blmC(ii,beamIndex).re;
		  double prod2=slmT(ii,msky).im*blmT(ii,beamIndex).im + slmG(ii,msky).im*blmG(ii,beamIndex).im + slmC(ii,msky).im*blmC(ii,beamIndex).im;
		  double prod4=slmT(ii,msky).re*blmT(ii,beamIndex).im + slmG(ii,msky).re*blmG(ii,beamIndex).im + slmC(ii,msky).re*blmC(ii,beamIndex).im;
		  double tmp_1=prod1 + prod2;
		  double tmp_2=prod3 - prod4;
		  double tmp_3=prod1 - prod2;
		  double tmp_4=-prod3 - prod4;
		  double xtmp_1=tmp_1 * dMatrixElementmskypos;
		  double xtmp_2=tmp_2 * dMatrixElementmskypos;
		  double xtmp_3=dsign*tmp_3 * dMatrixElementmskyneg;
		  double xtmp_4=dsign*tmp_4 * dMatrixElementmskyneg;
		  Cmm(msky+lmaxOut, beamIndex, lat).real() += xtmp_1;
		  Cmm(msky+lmaxOut, beamIndex, lat).imag() += xtmp_2;
		  if (msky!=0)
		    {
		      Cmm(-msky+lmaxOut, beamIndex, lat).real() += xtmp_3;
		      Cmm(-msky+lmaxOut, beamIndex, lat).imag() += xtmp_4;
		    }
		  // DEBUG
		  if ( beamIndex==0 && msky==0 && lat==0 )
		    {
		      std::cerr << "    " << slmT(ii,msky) << " " << slmG(ii,msky) << " "  << slmC(ii,msky) << std::endl;
		      std::cerr << "    " << blmT(ii,beamIndex) << " " << blmG(ii,beamIndex) << " "  << blmC(ii,beamIndex) << std::endl;
 		      std::cerr << ii << " " << dMatrixElementmskypos << " " << dMatrixElementmskyneg << " " << prod1 << " " << prod3 << " " << prod2 << " " << prod4 << " " << tmp_1 << " " << tmp_2 << " " << tmp_3 << " " << tmp_4 << " " << xtmp_1 << " " << xtmp_2 << " " << xtmp_3 << " " << xtmp_4 << std::endl;
		    }
		  // end DEBUG
		}
	    }
	}
    }

  // DEBUG
  //if ( corenum == 0 )
  {
    std::cerr << "Cmm: " << std::endl;
    for ( int i=0; i<3; ++i )
      for ( int j=0; j<3; ++j )
	for ( int k=0; k<3; ++k )
	  std::cerr << Cmm(i,j,k) << std::endl;
    std::cerr << std::endl;
  }
  // end DEBUG

  arr<double> cs(2*lmax+1), sn(2*lmax+1), cs0(binElements), sn0(binElements);
  for (long msky = -lmax; msky <= lmax; msky++)
      {
	double dPhi=-halfpi;
	cs[lmax+msky] =cos(msky*dPhi);
	sn[lmax+msky] =sin(msky*dPhi);
      }
  for (long msky=-lmaxOut; msky<=lmaxOut; msky++)
    {
      cs0[lmaxOut+msky]=cos(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
      sn0[lmaxOut+msky]=sin(2.*pi*(msky+lmaxOut)*lmaxOut/(2.*lmaxOut+1.));
    }
  todAnnulus_v3(tod1, Cmm, lmax, beammmax, lmaxOut, cs, sn, cs0, sn0, NThetaIndex1);
}

void ithetacalc(arr<long> &itheta0, arr<double> &outpntarr, long ntod, double inv_delta_theta, double theta0, long ioffset, long ntheta, long npoints)
{
  for (long ii=0; ii<ntod; ii++)
    {
      double beta=outpntarr[5*ii+1];
      beta=(beta>halfpi) ? pi-beta : beta;
      double frac = (beta-theta0)*inv_delta_theta;//Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
      itheta0[ii] = int (frac) - ioffset;//Note also that itheta0[ii] is always positive.
      if (itheta0[ii]>(ntheta-npoints)) itheta0[ii] = ntheta-npoints;
      if (itheta0[ii]<0) itheta0[ii] = 0;
    }
}

void conviqt_tod_loop_fast(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr2<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long order, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long beammmax, long npoints, long ntod, arr<double> &TODValue)
{
  for (int ii=lowerIndex[thetaIndex]; ii>=upperIndex[thetaIndex]; ii--)
    {
      arr<double> wgt1(max_order+1,0.);
      
      double frac = (outpntarr[5*ii+1]-theta0)*inv_delta_theta;//Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
      frac -= itheta0[ii];
      weight_ncm (frac,wgt1, order);
      
      arr<double> wgt2(max_order+1,0.);
      frac = outpntarr[5*ii]*inv_delta_phi - phioffset;
      frac = fmodulo (frac,double(nphi));
      int iphi0 = int (frac) - ioffset;
      frac -= iphi0;
      if (iphi0 >= nphi) iphi0-=nphi;
      if (iphi0 < 0) iphi0+=nphi;
      weight_ncm (frac,wgt2,order);
      
      double omega=outpntarr[5*ii+2]-halfpi;
      double sinomg=sin(omega), cosomg=cos(omega);
      arr<double> cosang(beammmax+1), sinang(beammmax+1);
      cosang[0]=1.; sinang[0]=0.;
      double weight = 2*wgt2[0]*wgt1[thetaIndex-(itheta0[ii]-itheta0[ntod-1])];
      for (long psiIndex=1; psiIndex<=beammmax; psiIndex++) 
	{
	  cosang[psiIndex]=weight*(cosang[psiIndex-1]*cosomg - sinang[psiIndex-1]*sinomg);
	  sinang[psiIndex]=weight*(sinang[psiIndex-1]*cosomg + cosang[psiIndex-1]*sinomg);
	}
      for (int phiIndex=0; phiIndex<npoints; ++phiIndex)
	{
	  //	  double weight = wgt2[phiIndex]*wgt1[thetaIndex-(itheta0[ii]-itheta0[ntod-1])];
	  long newphiIndex=iphi0+phiIndex;
	  if (newphiIndex>=nphi) newphiIndex-=nphi;
	  TODValue[ii]+=0.5*weight*TODAsym[newphiIndex][0].real();
	  for (long psiIndex=1; psiIndex<=beammmax; psiIndex++) 
	    {
	      TODValue[ii]+=( cosang[psiIndex]*TODAsym[newphiIndex][psiIndex].real() - sinang[psiIndex]*TODAsym[newphiIndex][psiIndex].imag() );
	    }
	}
    }
}

void conviqt_tod_loop_read(arr<long> &lowerIndex, arr<long> &upperIndex, arr2<xcomplex<double> > &TODAsym, long thetaIndex, long nphi, long beammmax, long npoints, arr<double> &TODValue, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight)
{
  for (int ii=lowerIndex[thetaIndex]; ii>=upperIndex[thetaIndex]; ii--)
    for (int phiIndex=0; phiIndex<npoints; ++phiIndex)
      {
	long newphiIndex=iphi0arr[ii]+phiIndex;
	if (newphiIndex>=nphi) newphiIndex-=nphi;
	for (long psiIndex=0; psiIndex<=beammmax; psiIndex++) TODValue[ii]+=cosweight[ii][psiIndex]*TODAsym[newphiIndex][psiIndex].real() - sinweight[ii][psiIndex]*TODAsym[newphiIndex][psiIndex].imag();
      }
}

void interpolTOD_arrTestcm_pol_fast(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long lmax, long lmaxOut, long beammmax, long order, arr<double> &TODValue, arr<double> &TODValue2, int corenum, arr<double> &betaIni_arr, long &NThetaIndex1, arr<int> &iphi0arr, arr<int> &iphi0arr2, arr2<double> &sinweight, arr2<double> &cosweight, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<long> &effM, arr<double> &dp1, int corechunk, long newcorenum, paramfile &par_file, bool pol, long coreupper)
{
  if ( CMULT_VERBOSITY > 1 ) cout << "Entered interpolTOD_arrTestcm_pol_fast"<< endl;
  long nphi=(2*lmaxOut+1);
  int npoints=order+1;

  arr2<xcomplex<double> > TODAsym(nphi,beammmax+1,0.), TODAsym2(nphi,beammmax+1,0.);
  double elapsed_secs=0.;
  long countdlm=0;
  if (corenum >= 2*corechunk) readinslm(corenum, newcorenum, coreupper, slmT, slmG, slmC, lmax, pol, par_file);
  for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex)
    {
      double beta=betaIni_arr[thetaIndex];
      clock_t begin = clock();
      conviqt_hemiscm_pol_fast(blmT, slmT, blmG, slmG, blmC, slmC, beammmax, lmax, lmaxOut, TODAsym, TODAsym2, beta, countdlm, effM, dp1, corenum);
      clock_t end = clock();
      elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
      conviqt_tod_loop_read(lowerIndex, upperIndex, TODAsym, thetaIndex, nphi, beammmax, npoints, TODValue, iphi0arr, sinweight, cosweight);
      conviqt_tod_loop_read(lowerIndex2, upperIndex2, TODAsym2, thetaIndex, nphi, beammmax, npoints, TODValue2, iphi0arr2, sinweight2, cosweight2);
    }
  if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "  elapsed_secs for conviqt2 = " << elapsed_secs << "   average = " << elapsed_secs/NThetaIndex1 << endl; 
  if ( CMULT_VERBOSITY > 1 ) cout << "Leaving interpolTOD_arrTestcm_pol_fast"<< endl;
}



void conviqt_tod_loop_fill(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr2<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long order, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long beammmax, long npoints, long ntod, arr<double> &TODValue, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight)
{
  for (int ii=lowerIndex[thetaIndex]; ii>=upperIndex[thetaIndex]; ii--)
    {
      arr<double> wgt1(max_order+1,0.);
      
      double frac = (outpntarr[5*ii+1]-theta0)*inv_delta_theta;//Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
      frac -= itheta0[ii];
      weight_ncm (frac,wgt1, order);
      
      arr<double> wgt2(max_order+1,0.);
      frac = outpntarr[5*ii]*inv_delta_phi - phioffset;
      frac = fmodulo (frac,double(nphi));
      int iphi0 = int (frac) - ioffset;
      frac -= iphi0;
      if (iphi0 >= nphi) iphi0-=nphi;
      if (iphi0 < 0) iphi0+=nphi;
      iphi0arr[ii]=iphi0;
      weight_ncm (frac,wgt2,order);
      
      double omega=outpntarr[5*ii+2]-halfpi;
      double sinomg=sin(omega), cosomg=cos(omega);
      arr<double> cosang(beammmax+1), sinang(beammmax+1);
      cosang[0]=1.; sinang[0]=0.;
      for (long psiIndex=1; psiIndex<=beammmax; psiIndex++) 
	{
	  cosang[psiIndex]=cosang[psiIndex-1]*cosomg - sinang[psiIndex-1]*sinomg;
	  sinang[psiIndex]=sinang[psiIndex-1]*cosomg + cosang[psiIndex-1]*sinomg;
	}
      cosang[0]=0.5;
      for (int phiIndex=0; phiIndex<npoints; ++phiIndex)
	{
	  double weight = wgt2[phiIndex]*wgt1[thetaIndex-(itheta0[ii]-itheta0[ntod-1])];
	  long newphiIndex=iphi0+phiIndex;
	  if (newphiIndex>=nphi) newphiIndex-=nphi;
	  for (long psiIndex=0; psiIndex<=beammmax; psiIndex++) 
	    {
	      cosweight[ii][psiIndex]=2.*weight*cosang[psiIndex];
	      sinweight[ii][psiIndex]=2.*weight*sinang[psiIndex];
	      TODValue[ii]+=cosweight[ii][psiIndex]*TODAsym[newphiIndex][psiIndex].real() - sinweight[ii][psiIndex]*TODAsym[newphiIndex][psiIndex].imag();
	    }
	}
    }
}

void conviqt_tod_loop_v4(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr3<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long order, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long beammmax, long npoints, long ntod, arr<double> &TODValue, long lat)
{
  for (int ii=lowerIndex[thetaIndex]; ii>=upperIndex[thetaIndex]; ii--)
    {
      arr<double> wgt1(max_order+1,0.);
      
      double frac = (outpntarr[5*ii+1]-theta0)*inv_delta_theta;//Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.
      frac -= itheta0[ii];
      weight_ncm (frac,wgt1, order);
      
      arr<double> wgt2(max_order+1,0.);
      frac = outpntarr[5*ii]*inv_delta_phi - phioffset;
      frac = fmodulo (frac,double(nphi));
      int iphi0 = int (frac) - ioffset;
      frac -= iphi0;
      if (iphi0 >= nphi) iphi0-=nphi;
      if (iphi0 < 0) iphi0+=nphi;
      weight_ncm (frac,wgt2,order);
      
      double omega=outpntarr[5*ii+2]+halfpi;
      double sinomg=sin(omega), cosomg=cos(omega);
      arr<double> cosang(beammmax+1), sinang(beammmax+1);
      cosang[0]=1.; sinang[0]=0.;
      for (long psiIndex=1; psiIndex<=beammmax; psiIndex++) 
	{
	  cosang[psiIndex]=cosang[psiIndex-1]*cosomg - sinang[psiIndex-1]*sinomg;
	  sinang[psiIndex]=sinang[psiIndex-1]*cosomg + cosang[psiIndex-1]*sinomg;
	}
      cosang[0]=0.5;
      for (int phiIndex=0; phiIndex<npoints; ++phiIndex)
	{
	  double weight = wgt2[phiIndex]*wgt1[thetaIndex-(itheta0[ii]-itheta0[ntod-1])];
	  long newphiIndex=iphi0+phiIndex;
	  if (newphiIndex>=nphi) newphiIndex-=nphi;
	  for (long psiIndex=0; psiIndex<=beammmax; psiIndex++) 
	    {
	      TODValue[ii]+=2.*weight*( cosang[psiIndex]*TODAsym(newphiIndex, psiIndex, lat).real() - sinang[psiIndex]*TODAsym(newphiIndex, psiIndex, lat).imag() );
	    }
	}
    }
}

void conviqt_tod_loop_pol_v5(arr<long> &lowerIndex, arr<long> &upperIndex, arr<double> &outpntarr, arr3<xcomplex<double> > &TODAsym, long thetaIndex, arr<long> &itheta0, long order, long max_order, double inv_delta_theta, double theta0, double inv_delta_phi, double phioffset, long nphi, long ioffset, long beammmax, long npoints, long ntod, arr<double> &TODValue, long lat)
{
  arr2<xcomplex<double> > conviqtarr(nphi,beammmax+1);
  for (long ii=0; ii<nphi; ii++)
    for (long jj=0; jj<beammmax+1; jj++)
      conviqtarr[ii][jj]=TODAsym(ii, jj, lat);
  for (int ii=lowerIndex[thetaIndex]; ii>=upperIndex[thetaIndex]; ii--)
    {
      arr<double> wgt1(max_order+1,0.);

      double frac = (outpntarr[5*ii+1]-theta0)*inv_delta_theta;//Note that the larger the ii the smaller frac is and the smaller itheta0[ii] is.                                    
      frac -= itheta0[ii];
      weight_ncm (frac,wgt1, order);

      arr<double> wgt2(max_order+1,0.);
      frac = outpntarr[5*ii]*inv_delta_phi - phioffset;
      frac = fmodulo (frac,double(nphi));
      int iphi0 = int (frac) - ioffset;
      frac -= iphi0;
      if (iphi0 >= nphi) iphi0-=nphi;
      if (iphi0 < 0) iphi0+=nphi;
      weight_ncm (frac,wgt2,order);

      double omega=outpntarr[5*ii+2]+halfpi;
      double sinomg=sin(omega), cosomg=cos(omega);
      arr<double> cosang(beammmax+1), sinang(beammmax+1);
      cosang[0]=1.; sinang[0]=0.;
      for (long psiIndex=1; psiIndex<=beammmax; psiIndex++)
        {
          cosang[psiIndex]=cosang[psiIndex-1]*cosomg - sinang[psiIndex-1]*sinomg;
          sinang[psiIndex]=sinang[psiIndex-1]*cosomg + cosang[psiIndex-1]*sinomg;
        }
      cosang[0]=0.5;
      for (int phiIndex=0; phiIndex<npoints; ++phiIndex)
        {
          double weight = wgt2[phiIndex]*wgt1[thetaIndex-(itheta0[ii]-itheta0[ntod-1])];
          long newphiIndex=iphi0+phiIndex;
          if (newphiIndex>=nphi) newphiIndex-=nphi;
          for (long psiIndex=0; psiIndex<=beammmax; psiIndex++)
            {
              TODValue[ii]+=2.*weight*( cosang[psiIndex]*conviqtarr[newphiIndex][psiIndex].real() - sinang[psiIndex]*conviqtarr[newphiIndex][psiIndex].imag() );
            }
        }
    }
}

void effMFiller(double beta, arr<long> &effM, long beammmax, long lmax)
{
  double eps=1.e-5;
  for (long ii=0; ii<=lmax; ii++)
    {
      long beamupper(beammmax);
      if (ii<beammmax) beamupper=ii;
      arr2<double> dp1;
      wignerCalc(ii, beamupper, beta, dp1);
      effM[ii]=1e6;
      for (long mm=ii; mm>=beamupper; mm--) if ( abs(dp1[ii+mm][2*beamupper]) > eps ) {effM[ii]=mm; break;}
      if (effM[ii]==1e6) effM[ii]=0;
      dp1.dealloc();
    }
}

void arrsizecounter(arr<long> &effM, long &numberofdlms, long beammmax, long lmax)
{
  for (long ii=0; ii<=lmax; ii++)
    {
      long beamupper(beammmax);
      if (ii<beammmax) beamupper=ii;
      long mskyUpper=(ii<effM[ii]) ? ii : effM[ii];
      for(long beamIndex = 0; beamIndex <= beamupper; beamIndex++)
	for (long msky = 0; msky <= mskyUpper; msky++) 
	  {
	    numberofdlms+=2;
	  }
    }
}

void interpolTOD_arrTestcm_pol_fill(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long lmax, long lmaxOut, long beammmax, arr<double> &outpntarr1, arr<double> &outpntarr2, long order, arr<double> &TODValue, arr<double> &TODValue2, long ntod, long ntod2, int corenum, int cores, arr<double> &betaIni_arr, long &NThetaIndex1, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight, arr<int> &iphi0arr2, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &effM, arr<double> &dp1)
{
  if ( CMULT_VERBOSITY > 1 ) cout << "Entered interpolTOD_arrTestcm_pol_fill"<< endl;
  double phi0=halfpi;
  long npsi=2*beammmax+1;
  long nphi=(2*lmaxOut+1);
  double dphi=2.*pi/nphi;
  double inv_delta_phi = 1./dphi;
  double phioffset = phi0/dphi;
  long ntheta=lmaxOut+1+21;
  double dtheta=-pi/(ntheta-21);
  double theta0=pi-10*dtheta;
  double inv_delta_theta = 1./dtheta;
  long max_order=19;
  int npoints=order+1;
  int ioffset = order/2;
  arr<long> itheta0(ntod), itheta0_2(ntod2);

  ithetacalc(itheta0, outpntarr1, ntod, inv_delta_theta, theta0, ioffset, ntheta, npoints);
  ithetacalc(itheta0_2, outpntarr2, ntod2, inv_delta_theta, theta0, ioffset, ntheta, npoints);

  NThetaIndex1=itheta0[0]-itheta0[ntod-1]+npoints;
  long NThetaIndex2=itheta0_2[0]-itheta0_2[ntod2-1]+npoints;

  if (itheta0[ntod-1] != itheta0_2[ntod2-1]) { if ( CMULT_VERBOSITY > 1 ) cout << "itheta0[ntod-1] = " << itheta0[ntod-1] << "   itheta0_2[ntod2-1] = " << itheta0_2[ntod2-1] << "   itheta0[0] = " << itheta0[0] << "   itheta0_2[0] = " << itheta0_2[0] << endl; planck_fail("itheta0[ntod-1] != itheta0_2[ntod2-1]");}
  if (NThetaIndex1 != NThetaIndex2) { if ( CMULT_VERBOSITY > 1 ) cout << "NThetaIndex1 = " << NThetaIndex1 << "   NThetaIndex2 = " << NThetaIndex2 << endl; planck_fail("NThetaIndex1 not equal to NThetaIndex2");}

  arr<long> inNThetaIndex1(1);
  arr<double> incorethetas;
  if (corenum!=cores-1) mpiMgr.sendRaw(&NThetaIndex1, 1, corenum+1);
  if (corenum != 0) mpiMgr.recvRaw (&inNThetaIndex1[0], 1, corenum-1);

  arr<double> corethetas(NThetaIndex1);
  if (corenum!=0) incorethetas.alloc(inNThetaIndex1[0]);
  for (long ii=0; ii<NThetaIndex1; ii++) corethetas[ii]=theta0+(itheta0[ntod-1]+ii)*dtheta;

  long basethetaIndex=NThetaIndex1;
  if (corenum!=cores-1) mpiMgr.sendRaw(&corethetas[0], NThetaIndex1, corenum+1);
  if (corenum != 0) 
    {
      mpiMgr.recvRaw (&incorethetas[0], inNThetaIndex1[0], corenum-1);
      for (long ii = 0; ii<NThetaIndex1; ii++) 
	if (corethetas[ii] <= incorethetas[0])
	  {
	    basethetaIndex=ii;
	    break;
	  }
    }
  if (NThetaIndex1>1) if ( CMULT_VERBOSITY > 1 ) cout << "corethetas[0] = " << corethetas[0] << "  corethetas[1] = " << corethetas[1] << endl;
  
  lowerIndex.alloc(NThetaIndex1);
  lowerIndex.fill(ntod);
  upperIndex.alloc(NThetaIndex1);
  upperIndex.fill(0);
  for (long jj=0; jj<NThetaIndex1; jj++)
    for (long ii=ntod-1; ii>=0; ii--) if (itheta0[ii]>itheta0[ntod-1]+jj-npoints) {lowerIndex[jj]=ii; break;}
  for (long jj=0; jj<NThetaIndex1; jj++)
    for (long ii=0; ii<ntod; ii++) if (itheta0[ii]<=itheta0[ntod-1]+jj) {upperIndex[jj]=ii; break;}//this yields the smallest ii where itheta0[ii]<=itheta0[ntod-1]+jj

  lowerIndex2.alloc(NThetaIndex1);
  lowerIndex2.fill(ntod2);
  upperIndex2.alloc(NThetaIndex1);
  upperIndex2.fill(0);
  for (long jj=0; jj<NThetaIndex1; jj++)
    for (long ii=ntod2-1; ii>=0; ii--) if (itheta0_2[ii]>itheta0_2[ntod2-1]+jj-npoints) {lowerIndex2[jj]=ii; break;}
  for (long jj=0; jj<NThetaIndex1; jj++)
    for (long ii=0; ii<ntod2; ii++) if (itheta0_2[ii]<=itheta0_2[ntod2-1]+jj) {upperIndex2[jj]=ii; break;}//this yields the smallest ii where itheta0_2[ii]<=itheta0_2[ntod-1]+jj

  arr2<xcomplex<double> > TODAsym(nphi,beammmax+1,0.), TODAsym2(nphi,beammmax+1,0.);

  double elapsed_secs=0.;
  betaIni_arr.alloc(NThetaIndex1);
  try {
    iphi0arr.alloc(ntod), sinweight.alloc(ntod,beammmax+1), cosweight.alloc(ntod,beammmax+1);
    iphi0arr2.alloc(ntod2), sinweight2.alloc(ntod2,beammmax+1), cosweight2.alloc(ntod2,beammmax+1);
  } catch ( bad_alloc & e ) {
    cerr << "interpolTOD_arrTestcm_pol_fill : Out of memory allocating " << (ntod+ntod2)*(2*beammmax+3)*8./1024/1024 << "MB for pointing weights" << endl;
    throw;
  }

  double betatmp=theta0+(itheta0[ntod-1]+NThetaIndex1-1)*dtheta;  
  effMFiller(betatmp, effM, beammmax, lmax);

  long numberofdlms=0;
  for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex) arrsizecounter(effM, numberofdlms, beammmax, lmax);
  dp1.alloc(numberofdlms);
  long dp1count=0;
  for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex)
    {
      double beta=theta0+(itheta0[ntod-1]+thetaIndex)*dtheta;
      betaIni_arr[thetaIndex]=beta;
      if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "   thetaIndex = " << thetaIndex << "   beta = " << beta << "   NThetaIndex1 = " << NThetaIndex1 << "  basethetaIndex = " << basethetaIndex << endl;
      clock_t begin = clock();
      conviqt_hemiscm_pol_fill(blmT, slmT, blmG, slmG, blmC, slmC, beammmax, lmax, lmaxOut, TODAsym, TODAsym2, beta, effM, dp1count, corenum, dp1);
      clock_t end = clock();
      elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;
      conviqt_tod_loop_fill(lowerIndex, upperIndex, outpntarr1, TODAsym, thetaIndex, itheta0, order, max_order, inv_delta_theta, theta0, inv_delta_phi, phioffset, nphi, ioffset, beammmax, npoints, ntod, TODValue, iphi0arr, sinweight, cosweight);
      conviqt_tod_loop_fill(lowerIndex2, upperIndex2, outpntarr2, TODAsym2, thetaIndex, itheta0_2, order, max_order, inv_delta_theta, theta0, inv_delta_phi, phioffset, nphi, ioffset, beammmax, npoints, ntod2, TODValue2, iphi0arr2, sinweight2, cosweight2);
    }
  if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "  elapsed_secs for conviqt = " << elapsed_secs << "   average = " << elapsed_secs/NThetaIndex1 << "  effM[lmax]= " << effM[lmax]<< endl; 
  if ( CMULT_VERBOSITY > 1 ) cout << "Leaving interpolTOD_arrTestcm_pol_fill"<< endl;
}


void itheta0SetUp(int npoints, int ioffset, long ntheta, double theta0, double inv_delta_theta, long nphi, long beammmax, arr<double> outpntarr, long ntod, long &NThetaIndex, arr<long> &itheta0, arr<long> &lowerIndex, arr<long> &upperIndex, arr3<xcomplex<double> > &TODAsym)
{
  itheta0.alloc(ntod);
  ithetacalc(itheta0, outpntarr, ntod, inv_delta_theta, theta0, ioffset, ntheta, npoints);
  NThetaIndex=itheta0[0]-itheta0[ntod-1]+npoints;
  lowerIndex.alloc(NThetaIndex);
  lowerIndex.fill(ntod);
  upperIndex.alloc(NThetaIndex);
  upperIndex.fill(0);
  for (long jj=0; jj<NThetaIndex; jj++)
    for (long ii=ntod-1; ii>=0; ii--) if (itheta0[ii]>itheta0[ntod-1]+jj-npoints) {lowerIndex[jj]=ii; break;}
  for (long jj=0; jj<NThetaIndex; jj++)
    for (long ii=0; ii<ntod; ii++) if (itheta0[ii]<=itheta0[ntod-1]+jj) {upperIndex[jj]=ii; break;}//this yields the smallest ii where itheta0[ii]<=itheta0[ntod-1]+jj
  try {
    TODAsym.alloc(nphi,beammmax+1, NThetaIndex);
  } catch ( bad_alloc & e ) {
    cerr << "itheta0SetUp :  Out of memory allocating " << nphi*(beammmax+1)*NThetaIndex*8./1024/1024 << "MB for TODasym" << endl;
    throw;
  }

  TODAsym.fill(0.);
}

void interpolTOD_arrTestcm_v4(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, long lmax, long lmaxOut, long beammmax, arr<double> &outpntarr1, arr<double> &outpntarr2, long order, arr<double> &TODValue, arr<double> &TODValue2, long ntod, long ntod2, int corenum)
{
  if ( CMULT_VERBOSITY > 1 ) cout << "Entered interpolTOD_arrTestcm_v4"<< endl;
  double phi0=halfpi;
  long npsi=2*beammmax+1;
  long nphi=(2*lmaxOut+1);
  double dphi=2.*pi/nphi;
  double inv_delta_phi = 1./dphi;
  double phioffset = phi0/dphi;
  long ntheta=lmaxOut+1+21;
  double dtheta=-pi/(ntheta-21);
  double theta0=pi-10*dtheta;
  double inv_delta_theta = 1./dtheta;
  long max_order=19;
  int npoints=order+1;
  int ioffset = order/2;
  arr<long> itheta0, itheta0_2;
  long NThetaIndex1(0), NThetaIndex2(0);
  arr<long> lowerIndex, lowerIndex2;
  arr<long> upperIndex, upperIndex2;
  arr3<xcomplex<double> > TODAsym, TODAsym2;
  if (ntod!=0) itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi, beammmax, outpntarr1, ntod, NThetaIndex1, itheta0, lowerIndex, upperIndex, TODAsym);
  if (ntod2!=0)itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi, beammmax, outpntarr2, ntod2, NThetaIndex2, itheta0_2, lowerIndex2, upperIndex2, TODAsym2);

  if ( CMULT_VERBOSITY > 1 ) cout << "DONE with NThetaIndex1 = " << NThetaIndex1 << "  and NThetaIndex2 = " << NThetaIndex2 << "  corenum = " << corenum << "  ntod = " << ntod << "  ntod2 = " << ntod2 << endl;
  bool fNThetaIndex=true;
  if (ntod!=0 && ntod2!=0)
    {
      if (itheta0[ntod-1] != itheta0_2[ntod2-1]) { if ( CMULT_VERBOSITY > 1 ) cout << "itheta0[ntod-1] = " << itheta0[ntod-1] << "   itheta0_2[ntod2-1] = " << itheta0_2[ntod2-1] << "   itheta0[0] = " << itheta0[0] << "   itheta0_2[0] = " << itheta0_2[0] << endl; fNThetaIndex=false;}
      if (NThetaIndex1 != NThetaIndex2) { if ( CMULT_VERBOSITY > 1 ) cout << "NThetaIndex1 = " << NThetaIndex1 << "   NThetaIndex2 = " << NThetaIndex2 << endl; fNThetaIndex=false;}
    }
  if (ntod==0 || ntod2==0) fNThetaIndex=false;

  arr<double> rthetas1, rthetas2;
  if (NThetaIndex1 !=0) {rthetas1.alloc(NThetaIndex1); for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex) rthetas1[thetaIndex]=theta0+(itheta0[ntod-1]+thetaIndex)*dtheta;}
  if (NThetaIndex2 !=0) {rthetas2.alloc(NThetaIndex2); for (int thetaIndex=0; thetaIndex<NThetaIndex2; ++thetaIndex) rthetas2[thetaIndex]=pi - (theta0+(itheta0_2[ntod2-1]+thetaIndex)*dtheta);}

  double elapsed_secs=0.;
  clock_t begin = clock();

  if (fNThetaIndex) 
    conviqt_hemiscm_v4(blmT, slmT, beammmax, lmax, lmaxOut, TODAsym, TODAsym2, NThetaIndex1, rthetas1);
  else
    {
      if (NThetaIndex1!=0) conviqt_hemiscm_single(blmT, slmT, beammmax, lmax, lmaxOut, TODAsym, NThetaIndex1, rthetas1);
      if (NThetaIndex2!=0) conviqt_hemiscm_single(blmT, slmT, beammmax, lmax, lmaxOut, TODAsym2, NThetaIndex2, rthetas2);
    }
  if ( CMULT_VERBOSITY > 1 ) cout << "After conviqt if" << endl;
  clock_t end = clock();
  elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;

  if (ntod>0) for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex) conviqt_tod_loop_v4(lowerIndex, upperIndex, outpntarr1, TODAsym, thetaIndex, itheta0, order, max_order, inv_delta_theta, theta0, inv_delta_phi, phioffset, nphi, ioffset, beammmax, npoints, ntod, TODValue, thetaIndex);
  if (ntod2>0) for (int thetaIndex=0; thetaIndex<NThetaIndex2; ++thetaIndex) conviqt_tod_loop_v4(lowerIndex2, upperIndex2, outpntarr2, TODAsym2, thetaIndex, itheta0_2, order, max_order, inv_delta_theta, theta0, inv_delta_phi, phioffset, nphi, ioffset, beammmax, npoints, ntod2, TODValue2, thetaIndex);

  double maxtod = 0., maxtod2 = 0.;
  if (ntod>0) for (long ii=0; ii<ntod; ii++) if (maxtod<abs(TODValue[ii])) maxtod=abs(TODValue[ii]);
  if (ntod2>0) for (long ii=0; ii<ntod2; ii++) if (maxtod2<abs(TODValue2[ii])) maxtod2=abs(TODValue2[ii]);

  if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "  elapsed_secs for conviqt = " << elapsed_secs << "   average = " << elapsed_secs << "  NThetaIndex1 = " << NThetaIndex1 << "  NThetaIndex2 = " << NThetaIndex2 << "  maxtod = " << maxtod << "  maxtod2 = " << maxtod2 << endl; 
  if ( CMULT_VERBOSITY > 1 ) cout << "Leaving interpolTOD_arrTestcm_v4"<< endl;
}

void interpolTOD_arrTestcm_pol_v4(Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, long lmax, long lmaxOut, long beammmax, arr<double> &outpntarr1, arr<double> &outpntarr2, long order, arr<double> &TODValue, arr<double> &TODValue2, long ntod, long ntod2, int corenum)
{
  if ( CMULT_VERBOSITY > 1 ) cout << "Entered interpolTOD_arrTestcm_pol_v4"<< endl;
  double phi0=halfpi;
  long npsi=2*beammmax+1;
  long nphi=(2*lmaxOut+1);
  double dphi=2.*pi/nphi;
  double inv_delta_phi = 1./dphi;
  double phioffset = phi0/dphi;
  long ntheta=lmaxOut+1+21;
  double dtheta=-pi/(ntheta-21);
  double theta0=pi-10*dtheta;
  double inv_delta_theta = 1./dtheta;
  long max_order=19;
  int npoints=order+1;
  int ioffset = order/2;
  arr<long> itheta0, itheta0_2;
  long NThetaIndex1(0), NThetaIndex2(0);
  arr<long> lowerIndex, lowerIndex2;
  arr<long> upperIndex, upperIndex2;
  arr3<xcomplex<double> > TODAsym, TODAsym2;
  if (ntod!=0) itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi, beammmax, outpntarr1, ntod, NThetaIndex1, itheta0, lowerIndex, upperIndex, TODAsym);
  if (ntod2!=0)itheta0SetUp(npoints, ioffset, ntheta, theta0, inv_delta_theta, nphi, beammmax, outpntarr2, ntod2, NThetaIndex2, itheta0_2, lowerIndex2, upperIndex2, TODAsym2);

  if ( CMULT_VERBOSITY > 1 ) cout << "DONE with NThetaIndex1 = " << NThetaIndex1 << "  and NThetaIndex2 = " << NThetaIndex2 << "  corenum = " << corenum << "  ntod = " << ntod << "  ntod2 = " << ntod2 << endl;
  bool fNThetaIndex=true;
  if (ntod!=0 && ntod2!=0)
    {
      if (itheta0[ntod-1] != itheta0_2[ntod2-1]) { if ( CMULT_VERBOSITY > 1 ) cout << "itheta0[ntod-1] = " << itheta0[ntod-1] << "   itheta0_2[ntod2-1] = " << itheta0_2[ntod2-1] << "   itheta0[0] = " << itheta0[0] << "   itheta0_2[0] = " << itheta0_2[0] << endl; fNThetaIndex=false;}
      if (NThetaIndex1 != NThetaIndex2) { if ( CMULT_VERBOSITY > 1 ) cout << "NThetaIndex1 = " << NThetaIndex1 << "   NThetaIndex2 = " << NThetaIndex2 << endl; fNThetaIndex=false;}
    }
  if (ntod==0 || ntod2==0) fNThetaIndex=false;

  arr<double> rthetas1, rthetas2;
  if (NThetaIndex1 !=0) {rthetas1.alloc(NThetaIndex1); for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex) rthetas1[thetaIndex]=theta0+(itheta0[ntod-1]+thetaIndex)*dtheta;}
  if (NThetaIndex2 !=0) {rthetas2.alloc(NThetaIndex2); for (int thetaIndex=0; thetaIndex<NThetaIndex2; ++thetaIndex) rthetas2[thetaIndex]=pi - (theta0+(itheta0_2[ntod2-1]+thetaIndex)*dtheta);}

  double elapsed_secs=0.;
  clock_t begin = clock();

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "right before conviqt calls: " << std::endl;
      for ( int i=0; i<10; ++i )
	std::cerr << TODAsym(i,0,0) << " " << TODAsym2(i,0,0) << " " << rthetas1[i] << " " << rthetas2[i] << std::endl;
      std::cerr << " fNThetaIndex = " << fNThetaIndex << " blmT(0,0) = " << blmT(0,0) << " slmT(0,0) = " << slmT(0,0) << std::endl;
      std::cerr << std::endl;
    }
  // end DEBUG

  if (fNThetaIndex) 
    conviqt_hemiscm_pol_v4(blmT, slmT, blmG, slmG, blmC, slmC, beammmax, lmax, lmaxOut, TODAsym, TODAsym2, NThetaIndex1, rthetas1);
  else
    {
      if (NThetaIndex1!=0) conviqt_hemiscm_pol_single(blmT, slmT, blmG, slmG, blmC, slmC, beammmax, lmax, lmaxOut, TODAsym, NThetaIndex1, rthetas1);
      if (NThetaIndex2!=0) conviqt_hemiscm_pol_single(blmT, slmT, blmG, slmG, blmC, slmC, beammmax, lmax, lmaxOut, TODAsym2, NThetaIndex2, rthetas2);
    }

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "right after conviqt call (1): " << std::endl;
      for ( int i=0; i<10; ++i )
	std::cerr << TODAsym(i,0,0) << " " << TODAsym2(i,0,0) << " " << rthetas1[i] << " " << rthetas2[i] << std::endl;
      std::cerr << std::endl;
    }
  // end DEBUG

  if ( CMULT_VERBOSITY > 1 ) cout << "After conviqt if" << endl;
  clock_t end = clock();
  elapsed_secs += double(end - begin) / CLOCKS_PER_SEC;

  if (ntod>0) for (int thetaIndex=0; thetaIndex<NThetaIndex1; ++thetaIndex) conviqt_tod_loop_pol_v5(lowerIndex, upperIndex, outpntarr1, TODAsym, thetaIndex, itheta0, order, max_order, inv_delta_theta, theta0, inv_delta_phi, phioffset, nphi, ioffset, beammmax, npoints, ntod, TODValue, thetaIndex);

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "right after conviqt call (2): " << std::endl;
      for ( int i=0; i<10; ++i )
	std::cerr << TODAsym(i,0,0) << " " << TODAsym2(i,0,0) << " " << rthetas1[i] << " " << rthetas2[i] << std::endl;
      std::cerr << std::endl;
    }
  // end DEBUG

if (ntod2>0) for (int thetaIndex=0; thetaIndex<NThetaIndex2; ++thetaIndex) conviqt_tod_loop_pol_v5(lowerIndex2, upperIndex2, outpntarr2, TODAsym2, thetaIndex, itheta0_2, order, max_order, inv_delta_theta, theta0, inv_delta_phi, phioffset, nphi, ioffset, beammmax, npoints, ntod2, TODValue2, thetaIndex);

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "right after conviqt call (3): " << std::endl;
      for ( int i=0; i<10; ++i )
	std::cerr << TODAsym(i,0,0) << " " << TODAsym2(i,0,0) << " " << rthetas1[i] << " " << rthetas2[i] << std::endl;
      std::cerr << std::endl;
    }
  // end DEBUG

  double maxtod = 0., maxtod2 = 0.;
  if (ntod>0) for (long ii=0; ii<ntod; ii++) if (maxtod<abs(TODValue[ii])) maxtod=abs(TODValue[ii]);
  if (ntod2>0) for (long ii=0; ii<ntod2; ii++) if (maxtod2<abs(TODValue2[ii])) maxtod2=abs(TODValue2[ii]);

  if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "  elapsed_secs for conviqt = " << elapsed_secs << "   average = " << elapsed_secs << "  NThetaIndex1 = " << NThetaIndex1 << "  NThetaIndex2 = " << NThetaIndex2 << "  maxtod = " << maxtod << "  maxtod2 = " << maxtod2 << endl; 
  if ( CMULT_VERBOSITY > 1 ) cout << "Leaving interpolTOD_arrTestcm_pol_v4"<< endl;
}

void arrFillingcm_v2(long ntod, arr<double> &timeTest_arr, arr<double> &outpntarrx, arr<double> &outpntarr, long offindex)
{
  timeTest_arr.alloc(ntod);
  outpntarrx.alloc(5*ntod);
  for (long ii=0; ii<ntod; ii++)
    {
      double theta=(outpntarr[5*(ii+offindex)+1]>halfpi) ? pi-outpntarr[5*(ii+offindex)+1]: outpntarr[5*(ii+offindex)+1];
      outpntarrx[5*ii]=outpntarr[5*(ii+offindex)];
      outpntarrx[5*ii+1]=theta;
      outpntarrx[5*ii+2]=outpntarr[5*(ii+offindex)+2];
      outpntarrx[5*ii+3]=outpntarr[5*(ii+offindex)+3];
      outpntarrx[5*ii+4]=outpntarr[5*(ii+offindex)+4];
    }
  hpsort_arrTheta(outpntarrx);
  for (long ii=0; ii<ntod; ii++) timeTest_arr[ii]=outpntarrx[5*ii+4];
}

void todgen_v4(paramfile &par_file, long ntod, long ntod2, arr<double> &todTest_arr, arr<double> &timeTest_arr, arr<double> &todTest_arr2, arr<double> &timeTest_arr2, arr<double> &outpntarr, long lmax, long beammmax, bool pol, Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, int corenum, long dorder)
{
  arr<double> outpntarr1, outpntarr2;
  if (ntod!=0) 
    {
      arrFillingcm_v2(ntod, timeTest_arr, outpntarr1, outpntarr, 0);
      todTest_arr.alloc(ntod);
      todTest_arr.fill(0);
    }
  if (ntod2!=0) 
    {
      arrFillingcm_v2(ntod2, timeTest_arr2, outpntarr2, outpntarr, ntod);
      todTest_arr2.alloc(ntod2);
      todTest_arr2.fill(0);
    }

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "Beginning of outpntarr1&2 and timetest_arr and timetest_arr2: " << std::endl;
      for ( int i=0; i<10; ++i )
	std::cerr << outpntarr1[i] << " " << timeTest_arr[i] << " " << outpntarr2[i] << " " << timeTest_arr2[i] << std::endl;
      std::cerr << std::endl;
    }
  // end DEBUG

  long lmaxOut = par_file.find<long>("lmaxOut", 3000);
  long order= par_file.find<long>("order", 5), offindex;
  order+=dorder;
  if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "  order = " << order << "  ntod = " << ntod << "  ntod2 = " << ntod2 << endl;
  focalplane_db fpdb(par_file);
  if (!pol)
    interpolTOD_arrTestcm_v4(blmT, slmT, lmax, lmaxOut, beammmax, outpntarr1, outpntarr2, order, todTest_arr, todTest_arr2, ntod, ntod2, corenum);
  else
    interpolTOD_arrTestcm_pol_v4(blmT, slmT, blmG, slmG, blmC, slmC, lmax, lmaxOut, beammmax, outpntarr1, outpntarr2, order, todTest_arr, todTest_arr2, ntod, ntod2, corenum);
}

void todgen_fill(paramfile &par_file, long ntod, long ntod2, arr<double> &todTest_arr, arr<double> &timeTest_arr, arr<double> &todTest_arr2, arr<double> &timeTest_arr2, arr<double> &outpntarr, long lmax, long beammmax, bool pol, Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, int corenum, int cores, arr<double> &betaIni_arr, long &NThetaIndex1, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight, arr<int> &iphi0arr2, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &effM, arr<double> &dp1, long dorder)
{
  arr<double> outpntarr1, outpntarr2;
  arrFillingcm_v2(ntod, timeTest_arr, outpntarr1, outpntarr, 0);
  arrFillingcm_v2(ntod2, timeTest_arr2, outpntarr2, outpntarr, ntod);
  todTest_arr.alloc(ntod);
  todTest_arr.fill(0);
  todTest_arr2.alloc(ntod2);
  todTest_arr2.fill(0);
  long lmaxOut = par_file.find<long>("lmaxOut", 3000);
  long order= par_file.find<long>("order", 5), offindex;
  order+=dorder;
  //  if (!pol)
  //    interpolTOD_arrTestcm(blmT, slmT, lmax, lmaxOut, beammmax, outpntarr, order, todTest_arr, ntod, offindex);
  //  else
  if (ntod == 0 || ntod2 == 0) { if ( CMULT_VERBOSITY > 1 ) cout << "ntod = " << ntod << "  ntod2 = " << ntod2 << "   corenum = " << corenum << endl; planck_fail("ntod or ntod2 is zero"); }
    interpolTOD_arrTestcm_pol_fill (blmT, slmT, blmG, slmG, blmC, slmC, lmax, lmaxOut, beammmax, outpntarr1, outpntarr2, order, todTest_arr, todTest_arr2, ntod, ntod2, corenum, cores, betaIni_arr, NThetaIndex1, lowerIndex, upperIndex, lowerIndex2, upperIndex2, iphi0arr, sinweight, cosweight, iphi0arr2, sinweight2, cosweight2, effM, dp1);
}

void todgen_fast(paramfile &par_file, long ntod, long ntod2, arr<double> &todTest_arr, arr<double> &todTest_arr2, long lmax, long beammmax, bool pol, Alm<xcomplex<float> > &blmT, Alm<xcomplex<float> > &slmT, Alm<xcomplex<float> > &blmG, Alm<xcomplex<float> > &slmG, Alm<xcomplex<float> > &blmC, Alm<xcomplex<float> > &slmC, int corenum, arr<double> &betaIni_arr, long &NThetaIndex1, arr<long> &lowerIndex, arr<long> &upperIndex, arr<long> &lowerIndex2, arr<long> &upperIndex2, arr<int> &iphi0arr, arr2<double> &sinweight, arr2<double> &cosweight, arr<int> &iphi0arr2, arr2<double> &sinweight2, arr2<double> &cosweight2, arr<long> &effM, arr<double> &dp1, long dorder, int corechunk, long newcorenum, long coreupper)
{
  todTest_arr.alloc(ntod);
  todTest_arr.fill(0);
  todTest_arr2.alloc(ntod2);
  todTest_arr2.fill(0);
  long lmaxOut = par_file.find<long>("lmaxOut", 3000);
  long order= par_file.find<long>("order", 5), offindex;
  order+=dorder;
  //  if (!pol)
    //    interpolTOD_arrTestcm(blmT, slmT, lmax, lmaxOut, beammmax, outpntarr, order, todTest_arr, ntod, offindex);
  //  else
  if (corenum >= 1*corechunk && corenum < 2*corechunk) readinslm(corenum, newcorenum, coreupper, slmT, slmG, slmC, lmax, pol, par_file);
  interpolTOD_arrTestcm_pol_fast (blmT, slmT, blmG, slmG, blmC, slmC, lmax, lmaxOut, beammmax, order, todTest_arr, todTest_arr2, corenum, betaIni_arr, NThetaIndex1, iphi0arr, iphi0arr2, sinweight, cosweight, sinweight2, cosweight2, lowerIndex, upperIndex, lowerIndex2, upperIndex2, effM, dp1, corechunk, newcorenum, par_file, pol, coreupper);
}

void fillingBetaSeg ( arr < double > & pntarr, long & arrsize, double ratiodeltas, int cores, int corenum, arr<double> &corethetaarr, arr<int> &inBetaSeg)
{
  if ( CMULT_VERBOSITY > 1 ) cout << "Entered fillingBetaSeg in core = " << corenum << endl;
  for (long ii=0; ii<arrsize; ii++) 
    {
      double theta=pntarr[5*ii+1];
      if (theta >= halfpi) theta=pi-theta;
      long ratiobeta;
      ratiobeta=theta/ratiodeltas;
      if (theta>=corethetaarr[ratiobeta])
        ratiobetacalcgreatercm(ratiobeta, theta, cores, corethetaarr);
      else
        ratiobetacalcsmallercm(ratiobeta, theta, cores, corethetaarr);
      inBetaSeg[ratiobeta]+=5;
    }
  if ( CMULT_VERBOSITY > 1 ) cout << "Leaving fillingBetaSeg in core = " << corenum << endl;
}

void deltaTheta2(long iival, double thetaini, arr<double> &dbeta)
{
  dbeta.alloc(iival+1);
  double dy=(1.-thetaini)/2.;
  //  double theta2=pi-abs(asin(1-dy));
  double yshift=-0.7;
  double theta2=pi-abs(asin(1-dy))+yshift;
  //  double theta1=pi-abs(asin(dy));
  double theta1=pi-abs(asin(sin(theta2)-thetaini));
  double dt=(theta2-theta1)/(iival+1.);
  double dbeta0=sin(theta2);
  if ( CMULT_VERBOSITY > 1 ) cout << "dy = " << dy << "   theta1 = " << theta1 << "   theta2 = " << theta2 << "   dbeta0 = " << dbeta0 << "  dt = " << dt << "   iival = " << iival <<"  sin(theta2-(iival+1-iival)*dt) = " << sin(theta2-(iival+1-iival)*dt) <<"  sin(theta2-(iival+1-(iival-1) )*dt) = " << sin(theta2-(iival+1-(iival-1) )*dt) << endl;
  for (long ii=iival; ii>=0; ii--)
    {
      dbeta[ii]=dbeta0-sin(theta2-(iival+1-ii)*dt);
      dbeta0=sin(theta2-(iival+1-ii)*dt);
      //      if (ii>=iival-2) cout << "ii = " << ii << "  dbeta0 = " << dbeta0 << "  dbeta[" << ii <<  "] = " << dbeta[ii] << endl;
    }
}

void preReorderingStep(long ntod, long ntod2, arr<double> &todAll, arr<double> &todTest_arr, arr<double> &todTest_arr2)
{
  if ((ntod+ntod2)!=0) todAll.alloc(ntod+ntod2);
  for (long ii=0; ii<ntod; ii++) todAll[ii]=todTest_arr[ii];
  if (ntod!=0) todTest_arr.dealloc();
  for (long ii=0; ii<ntod2; ii++) todAll[ntod+ii]=todTest_arr2[ii];
  if (ntod2!=0) todTest_arr2.dealloc();
}

void tod_calc4 (paramfile &par_file, arr<double> &pntarr, double &t_read)
{
  long totsize=pntarr.size()/5;
  int cores=mpiMgr.num_ranks(), corenum=mpiMgr.rank();
  bool pol = par_file.find<bool>("polarisation",false);
  long lmax = par_file.find<long>("lmax", 5000);
  long beammmax = par_file.find<long>("beammmax", 14);
  long Nbetafac = par_file.find<long>("Nbetafac", 2000);
  long Nbeta=cores*Nbetafac;
  long MCSamples = par_file.find<long>("MCSamples", 0);

  long ntod=0, ntod2=0, outBetaSegSize;
  const string det_id = par_file.find<string>("detector_id");
  bool galactic = par_file.find<bool>("galactic", true);
  arr<int> inBetaSeg, inBetaSegFir, inBetaSegSec, outBetaSeg;
  long arrsize, countFirstHalf, countSecondHalf;
  inBetaSeg.alloc(cores);
  inBetaSeg.fill(0);
  inBetaSegFir.alloc(cores);
  inBetaSegFir.fill(0);
  inBetaSegSec.alloc(cores);
  inBetaSegSec.fill(0);
  arr<int> inBetaSegAcc, outBetaSegAcc;
  arr<double> pntarr2, outpntarr;
  arr<int> inOffset, outOffset;
  double deltabeta=pi/Nbeta;
  double ratiodeltas=halfpi/cores;
  double dtheta, newtheta, thetaini=halfpi;
  double countdeltatheta=0.;
  double coretheta, coredtheta;
  arr<double> corethetaarr(cores), coredthetaarr(cores);
  arr<double> dbeta;
  for (long ii=cores-1; ii>=0; ii--)
    {
      if (thetaini>0.18) 
	thetaDeltaThetacm(thetaini, cores, ii, newtheta, dtheta);
      else
	{
	  //	  dtheta=newtheta/(ii+1.);
	  long iival=ii;
	  deltaTheta2(iival, thetaini, dbeta);
	  for (long jj=ii; jj>=0; jj--) 
	    {
	      //	      newtheta-=dtheta;
	      newtheta-=dbeta[jj];
	      corethetaarr[jj]=newtheta;
	      coredthetaarr[jj]=dbeta[jj];
	      thetaini=newtheta;
	    }
	  newtheta=0.;
	  corethetaarr[0]=0.;
	  break;
	}
      //      thetaDeltaThetacm(thetaini, cores, ii, newtheta, dtheta);
      if (corenum==ii) {coretheta=newtheta; coredtheta=dtheta;}
      corethetaarr[ii]=newtheta;
      coredthetaarr[ii]=dtheta;
      thetaini=newtheta;
      countdeltatheta+=sin(newtheta+dtheta/2.)*dtheta;
    }
  fillingBetaSeg ( pntarr, totsize, ratiodeltas, cores, corenum, corethetaarr, inBetaSeg);
  mpiMgr.all2all(inBetaSeg, outBetaSeg);

  todRedistribution5cm (pntarr, inBetaSeg, outBetaSeg, inBetaSegAcc, outBetaSegAcc, outBetaSegSize, cores, pntarr2, totsize, inOffset, outOffset, ratiodeltas, corethetaarr);
  arr<double> todtmp;
  if (totsize>0)
    {
      todtmp.alloc(totsize);
      for (long ii=0; ii<totsize; ii++) todtmp[ii]=pntarr[5*ii+3];
      pntarr.dealloc();
    }
  if (totsize !=0 || outBetaSegSize != 0) mpiMgr.all2allv(pntarr2, inBetaSeg, inOffset, outpntarr, outBetaSeg, outOffset, outBetaSegSize);
  if (totsize>0) pntarr2.dealloc();
  
  long nthetaCount(0);
  long nthetaCount2(0);
  long csmall(0), cbig(0);
  for (long ii=0; ii<outBetaSegSize/5; ii++)
    {
      if (outpntarr[5*ii+1] < halfpi)
	{
	  nthetaCount++; 
	  if (outpntarr[5*ii+1] <= 0) csmall++; 
	}
      else 
	{
	  nthetaCount2++;
	  if (outpntarr[5*ii+1] >= pi) cbig++; 
	}
      outpntarr[5*ii+3]=ii*1.;
    }
  ntod=nthetaCount; 
  ntod2=nthetaCount2;
  if (csmall!=0 || cbig!=0) if ( CMULT_VERBOSITY > 1 ) cout << "WARNING: PRODUCTION OF TOD AT THETA=0 OR PI IN corenum = " << corenum << ".   # OF TOD AT THETA=0: " << csmall << "  # OF TOD AT THETA=PI: " << cbig << endl;
  arr<double> todTest_arr, todTest_arr2, timeTest_arr, timeTest_arr2;
  if (outBetaSegSize!=0) hpsort_arrTheta(outpntarr);
  long dorder=0;
  //  if (outpntarr[1]<0.05) dorder=2;

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "Beginning of outpntarr: ";
      for ( int i=0; i<10; ++i )
	std::cerr << outpntarr[i] << ", ";
      std::cerr << std::endl;
    }
  // end DEBUG

  long ncores=125;
  long coreratio=corenum/ncores;
  long newcorenum=ncores*coreratio;
  long coreupper=(newcorenum+ncores>=cores) ? cores : newcorenum+ncores;

  Alm<xcomplex<float> > blmT(lmax,beammmax), blmG(lmax,beammmax), blmC(lmax,beammmax);

  double t = MPI_Wtime();
  readinblm(corenum, newcorenum, coreupper, blmT, blmG, blmC, lmax, beammmax, pol, par_file);
  t_read += MPI_Wtime() - t;

  arr<int> iphi0arr, iphi0arr2;
  arr2<double> sinweight, cosweight, sinweight2, cosweight2;
  arr<long> lowerIndex, upperIndex, lowerIndex2, upperIndex2;
  arr<double> betaIni_arr;
  long NThetaIndex1;
  arr<long> effM(lmax+1);
  arr<double> dp1;
  if (MCSamples==0)
    {
      clock_t begin = clock();

      Alm<xcomplex<float> > slmInputT(lmax,lmax), slmInputG(lmax,lmax), slmInputC(lmax,lmax);
      double t = MPI_Wtime();
      readinslm(corenum, newcorenum, coreupper, slmInputT, slmInputG, slmInputC, lmax, pol, par_file);
      t_read += MPI_Wtime() - t;
      if (ntod!=0 || ntod2!=0) 
	todgen_v4(par_file, ntod, ntod2, todTest_arr, timeTest_arr, todTest_arr2, timeTest_arr2, outpntarr, lmax, beammmax, pol, blmT, slmInputT, blmG, slmInputG, blmC, slmInputC, corenum, dorder);
      else
	if ( CMULT_VERBOSITY > 1 ) cout << "no data to process in corenum = " << corenum << endl;
      // DEBUG
      if ( corenum == 0 )
	{
	  std::cerr << "Beginning of outpntarr (v2): ";
	  for ( int i=0; i<10; ++i )
	    std::cerr << outpntarr[i] << ", ";
	  std::cerr << std::endl;
	}
      if ( corenum == 0 )
	{
	  std::cerr << "Beginning of todtest and timetest: " << std::endl;
	  for ( int i=0; i<10; ++i )
	    std::cerr << todTest_arr[i] << " " << timeTest_arr[i] << " " << todTest_arr2[i] << " " << timeTest_arr2[i] << std::endl;
	  std::cerr << std::endl;
	}
      // end DEBUG

      //    outpntarr.dealloc();
      clock_t end = clock();
      double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
      
      arr<double> tarr(cores,elapsed_secs), toutarr;
      arr<long> arrcores(cores); for (long ii=0; ii<cores; ii++) arrcores[ii]=ii;
      arr<long> outbetasegsizearr(cores,outBetaSegSize), obt;
      mpiMgr.all2all(tarr,toutarr);
      mpiMgr.all2all(outbetasegsizearr,obt);
      double maxtint=0;
      double mintint=1.e30;
      long corevaluemax, corevaluemin;
      for (long ii=0; ii<cores; ii++)
	{
	  if (toutarr[ii] > maxtint) { maxtint=toutarr[ii]; corevaluemax=ii;}
	  if (toutarr[ii] < mintint) { mintint=toutarr[ii]; corevaluemin=ii;}
	}
      if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "   elapsed_secs = " << elapsed_secs << endl;
      if (corenum==0) if ( CMULT_VERBOSITY > 1 ) cout << "maxtint = " << maxtint << "   mintint = " << mintint << "   corevaluemax = " << corevaluemax << "   corevaluemin = " << corevaluemin << endl;
      hpsort_DL(toutarr, arrcores);
      if (corenum==0) for (long ii=0; ii<cores; ii++) if ( CMULT_VERBOSITY > 1 ) cout << "toutarr = " << toutarr[ii] <<"   corenum = " << arrcores[ii] <<"   outBetaSegSize = " << obt[arrcores[ii]] << endl;
    }
  else
    {
      {
	clock_t begin = clock();
	
	Alm<xcomplex<float> > slmInputT(lmax,lmax), slmInputG(lmax,lmax), slmInputC(lmax,lmax);
	readinslm(corenum, newcorenum, coreupper, slmInputT, slmInputG, slmInputC, lmax, pol, par_file);
	todgen_fill(par_file, ntod, ntod2, todTest_arr, timeTest_arr, todTest_arr2, timeTest_arr2, outpntarr, lmax, beammmax, pol, blmT, slmInputT, blmG, slmInputG, blmC, slmInputC, corenum, cores, betaIni_arr, NThetaIndex1, lowerIndex, upperIndex, lowerIndex2, upperIndex2, iphi0arr, sinweight, cosweight, iphi0arr2, sinweight2, cosweight2, effM, dp1, dorder);
	clock_t end = clock();
	//    outpntarr.dealloc();
	double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
	
	arr<double> tarr(cores,elapsed_secs), toutarr;
	arr<long> arrcores(cores); for (long ii=0; ii<cores; ii++) arrcores[ii]=ii;
	arr<long> outbetasegsizearr(cores,outBetaSegSize/5), obt;
	mpiMgr.all2all(tarr,toutarr);
	mpiMgr.all2all(outbetasegsizearr,obt);
	double maxtint=0;
	double mintint=1.e30;
	long corevaluemax, corevaluemin;
	for (long ii=0; ii<cores; ii++)
	  {
	    if (toutarr[ii] > maxtint) { maxtint=toutarr[ii]; corevaluemax=ii;}
	    if (toutarr[ii] < mintint) { mintint=toutarr[ii]; corevaluemin=ii;}
	  }
	if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "   elapsed_secs = " << elapsed_secs << endl;
	if (corenum==0) if ( CMULT_VERBOSITY > 1 ) cout << "maxtint = " << maxtint << "   mintint = " << mintint << "   corevaluemax = " << corevaluemax << "   corevaluemin = " << corevaluemin << endl;
	hpsort_DL(toutarr, arrcores);
	if (corenum==0) for (long ii=0; ii<cores; ii++) if ( CMULT_VERBOSITY > 1 ) cout << "toutarr = " << toutarr[ii] <<"   corenum = " << arrcores[ii] <<"   outBetaSegSize = " << obt[arrcores[ii]] << endl;
	todTest_arr.dealloc(); todTest_arr2.dealloc();
      }
      int corechunk=cores/3;
      clock_t begin = clock();
      for (long mc=0; mc<MCSamples; mc++)
	{

	  Alm<xcomplex<float> > slmInputT(lmax,lmax), slmInputG(lmax,lmax), slmInputC(lmax,lmax);
	  if (corenum >= 0*corechunk && corenum < 1*corechunk) readinslm(corenum, newcorenum, coreupper, slmInputT, slmInputG, slmInputC, lmax, pol, par_file);
	  
	  todgen_fast(par_file, ntod, ntod2, todTest_arr, todTest_arr2, lmax, beammmax, pol, blmT, slmInputT, blmG, slmInputG, blmC, slmInputC, corenum, betaIni_arr, NThetaIndex1, lowerIndex, upperIndex, lowerIndex2, upperIndex2, iphi0arr, sinweight, cosweight, iphi0arr2, sinweight2, cosweight2, effM, dp1, dorder, corechunk, newcorenum, coreupper);
	}
      iphi0arr.dealloc(); iphi0arr2.dealloc();
      sinweight.dealloc(); cosweight.dealloc(); sinweight2.dealloc(); cosweight2.dealloc();
      lowerIndex.dealloc(); upperIndex.dealloc(); lowerIndex2.dealloc(); upperIndex2.dealloc();
      betaIni_arr.dealloc();

      clock_t end = clock();
      double elapsed_secs2 = double(end - begin) / CLOCKS_PER_SEC;
      arr<double> tarr(cores,elapsed_secs2), toutarr;
      arr<long> arrcores(cores); for (long ii=0; ii<cores; ii++) arrcores[ii]=ii;
      arr<long> outbetasegsizearr(cores,outBetaSegSize/5), obt;
      mpiMgr.all2all(tarr,toutarr);
      mpiMgr.all2all(outbetasegsizearr,obt);
      double maxtint=0;
      double mintint=1.e30;
      long corevaluemax, corevaluemin;
      for (long ii=0; ii<cores; ii++)
	{
	  if (toutarr[ii] > maxtint) { maxtint=toutarr[ii]; corevaluemax=ii;}
	  if (toutarr[ii] < mintint) { mintint=toutarr[ii]; corevaluemin=ii;}
	}
      if ( CMULT_VERBOSITY > 1 ) cout << "corenum = " << corenum << "   elapsed_secs2 = " << elapsed_secs2 << endl;
      if (corenum==0) if ( CMULT_VERBOSITY > 1 ) cout << "maxtint2 = " << maxtint << "   mintint2 = " << mintint << "   corevaluemax2 = " << corevaluemax << "   corevaluemin2 = " << corevaluemin << endl;
      hpsort_DL(toutarr, arrcores);
      if (corenum==0) for (long ii=0; ii<cores; ii++) if ( CMULT_VERBOSITY > 1 ) cout << "toutarr2 = " << toutarr[ii] <<"   corenum = " << arrcores[ii] <<"   outBetaSegSize = " << obt[arrcores[ii]] << endl;
    }

  if (outBetaSegSize!=0) hpsort_arrTOD(outpntarr);
  if (totsize !=0 || outBetaSegSize != 0) mpiMgr.all2allv(outpntarr, outBetaSeg, outOffset, pntarr, inBetaSeg, inOffset, 5*totsize);

  if (outBetaSegSize!=0) outpntarr.dealloc();
  if (totsize!=0) hpsort_arrTime(pntarr);

  arr<double> todAll;
  preReorderingStep(ntod, ntod2, todAll, todTest_arr, todTest_arr2);
  arr<double> timeAll;
  preReorderingStep(ntod, ntod2, timeAll, timeTest_arr, timeTest_arr2);

  if (outBetaSegSize!=0) hpsort_DDcm(timeAll,todAll);
  for (long ii=0; ii<cores; ii++)
    {
      inBetaSeg[ii]=inBetaSeg[ii]/5;
      outBetaSeg[ii]=outBetaSeg[ii]/5;
      inOffset[ii]=inOffset[ii]/5;
      outOffset[ii]=outOffset[ii]/5;
    }
  arr<double> outtodarr, outnumarr;
  if (totsize !=0 || outBetaSegSize != 0) mpiMgr.all2allv(todAll, outBetaSeg, outOffset, outtodarr, inBetaSeg, inOffset, totsize);
  if (totsize !=0 || outBetaSegSize != 0) mpiMgr.all2allv(timeAll, outBetaSeg, outOffset, outnumarr, inBetaSeg, inOffset, totsize);

  // DEBUG
  if ( corenum == 0 )
    {
      std::cerr << "Beginning of todAll and timeAll: " << std::endl;
      for ( int i=0; i<10; ++i )
	std::cerr << todAll[i] << " " << timeAll[i] << " " << outtodarr[i] << " " << outnumarr[i] << std::endl;
      std::cerr << std::endl;
    }
  // end DEBUG
  
  if (outBetaSegSize!=0) todAll.dealloc();
  if (outBetaSegSize!=0) timeAll.dealloc();

  double maxtodall(0), mintodall(1e10);
  for (long ii=0; ii< totsize; ii++)
    {
      maxtodall=(maxtodall < outtodarr[ii]) ? outtodarr[ii] : maxtodall;
      mintodall=(mintodall > outtodarr[ii]) ? outtodarr[ii] : mintodall;
    }
  if (totsize > 0) if ( CMULT_VERBOSITY > 1 ) cout << "maxtodall = " << maxtodall  << "   mintodall = " << mintodall << "   outtodarr.size()/totsize = " << outtodarr.size()/totsize*1. << endl;

  focalplane_db fpdb(par_file);
  double calibration=1.;
  if (par_file.find<bool> ("calibrate_signal")) calibration=2./(1.+fpdb.getValue<double>(det_id,"epsilon"));
  if (totsize > 0)
    {
      hpsort_DDcm(outnumarr,outtodarr);
      outnumarr.dealloc();
      double maxtoddiff=0.;
      for (long ii=0; ii<totsize; ii++) 
	{
	  pntarr[5*ii+3]=calibration*outtodarr[ii];
	  maxtoddiff=(abs(todtmp[ii]-pntarr[5*ii+3]) > maxtoddiff) ? abs(todtmp[ii]-pntarr[5*ii+3]) : maxtoddiff;
	  if (ii%100000 == 0) if ( CMULT_VERBOSITY > 1 ) cout << "todtmp[ii] = " << todtmp[ii] << "  pntarr[5*ii+3] = " << pntarr[5*ii+3] << "  difference = " << abs(todtmp[ii]-pntarr[5*ii+3]) << endl;
	}
      if ( CMULT_VERBOSITY > 1 ) cout << "  corenum = " << corenum << "   maxtoddiff = " << maxtoddiff << "   calibration = " << calibration << endl;
    }
}



int cmult_module (int argc, const char **argv)
  {
  module_startup ("cmult", argc, const_cast<const char **>(argv), 2,
    "<init object>", mpiMgr.master());
  iohandle_current::Manager mng (argv[1]);
  paramfile par_file (mng.getParams(mpiMgr.master()));

  return 0;
  }
