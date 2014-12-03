#include <stdexcept>
#include <complex>
#include <algorithm>
//#include <toast_mpi.hpp>
//#include <toast.hpp>
#include <cmath>
#include <cstdlib>
#include <time.h>
#include "planck_rng.h"
#include "lsconstants.h"
#include "arr.h"
#include "cxxutils.h"
#include "focalplane_db.h"
#include "paramfile.h"
#include "mpi_support.h"
#include "iohandle_current.h"
#include "repointing.h"
#include "levels_facilities.h"
#include "xcomplex.h"
#include "trafos.h"
#include "walltimer.h"
#include "focalplane_db.h"

#include <healpix_map.h>
#include <healpix_map_fitsio.h>
#include <healpix_base.h>
#include <healpix_tables.h>
#include <pointing.h>

#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>

#include "injector.hpp"

using namespace std;

int bin_map( paramfile & par_file, arr<double> pntarr, int id ) {
  string file_binmap = par_file.find<string>("binmap_file","");
  string file_hitmap = par_file.find<string>("hitmap_file","");

  if ( file_binmap.size() == 0 && file_hitmap.size() == 0 ) return 0;

  int nside = par_file.find<int>("binmap_nside",64);
  long npix = 12*nside*nside;

  Healpix_Map<double> binmap;
  binmap.SetNside( nside, RING );
  binmap.fill( 0.0 );

  Healpix_Map<int> hitmap;
  hitmap.SetNside( nside, RING );
  hitmap.fill( 0 );

  //T_Healpix_Base<long> hpx;
  //hpx.SetNside( nside, RING );

  for ( int i = 0; i < pntarr.size(); i += 5 ) {
    double phi = pntarr[i];
    double theta = pntarr[i+1];
    pointing ang( theta, phi );
    long pix = binmap.ang2pix( ang );

    binmap[pix] += pntarr[i+3];
    hitmap[pix]++;
  }

  if ( id != 0 ) {
    MPI_Reduce( &(hitmap[0]), NULL, npix, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &(binmap[0]), NULL, npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  } else {
    MPI_Reduce( MPI_IN_PLACE, &(hitmap[0]), npix, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &(binmap[0]), npix, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

    for ( int i = 0; i < npix; ++i ) if ( hitmap[i] != 0 ) binmap[i] /= hitmap[i];

    int err;
    if ( file_hitmap.size() != 0 ) {
      err = remove( file_hitmap.c_str() );
      write_Healpix_map_to_fits( file_hitmap, hitmap, PLANCK_INT32 );
    }
    if ( file_binmap.size() != 0 ) {
      err = remove( file_binmap.c_str() );
      write_Healpix_map_to_fits( file_binmap, binmap, PLANCK_FLOAT64 );
    }

    cout << "Wrote the binned map to " << file_binmap << endl;
  }

  return 0;
}

// This functor sums up the number of samples assigned to this process
// FOR EACH CHANNEL.
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
  
  toast_chunk_process ( bool galactic, string const & det_id, bool apply_flags, arr < double > * pntarr, double psi_pol ) 
  {
    det_ = det_id;
    gtype_ = toast::GEOM_THETAPHI;
    
    if ( galactic ) 
      {
        coord_ = toast::PNTG_GALACTIC;
      } 
    else 
      {
        coord_ = toast::PNTG_ECLIPTIC;
      }

    apply_flags_ = apply_flags;
    psi_pol_ = psi_pol;
    
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
  bool apply_flags_;
  double psi_pol_;
  
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
      if ( (*itchan)->name().find(det_) != string::npos ) {
        // this is the channel we want!
        std::vector < toast::channel_p > chan(1);
        chan[0] = (*itchan);
        double psipol = (*itchan)->query_detector()->response_angle();
        
        // read channel pointing
        
        tele->channel_geom_pointing ( chnk->observation(), chnk->first(), chnk->last(), chan, gtype_, coord_, chanpnt, pntflags );
        
        // read channel data and flags          
        
        (*itchan)->read ( chnk->observation(), "DEFAULT", chnk->first(), chnk->last(), chandata, flags );
        
        // read timestamps
        
        chnk->observation()->global_times ( chnk->first(), len, timestamps );
        
        int64_t chunkoff;
        long outoff;
        
        for ( int64_t i = 0; i < len; ++i ) {
          if ( ( flags[i] | pntflags[i] ) == 0 || !apply_flags_ ) {
            // no flags, data is good
            chunkoff = 3 * i;
            outoff = 5 * good_;
            
            (*out_)[ outoff ] = (chanpnt[0])( chunkoff + 1 );
            (*out_)[ outoff + 1 ] = (chanpnt[0])( chunkoff );
            (*out_)[ outoff + 2 ] = (chanpnt[0])( chunkoff + 2 ) - psi_pol_;
            (*out_)[ outoff + 3 ] = chandata(i);
            (*out_)[ outoff + 4 ] = timestamps(i);
            ++good_;
          }
        }
        
      }
    }
    //cout << "good_ = " << good_ << "  len = " << len << endl; 
  }

};
*/
void ratiobetacalcgreatercm(long &ratiobeta, double theta, int cores, arr<double> &corethetaarr);

void ratiobetacalcsmallercm(long &ratiobeta, double theta, int cores, arr<double> &corethetaarr);

void thetaDeltaThetacm(double thetaini, int cores, int corenum, double &theta, double &deltatheta);
/*
void toastReadData ( arr < double > & pntarr, const string det_id, paramfile & par_file, bool galactic, bool apply_flags, double ratiodeltas, int cores, arr<double> &corethetaarr, int corenum, double psi_pol, string toast_distribution )
{
  // Read run file
  string toast_run_file=par_file.find<string>("toast_run_file");

  //  toast::mpi_run conf ( comm );
  toast::mpi_run conf ( MPI_COMM_WORLD );
  conf.mpi_read_or_load ( 0, toast_run_file );

  // Distribute data

  toast::dist_type dtype;
  if ( toast_distribution.find("obs") != string::npos || toast_distribution.find("OBS") != string::npos ) {
    dtype = toast::DIST_OBS;
  } else if ( toast_distribution.find("int") != string::npos || toast_distribution.find("INT") != string::npos ) {
    dtype = toast::DIST_INTERVAL;
  } else {
    throw runtime_error( string("Unknown TOAST distribution: ") + toast_distribution );
  }

  toast::mpi_distribution_p dist ( new toast::mpi_distribution ( conf, dtype, 1 ) );

  // Compute number of samples

  toast_samples samples;

  dist->exec ( samples );

  // allocate output array

  long totsize = samples.total_;
  long totnchunk = samples.nchunk_;
  long totsize_global=0;
  long totnchunk_global=0;
  MPI_Reduce( &totsize, &totsize_global, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );
  MPI_Reduce( &totnchunk, &totnchunk_global, 1, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD );

  if ( corenum == 0 ) {
    if ( totsize_global * 40. <  pow(2.,30.) ) 
      cout << "Total number of samples: " << totsize_global << " ( " << totsize_global * 40. / pow(2.,20.) << " MB ) in " << totnchunk_global << " chunks." << endl;
    else
      cout << "Total number of samples: " << totsize_global << " ( " << totsize_global * 40. / pow(2.,30.) << " GB ) in " << totnchunk_global << " chunks." << endl;
  }

  if ( totsize == 0 ) {
    cout << corenum << " : No samples to read from TOAST." << endl;
    return;
  }

  //cout << "internal totsize = "<< totsize << "  in MB = " << 5*totsize*8./1000000 << "MB in corenum = " << corenum << endl;

  arr<double> pntarrtmp( 5 * totsize );

  // Create functor for iterating over all local data chunks

  toast_chunk_process tproc ( galactic, det_id, apply_flags, &pntarrtmp, psi_pol );

  // Process data

  dist->exec ( tproc );

  long arrsize = tproc.good_;
  pntarr.alloc ( 5 * arrsize );

  double oldtheta=pntarrtmp[5*(arrsize-1)+1];
  double oldphi=pntarrtmp[5*(arrsize-1)];
  for (long ii=0; ii<arrsize; ii++) 
    {
      double theta=pntarrtmp[5*ii+1];
      if (theta >= halfpi) theta=pi-theta;
      long ratiobeta;
      ratiobeta=theta/ratiodeltas;
      if (theta>=corethetaarr[ratiobeta])
	ratiobetacalcgreatercm(ratiobeta, theta, cores, corethetaarr);
      else
	ratiobetacalcsmallercm(ratiobeta, theta, cores, corethetaarr);
      pntarr[5*ii]=pntarrtmp[5*ii];
      pntarr[5*ii+1]=pntarrtmp[5*ii+1];
      pntarr[5*ii+2]=pntarrtmp[5*ii+2];
      pntarr[5*ii+3]=pntarrtmp[5*ii+3];
      pntarr[5*ii+4]=pntarrtmp[5*ii+4];
    }
}
*/
void deltaTheta2(long iival, double thetaini, arr<double> &dbeta);

void tod_calc4 (paramfile &par_file, arr<double> &pntarr, double &t_read);

void IScalm2TODEx (paramfile &par_file)
{
  double t_init, t_read_tod, t_read_alm=0, t_convolve, t_bin, t_init_write, t_write=0, t0, t1, t2, t_tot;

  t0 = t1 = MPI_Wtime();

  int cores=mpiMgr.num_ranks(), corenum=mpiMgr.rank();

  // Parse the convolution parameters

  const string det_id = par_file.find<string>("detector_id");
  string outfile = par_file.find<string>("outSky");
  bool galactic = par_file.find<bool>("galactic", true);
  bool apply_flags = par_file.find<bool>("apply_flags", false);
  string beam_frame = par_file.find<string>("beam_frame"); // No default, force user to choose
  double ratiodeltas=halfpi/cores;
  double dtheta, newtheta, thetaini=halfpi;
  arr<double> corethetaarr(cores);
  arr<double> dbeta;

  // Parse the injection parameters

  string eff_dir = ""; // par_file.find<string>("eff_dir");
  string file_pattern = ""; // par_file.find<string>("file_pattern",".*");
  string det_pattern = ""; // par_file.find<string>("det_pattern",".*"+det_id+".*");
  bool add = false; // par_file.find<bool>("add",false);
  double scale;
  try {
    scale = 1; // par_file.find<double>("scale",1.0);
  } catch (...) {
    if ( corenum == 0 ) cerr << "ERROR: failed to parse the scale factor" << endl;
    throw;
  }

  // Determine the beam frame and possible polarization angle correction

  double psi_pol;
  focalplane_db fpdb(par_file);
  if ( boost::iequals(beam_frame, "Pxx") ) {
    // Beam is defined with respect to the polarization basis of the horn:
    // TOAST psi angle includes psi_uv+psi_pol -> subtract psi_pol according to FPDB but leave psi_uv
    psi_pol = fpdb.getValue<double>(det_id, "psi_pol") * degr2rad;
  } else if ( boost::iequals(beam_frame, "Dxx") ) {
    // Beam is defined in the focal plane frame. Must subtract both psi_pol and psi_uv from TOAST pointing
    psi_pol = (fpdb.getValue<double>(det_id, "psi_uv") + fpdb.getValue<double>(det_id, "psi_pol") ) * degr2rad;
  } else
    throw runtime_error( string("Unknown beam frame: ") + beam_frame );

    // Set up the convolution
    
  for (long ii=cores-1; ii>=0; ii--)
    {
      thetaDeltaThetacm(thetaini, cores, ii, newtheta, dtheta);
      corethetaarr[ii]=newtheta;
      thetaini=newtheta;
    }

  for (long ii=cores-1; ii>=0; ii--)
    {
      if (thetaini>0.18) 
	thetaDeltaThetacm(thetaini, cores, ii, newtheta, dtheta);
      else
	{
	  long iival=ii;
	  deltaTheta2(iival, thetaini, dbeta);
	  for (long jj=ii; jj>=0; jj--) 
	    {
	      newtheta-=dbeta[jj];
	      corethetaarr[jj]=newtheta;
	      thetaini=newtheta;
	    }
	  newtheta=0.;
	  corethetaarr[0]=0.;
	  break;
	}
      corethetaarr[ii]=newtheta;
      thetaini=newtheta;
    }

  t_init = (t2=MPI_Wtime()) - t1;

  arr<double> pntarr;
  
  //toastReadData ( pntarr, det_id, par_file, galactic, apply_flags, ratiodeltas, cores, corethetaarr, corenum, psi_pol, par_file.find<string>("toast_distribution","OBSERVATION") );

  // Populate the pointing array

  long ntheta = 3;
  long nphi = 3;
  long npsi = 3;
  long nsamp = ntheta*nphi*npsi;

  pntarr.alloc( 5 * nsamp );

  long row=0;
  for ( long itheta=0; itheta < ntheta; ++itheta ) {
    double theta = itheta * (pi / (double)ntheta);
    for ( long iphi=0; iphi < nphi; ++iphi ) {
      double phi = iphi * (twopi / (double)nphi);
      for ( long ipsi=0; ipsi < npsi; ++ipsi ) {
	double psi = ipsi * (pi / (double)npsi);
	pntarr[row*5+0] = phi; // longitude
	pntarr[row*5+1] = theta; // latitude
	pntarr[row*5+2] = psi; // position angle
	pntarr[row*5+3] = 0; // TOD
	pntarr[row*5+4] = row; // time
	++row;
      }
    }
  }

  t_read_tod = (t1=MPI_Wtime()) - t2;

  MPI_Barrier( MPI_COMM_WORLD );

  long ngood = pntarr.size() / 5;
  MPI_Allreduce( MPI_IN_PLACE, &ngood, 1, MPI_LONG, MPI_SUM, MPI_COMM_WORLD );
  if ( ngood == 0 ) {
    if ( corenum == 0 ) cout << "WARNING: There were no unflagged samples to simulate for " << det_id << endl;
    return;
  }

  // Run the first (possibly only) beam component

  if ( par_file.param_present( string("infile_beam1") ) ) {
    par_file.setParam<string>( "infile_beam", par_file.find<string>("infile_beam1") );
    par_file.setParam<long>( "lmax", par_file.find<long>("lmax1") );
    par_file.setParam<long>( "beammmax", par_file.find<long>("beammmax1") );
  }

  if ( corenum == 0 ) cout << "Convolving beam # " << 1 << " = " << par_file.find<string>("infile_beam") << endl;

  try {
    tod_calc4 (par_file, pntarr, t_read_alm);
  } catch ( bad_alloc & e ) {
    cerr << corenum << " : Out of memory calling tod_calc4 for beam # 1 " << endl;
    throw;
  }

  if ( par_file.param_present( string("infile_beam2") ) ) {

    // Run additional beam components

    arr<double> todtmp;
    long totsize = (long) pntarr.size() / (long) 5;
    try {
      todtmp.alloc(totsize);
    } catch ( bad_alloc & e ) {
      cerr << corenum << " : Out of memory allocating " << totsize*8./1024/1024 << "MB for temporary TOD" << endl;
      throw;
    }

    for (long ii=0; ii<totsize; ii++) todtmp[ii] = pntarr[5*ii+3];

    for ( int ibeam=2; ibeam<100; ++ibeam ) {
      string sbeam( to_string(ibeam) );
      if ( par_file.param_present( string("infile_beam")+sbeam ) ) {
        par_file.setParam<string>( "infile_beam", par_file.find<string>(string("infile_beam")+sbeam) );
        par_file.setParam<long>( "lmax", par_file.find<long>(string("lmax")+sbeam) );
        par_file.setParam<long>( "beammmax", par_file.find<long>(string("beammmax")+sbeam) );
      } else
        break;

      if ( corenum == 0 ) cout << "Convolving beam # " << ibeam << " = " << par_file.find<string>("infile_beam") << endl;
      
      try {
        tod_calc4 (par_file, pntarr, t_read_alm);
      } catch ( bad_alloc & e ) {
        cerr << corenum << " : Out of memory calling tod_calc4 for beam #  " << ibeam << endl;
        throw;
      }

      for (long ii=0; ii<totsize; ii++) todtmp[ii] += pntarr[5*ii+3];
    }

    for (long ii=0; ii<totsize; ii++) pntarr[5*ii+3] = todtmp[ii];

    todtmp.dealloc();
  }

  t_convolve = (t2=MPI_Wtime()) - t1 - t_read_alm;

  if ( corenum == 0 ) {
    std::cout << "Convolved TOD:" << std::endl;
    for ( long row=0; row < nsamp; ++row ) {
      std::cout << pntarr[row*5+4] << " " << pntarr[row*5+3] << std::endl;
    }
  }

  // Inject the data

  MPI_Barrier( MPI_COMM_WORLD );

  if ( corenum == 0 ) cout << endl << " ***** Injecting the data ***** " << endl << endl;

  // DEBUG : set the data to zero
  //for ( long i=0; i<pntarr.size(); i+=5 ) pntarr[i+3] = 0;
  // DEBUG : replace data with time stamps to identify where the problems are
  //for ( long i=0; i<pntarr.size(); i+=5 ) pntarr[i+3] = pntarr[i+4];
  // END DEBUG

  t1 = MPI_Wtime();

  //injector inj( eff_dir, file_pattern, det_pattern, scale, add, cores, corenum, 0 );

  t_init_write = (t2=MPI_Wtime()) - t1;

  //inj.inject( pntarr );

  t_write = (t1=MPI_Wtime()) - t2;

  MPI_Barrier( MPI_COMM_WORLD );

  // Optionally produce a binned map

  t1 = MPI_Wtime();

  //bin_map( par_file, pntarr, corenum );

  t_bin = (t2=MPI_Wtime()) - t1;

  MPI_Barrier( MPI_COMM_WORLD );

  // Report timing

  t_tot = MPI_Wtime() - t0;

  if ( corenum != 0 ) {
    MPI_Reduce( &t_init, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &t_read_tod, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &t_read_alm, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &t_convolve, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &t_bin, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &t_init_write, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( &t_write, NULL, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
  } else {
    MPI_Reduce( MPI_IN_PLACE, &t_init, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &t_read_tod, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &t_read_alm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &t_convolve, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &t_bin, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &t_init_write, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    MPI_Reduce( MPI_IN_PLACE, &t_write, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
    
    cout << endl;
    cout << "Total time:       " << t_tot << " s" << endl;
    cout << " -    IScalm init " << t_init/cores << " s" << endl;
    cout << " -     TOAST read " << t_read_tod/cores << " s" << endl;
    cout << " -       alm read " << t_read_alm/cores << " s" << endl;
    cout << " -IScalm convolve " << t_convolve/cores << " s" << endl;
    cout << " -  injector init " << t_init_write/cores << " s" << endl;
    cout << " - injector write " << t_write/cores << " s" << endl;
    cout << " -        bin map " << t_bin/cores << " s" << endl;
    cout << " -      imbalance " << t_tot - (t_init+t_read_tod+t_read_alm+t_convolve+t_bin+t_init_write+t_write)/cores << " s" << endl;
    cout << endl;
  }
}

int IScalm2TOD_module(int argc, const char **argv)
{
  module_startup ("IScalm2TOD", argc, const_cast<const char **>(argv), 2,
		  "<init object>", mpiMgr.master());
  
  if ( ! boost::filesystem::exists( argv[1] ) )
    throw runtime_error( string("Parameter file does not exist: ") + string( argv[1] ) );
       
  iohandle_current::Manager mng (argv[1]);
  paramfile par_file (mng.getParams(mpiMgr.master()));
  
  IScalm2TODEx(par_file);
  
  return 0;
}
