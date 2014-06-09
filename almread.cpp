#include "almread.hpp"

int main( int argc, char **argv ) {

  int err;
  //err = MPI_Init( &argc, &argv ); Interestingly, MPI_Init fails if the code instantiates a single Alm object...
  //if (err != MPI_SUCCESS) throw runtime_error( "MPI Initialization failed" );    

  int ntask, id;
  err = MPI_Comm_size(MPI_COMM_WORLD, &ntask);
  err = MPI_Comm_rank(MPI_COMM_WORLD, &id);

  if ( argc != 4 ) {
    if ( id == 0 ) cout << "Usage: " << argv[0] << " <alm in (fits)> <alm out (binary)> <nread_task>" << endl;
    err = MPI_Finalize();
    return 0;
  }

  string infile_sky( argv[1] );
  string outfile_sky( argv[2] );
  int nread_task = min( ntask, atoi(argv[3]) );

  int groupsize = ceil( ntask / (double)nread_task );
  MPI_Comm iogroup;
  int ntask_group, id_group;
  if ( nread_task > 1 ) {
    int igroup = id / groupsize;
    err = MPI_Comm_split( MPI_COMM_WORLD, igroup, id, &iogroup );
    err = MPI_Comm_size(iogroup, &ntask_group);
    err = MPI_Comm_rank(iogroup, &id_group);
  } else {
    iogroup = MPI_COMM_WORLD;
    ntask_group = ntask;
    id_group = id;
  }

  if ( id == 0 ) cout << "almread started with " << ntask << " MPI processes organized into " << nread_task << " I/O groups " << endl;

  long lmax=6000;
  Alm<xcomplex<float> > slmInputT(lmax,lmax);
  Alm<xcomplex<float> > slmInputG(lmax,lmax);
  Alm<xcomplex<float> > slmInputC(lmax,lmax);
  bool pol;

  //string infile_sky = "ffp7_fg_030_alm.fits";
  //string outfile_sky = "ffp7_fg_030_alm.dat";
  //long slmsize_OTF = (lmax+1)*(lmax+2);
  long slmsize = Alm<xcomplex<float> >::Num_Alms(lmax, lmax);

  if ( id == 0 ) cout << "Reading "<< infile_sky << ", lmax = " << lmax << endl;

  MPI_Barrier( MPI_COMM_WORLD );

  double t1 = MPI_Wtime();

  if ( id_group == 0 ) {  
    read_Alm_from_dmc( infile_sky, slmInputT, slmInputG, slmInputC, lmax, lmax );
  }

  MPI_Barrier( MPI_COMM_WORLD );

  double t2 = MPI_Wtime();

  if ( ntask_group > 1 ) {
    err = MPI_Bcast( &(slmInputT(0,0).re), slmsize*2, MPI_FLOAT, 0, iogroup );
    err = MPI_Bcast( &(slmInputG(0,0).re), slmsize*2, MPI_FLOAT, 0, iogroup );
    err = MPI_Bcast( &(slmInputC(0,0).re), slmsize*2, MPI_FLOAT, 0, iogroup );
  }

  MPI_Barrier( MPI_COMM_WORLD );

  double t3 = MPI_Wtime();

  cout << id << " : FITS read time = " << t2 - t1 << "s, Bcast time = " << t3 - t2 << "s. Tot = " << t3 - t1 << endl;

  long filesize = sizeof(xcomplex<float>) / sizeof(char) * slmsize;

  if ( id == 0 ) {    
    cout << "Writing " << 3*filesize/pow(2.,20.) << "MB to "<< outfile_sky << endl;
    double t1 = MPI_Wtime();  

    write_alm_bin<float>( outfile_sky, slmInputT, slmInputG, slmInputC, lmax, lmax );

    double t2 = MPI_Wtime();
    cout << "Binary write time : " << t2 - t1 << "s " << endl;
  }

  MPI_Barrier( MPI_COMM_WORLD );

  if ( id == 0 ) cout << "Reading " << outfile_sky << endl;

  Alm<xcomplex<float> > slmInputT_v2(lmax,lmax);
  Alm<xcomplex<float> > slmInputG_v2(lmax,lmax);
  Alm<xcomplex<float> > slmInputC_v2(lmax,lmax);
 
  MPI_Barrier( MPI_COMM_WORLD );

  t1 = MPI_Wtime();

  if ( id_group == 0 ) {  
    read_alm_bin<float>( outfile_sky, slmInputT_v2, slmInputG_v2, slmInputC_v2, lmax, lmax );
  }

  MPI_Barrier( MPI_COMM_WORLD );

  t2 = MPI_Wtime();

  if ( ntask_group > 1 ) {
    err = MPI_Bcast( &(slmInputT_v2(0,0).re), slmsize*2, MPI_FLOAT, 0, iogroup );
    err = MPI_Bcast( &(slmInputG_v2(0,0).re), slmsize*2, MPI_FLOAT, 0, iogroup );
    err = MPI_Bcast( &(slmInputC_v2(0,0).re), slmsize*2, MPI_FLOAT, 0, iogroup );
  }

  MPI_Barrier( MPI_COMM_WORLD );

  t3 = MPI_Wtime();

  cout << id << " : Binary read time = " << t2 - t1 << "s, Bcast time = " << t3 - t2 << "s. Tot = " << t3 - t1 << endl;

  // Compute the RMS difference between the two versions of the s_lm
  
  double rmsT = compare_alm<float>( slmInputT, slmInputT_v2 );
  double rmsG = compare_alm<float>( slmInputG, slmInputG_v2 );
  double rmsC = compare_alm<float>( slmInputC, slmInputC_v2 );

  cout << id << " : Relative errors in slmT = " << rmsT << ", slmG = " << rmsG << ", slmC = " << rmsC << endl;

  err = MPI_Barrier(MPI_COMM_WORLD);
  err = MPI_Finalize();

  return 0;

}
