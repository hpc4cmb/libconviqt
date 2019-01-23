#ifndef __CCONVIQT_H__

#define __CCONVIQT_H__

#ifdef __cplusplus
extern "C" {
#endif

  const int LMAXMAX=100000;

  // Beam

  void *conviqt_beam_new();

  int conviqt_beam_del(void *ptr);

  int conviqt_beam_read(void *ptr, long beamlmax, long beammmax, char pol,
			char *infile_beam, MPI_Comm comm);

  int conviqt_beam_lmax(void *ptr);

  int conviqt_beam_mmax(void *ptr);

  double conviqt_beam_normalize(void *ptr);

  // Sky

  void *conviqt_sky_new();

  int conviqt_sky_del(void *ptr);

  int conviqt_sky_read(void *ptr, long skylmax, char pol, char *infile_sky,
		       double fwhm_deconv_sky, MPI_Comm comm);

  int conviqt_sky_lmax(void *ptr);

  int conviqt_sky_remove_monopole(void *ptr);

  int conviqt_sky_remove_dipole(void *ptr);

  // Detector

  void *conviqt_detector_new();

  void *conviqt_detector_new_with_id(char *det_id);

  int conviqt_detector_del(void *ptr);

  int conviqt_detector_set_epsilon(void *ptr, double epsilon);

  int conviqt_detector_get_epsilon(void *ptr, double *epsilon);

  int conviqt_detector_get_id(void *ptr, char *det_id);

  // Pointing

  void *conviqt_pointing_new();

  int conviqt_pointing_del(void *ptr);

  int conviqt_pointing_alloc(void *ptr, long n);

  double *conviqt_pointing_data(void *ptr);

  // Convolver

  void *conviqt_convolver_new(void *skyptr, void *beamptr, void *detptr, char pol,
			      long lmax, long beammmax, long order, MPI_Comm comm);

  int conviqt_convolver_del(void *ptr);

  int conviqt_convolver_set_sky(void *cnvptr, void *skyptr);

  int conviqt_convolver_set_beam(void *cnvptr, void *beamptr);

  int conviqt_convolver_set_detector(void *cnvptr, void *detptr);

  int conviqt_convolver_convolve(void *cnvptr, void *pntptr, char calibrate);

#ifdef __cplusplus
}
#endif

#endif
