#ifndef LEVELS_FACILITIES_H
#define LEVELS_FACILITIES_H

int syn_alm_cxx_module (int argc, const char **argv);
int alm2map_cxx_module (int argc, const char **argv);
int anafast_cxx_module (int argc, const char **argv);
int totalconvolve_cxx_module (int argc, const char **argv);
int alm2grid_module (int argc, const char **argv);
int calc_powspec_module (int argc, const char **argv);
int hotspots_cxx_module (int argc, const char **argv);
int map2tga_module (int argc, const char **argv);
int median_filter_cxx_module (int argc, const char **argv);
int mult_alm_module (int argc, const char **argv);
int smoothing_cxx_module (int argc, const char **argv);
int udgrade_cxx_module (int argc, const char **argv);
int udgrade_harmonic_cxx_module (int argc, const char **argv);
int conviqt_v3_module (int argc, const char **argv);
int conviqt_v4_module (int argc, const char **argv);
int deconviqt_module (int argc, const char **argv);
int multimod_module (int argc, const char **argv);
int skymixer3_module (int argc, const char **argv);
int almmixer_module (int argc, const char **argv);

int LevelS_module (const char *module, int argc, const char **argv);

#endif
