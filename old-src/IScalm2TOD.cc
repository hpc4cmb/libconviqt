#include "levels_facilities.h"
#include "error_handling.h"

int IScalm2TOD_module(int argc, const char **argv);

int main (int argc, const char **argv)
  {
  PLANCK_DIAGNOSIS_BEGIN
  IScalm2TOD_module (argc, argv);
  PLANCK_DIAGNOSIS_END
  }
