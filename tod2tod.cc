#include "levels_facilities.h"
#include "error_handling.h"

int tod2tod_module(int argc, const char **argv);

int main (int argc, const char **argv)
  {
  PLANCK_DIAGNOSIS_BEGIN
  tod2tod_module (argc, argv);
  PLANCK_DIAGNOSIS_END
  }
