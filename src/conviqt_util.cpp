#include "conviqt.hpp"

void sift_down_DL( arr<double> &ra, arr<long> &brr, const int l, const int r ) //FROM P-340 OF NR
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

  
void hpsort_DL( arr<double> &ra, arr<long> &brr ) // FROM P-340 OF NR
{
  int i;
  int n = ra.size();
  for(i=n/2-1; i>=0; i--)
    sift_down_DL(ra,brr,i,n-1);
  for (i=n-1; i>0; i--)
    {
      std::swap(ra[0], ra[i]);
      std::swap(brr[0], brr[i]);
      sift_down_DL(ra,brr,0,i-1);
    }
}

void sift_down_arrTheta( arr<double> &ra, const int l, const int r ) // FROM P-340 OF NR
{
  int j, jold;
  double a,c,d,f;
  double b;
  a = ra[5*l+1];
  b = ra[5*l+0];
  f = ra[5*l+4];
  c = ra[5*l+2];
  d = ra[5*l+3];
  jold = l;
  j = 2*l + 1;
  while(j<=r)
    {
      if (j<r && ra[5*j+1] < ra[5*(j+1)+1]) j++;
      if (a>= ra[5*j+1]) break;
      ra[5*jold+1] = ra[5*j+1];
      ra[5*jold+2] = ra[5*j+2];
      ra[5*jold+0] = ra[5*j+0];
      ra[5*jold+4] = ra[5*j+4];
      ra[5*jold+3] = ra[5*j+3];
      jold = j;
      j = 2*j + 1;
    }
  ra[5*jold+1] = a;
  ra[5*jold+2] = c;
  ra[5*jold+0] = b;
  ra[5*jold+4] = f;
  ra[5*jold+3] = d;
}


void hpsort_arrTheta( arr<double> &ra ) // FROM P-340 OF NR
{
  int i;
  int n = ra.size()/5;
  for(i=n/2-1; i>=0; i--)
    sift_down_arrTheta(ra,i,n-1);
  for (i=n-1; i>0; i--)
    {
      std::swap(ra[1], ra[5*i+1]);
      std::swap(ra[0], ra[5*i+0]);
      std::swap(ra[4], ra[5*i+4]);
      std::swap(ra[2], ra[5*i+2]);
      std::swap(ra[3], ra[5*i+3]);
      sift_down_arrTheta(ra,0,i-1);
    }
}

void sift_down_arrTOD( arr<double> &ra, const int l, const int r ) // FROM P-340 OF NR
{
  int j, jold;
  double a, c, d, f;
  double b;
  a = ra[5*l+3];
  b = ra[5*l+0];
  f = ra[5*l+4];
  c = ra[5*l+2];
  d = ra[5*l+1];
  jold = l;
  j = 2*l + 1;
  while(j<=r)
    {
      if (j<r && ra[5*j+3] < ra[5*(j+1)+3]) j++;
      if (a>= ra[5*j+3]) break;
      ra[5*jold+3] = ra[5*j+3];
      ra[5*jold+2] = ra[5*j+2];
      ra[5*jold+0] = ra[5*j+0];
      ra[5*jold+4] = ra[5*j+4];
      ra[5*jold+1]  =ra[5*j+1];
      jold = j;
      j = 2*j+1;
    }
  ra[5*jold+3] = a;
  ra[5*jold+2] = c;
  ra[5*jold+0] = b;
  ra[5*jold+4] = f;
  ra[5*jold+1] = d;
}


void hpsort_arrTOD( arr<double> &ra ) // FROM P-340 OF NR
{
  int i;
  int n = ra.size()/5;
  for(i=n/2-1; i>=0; i--)
    sift_down_arrTOD(ra,i,n-1);
  for (i=n-1; i>0; i--)
    {
      std::swap(ra[3], ra[5*i+3]);
      std::swap(ra[0], ra[5*i+0]);
      std::swap(ra[4], ra[5*i+4]);
      std::swap(ra[2], ra[5*i+2]);
      std::swap(ra[1], ra[5*i+1]);
      sift_down_arrTOD( ra, 0, i-1 );
    }
}


void sift_down_arrTime( arr<double> &ra, const int l, const int r ) // FROM P-340 OF NR
{
  int j, jold;
  double a, c, d, f;
  double b;
  a = ra[5*l+4];
  b = ra[5*l+0];
  f = ra[5*l+1];
  c = ra[5*l+2];
  d = ra[5*l+3];
  jold = l;
  j = 2*l + 1;
  while(j<=r)
    {
      if (j<r && ra[5*j+4] < ra[5*(j+1)+4]) j++;
      if (a>= ra[5*j+4]) break;
      ra[5*jold+4] = ra[5*j+4];
      ra[5*jold+2] = ra[5*j+2];
      ra[5*jold+0] = ra[5*j+0];
      ra[5*jold+1] = ra[5*j+1];
      ra[5*jold+3] = ra[5*j+3];
      jold = j;
      j = 2*j+1;
    }
  ra[5*jold+4] = a;
  ra[5*jold+2] = c;
  ra[5*jold+0] = b;
  ra[5*jold+1] = f;
  ra[5*jold+3] = d;
}

  
void hpsort_arrTime( arr<double> &ra ) // FROM P-340 OF NR
{
  int i;
  int n = ra.size()/5;
  for(i=n/2-1; i>=0; i--)
    sift_down_arrTime(ra, i, n-1);
  for (i=n-1; i>0; i--)
    {
      std::swap(ra[4], ra[5*i+4]);
      std::swap(ra[0], ra[5*i+0]);
      std::swap(ra[1], ra[5*i+1]);
      std::swap(ra[2], ra[5*i+2]);
      std::swap(ra[3], ra[5*i+3]);
      sift_down_arrTime( ra, 0, i-1 );
    }
}


void sift_down_DDcm( arr<double> &ra, arr<double> &brr, const int l, const int r ) // FROM P-340 OF NR
{
  int j, jold;
  double a;
  double b;
  a = ra[l];
  b = brr[l];
  jold = l;
  j = 2*l+1;
  while(j<=r)
    {
      if (j<r && ra[j] < ra[j+1]) j++;
      if (a>= ra[j]) break;
      ra[jold]=ra[j];
      brr[jold]=brr[j];
      jold = j;
      j = 2*j + 1;
    }
  ra[jold] = a;
  brr[jold] = b;
}


void hpsort_DDcm( arr<double> &ra, arr<double> &brr ) //FROM P-340 OF NR
{
  int i;
  int n = ra.size();
  for(i=n/2-1; i>=0; i--)
    sift_down_DDcm( ra, brr, i, n-1 );
  for (i=n-1; i>0; i--)
    {
      std::swap(ra[0], ra[i]);
      std::swap(brr[0], brr[i]);
      sift_down_DDcm( ra, brr,0, i-1 );
    }
}

