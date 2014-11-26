#include "conviqt.hpp"

double xpow (int expo, double val)
  { return (expo&1) ? -val : val; }

double wignerCalc00_halfpi( tsize n )
{
  if (n&1) return 0.;
  double rec=1.;
  for (tsize k=1; k<n; k+=2)
    rec *= -(k/(k+1.));
  return rec;
}

double wignerCalc00( double theta, tsize n )
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

void wignerCalcGeneral( tsize n, tsize mmax, double theta, arr2<double> &d )
{
  
  try {
    d.fast_alloc(2*n+1,2*mmax+1);
  } catch ( std::bad_alloc & e ) {
    std::cerr << " Out of memory allocating " << (2*n+1)*(2*mmax+1)*8./1024/1024 << "MB for wignerCalcGeneral" << std::endl;
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
	  flag1 = true;
      d[n-i][mmax-i] = -sqtab[n-(i-1)]/sqtab[n-i]*d[n-i+2][mmax-i];
	}
      if (abs(d00_2) < 1e-14)
      {
	flag2 = true;
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

void wignerCalcHalfpi( tsize n, tsize mmax, arr2<double> &d )
{
  
  try {
    d.fast_alloc(2*n+1,2*mmax+1);
  } catch ( std::bad_alloc & e ) {
    std::cerr << " Out of memory allocating " << (2*n+1)*(2*mmax+1)*8./1024/1024 << "MB for wignerCalcHalfpi" << std::endl;
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

void wignerCalc( tsize n, tsize mmax, double theta, arr2<double> &d )
{
  try {
    d.fast_alloc (2*n+1,2*mmax+1);
  } catch ( std::bad_alloc & e ) {
    std::cerr << " Out of memory allocating " << (2*n+1)*(2*mmax+1)*8./1024/1024 << "MB for wignerCalc" << std::endl;
    throw;
  }
  
  tdiff mm = mmax;
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
      double anglediff = 2.*twopi;
      while( abs(anglediff) > twopi )
	{
	  anglefac++;
	  anglediff = theta - thetasign*anglefac*twopi;
	}
      theta = theta - thetasign*anglefac*twopi;
    }
  if ( (abs_approx(theta, 0., 1e-14)) || (abs_approx(fabs(theta), twopi, conv_acc) ))
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
  if ( abs_approx(theta, halfpi, 1e-14) )
    {
      wignerCalcHalfpi( n, mmax, d );
      return;
    }
  wignerCalcGeneral( n, mmax, theta, d );
}
