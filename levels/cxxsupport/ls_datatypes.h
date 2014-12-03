/*
 *  This file is part of libcxxsupport.
 *
 *  libcxxsupport is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  libcxxsupport is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with libcxxsupport; if not, write to the Free Software
 *  Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
 */

/*
 *  libcxxsupport is being developed at the Max-Planck-Institut fuer Astrophysik
 *  and financially supported by the Deutsches Zentrum fuer Luft- und Raumfahrt
 *  (DLR).
 */

/*! \file datatypes.h
 *  This file defines various platform-independent data types.
 *  If any of the requested types is not available, compilation aborts
 *  with an error (unfortunately a rather obscure one).
 *
 *  Copyright (C) 2004-2011 Max-Planck-Society
 *  \author Martin Reinecke
 */

#ifndef LEVELS_DATATYPES_H
#define LEVELS_DATATYPES_H

#include <string>
#include <cstddef>
#include "ls_error_handling.h"

namespace levels {

// Template magic to select the proper data types. These templates
// should not be used outside this file.

template <typename T, bool equalSize> struct sizeChooserHelper__
  { typedef void TYPE; };

template <typename T> struct sizeChooserHelper__<T,true>
  { typedef T TYPE; };

template <typename T1, typename T2, typename T3> struct sizeChooserHelper2__
  { typedef T1 TYPE; };

template <typename T2, typename T3> struct sizeChooserHelper2__ <void, T2, T3>
  { typedef T2 TYPE; };

template <typename T3> struct sizeChooserHelper2__ <void, void, T3>
  { typedef T3 TYPE; };

template <> struct sizeChooserHelper2__ <void, void, void>
  { };

template <int sz, typename T1, typename T2=char, typename T3=char>
  struct sizeChooser__
  {
  typedef typename sizeChooserHelper2__
    <typename sizeChooserHelper__<T1,sizeof(T1)==sz>::TYPE,
     typename sizeChooserHelper__<T2,sizeof(T2)==sz>::TYPE,
     typename sizeChooserHelper__<T3,sizeof(T3)==sz>::TYPE >::TYPE TYPE;
  };

typedef signed char int8;
typedef unsigned char uint8;

typedef sizeChooser__<2, short, int>::TYPE
  int16;
typedef sizeChooser__<2, unsigned short, unsigned int>::TYPE
  uint16;

typedef sizeChooser__<4, int, long, short>::TYPE
  int32;
typedef sizeChooser__<4, unsigned int, unsigned long, unsigned short>::TYPE
  uint32;

typedef sizeChooser__<8, long, long long>::TYPE
  int64;
typedef sizeChooser__<8, unsigned long, unsigned long long>::TYPE
  uint64;

typedef sizeChooser__<4, float, double>::TYPE
  float32;
typedef sizeChooser__<8, double, long double>::TYPE
  float64;

/*! unsigned integer type which should be used for array sizes */
typedef std::size_t tsize;
/*! signed integer type which should be used for relative array indices */
typedef std::ptrdiff_t tdiff;

/*! mapping of Planck data types to integer constants */
enum PDT {
       LEVELS_INT8    =  0,
       LEVELS_UINT8   =  1,
       LEVELS_INT16   =  2,
       LEVELS_UINT16  =  3,
       LEVELS_INT32   =  4,
       LEVELS_UINT32  =  5,
       LEVELS_INT64   =  6,
       LEVELS_UINT64  =  7,
       LEVELS_FLOAT32 =  8,
       LEVELS_FLOAT64 =  9,
       LEVELS_BOOL    = 10,
       LEVELS_STRING  = 11,
       LEVELS_INVALID = -1 };

/*! Returns the \a PDT constant associated with \a T. */
template<typename T> inline PDT planckType()
  { levels_fail(T::UNSUPPORTED_DATA_TYPE); }
template<> inline PDT planckType<int8>       () { return LEVELS_INT8;   }
template<> inline PDT planckType<uint8>      () { return LEVELS_UINT8;  }
template<> inline PDT planckType<int16>      () { return LEVELS_INT16;  }
template<> inline PDT planckType<uint16>     () { return LEVELS_UINT16; }
template<> inline PDT planckType<int32>      () { return LEVELS_INT32;  }
template<> inline PDT planckType<uint32>     () { return LEVELS_UINT32; }
template<> inline PDT planckType<int64>      () { return LEVELS_INT64;  }
template<> inline PDT planckType<uint64>     () { return LEVELS_UINT64; }
template<> inline PDT planckType<float32>    () { return LEVELS_FLOAT32;}
template<> inline PDT planckType<float64>    () { return LEVELS_FLOAT64;}
template<> inline PDT planckType<bool>       () { return LEVELS_BOOL;   }
template<> inline PDT planckType<std::string>() { return LEVELS_STRING; }

/*! Returns the size (in bytes) of the Planck data type \a type. */
inline int type2size (PDT type)
  {
  switch (type)
    {
    case LEVELS_INT8   :
    case LEVELS_UINT8  :
    case LEVELS_BOOL   :
    case LEVELS_STRING : return 1;
    case LEVELS_INT16  :
    case LEVELS_UINT16 : return 2;
    case LEVELS_INT32  :
    case LEVELS_UINT32 :
    case LEVELS_FLOAT32: return 4;
    case LEVELS_INT64  :
    case LEVELS_UINT64 :
    case LEVELS_FLOAT64: return 8;
    default:
      levels_fail ("type2size: unsupported data type");
    }
  }

/*! Converts the string \a type to a \a PDT. */
inline PDT string2type(const std::string &type)
  {
  if (type=="FLOAT64") return LEVELS_FLOAT64;
  if (type=="FLOAT32") return LEVELS_FLOAT32;
  if (type=="INT8")    return LEVELS_INT8;
  if (type=="UINT8")   return LEVELS_UINT8;
  if (type=="INT16")   return LEVELS_INT16;
  if (type=="UINT16")  return LEVELS_UINT16;
  if (type=="INT32")   return LEVELS_INT32;
  if (type=="UINT32")  return LEVELS_UINT32;
  if (type=="INT64")   return LEVELS_INT64;
  if (type=="UINT64")  return LEVELS_UINT64;
  if (type=="BOOL")    return LEVELS_BOOL;
  if (type=="STRING")  return LEVELS_STRING;
  levels_fail ("string2type: unknown data type '"+type+"'");
  }

/*! Converts the Planck data type \a type to a C string. */
inline const char *type2string (PDT type)
  {
  switch (type)
    {
    case LEVELS_INT8   : return "INT8";
    case LEVELS_UINT8  : return "UINT8";
    case LEVELS_INT16  : return "INT16";
    case LEVELS_UINT16 : return "UINT16";
    case LEVELS_INT32  : return "INT32";
    case LEVELS_UINT32 : return "UINT32";
    case LEVELS_INT64  : return "INT64";
    case LEVELS_UINT64 : return "UINT64";
    case LEVELS_FLOAT32: return "FLOAT32";
    case LEVELS_FLOAT64: return "FLOAT64";
    case LEVELS_BOOL   : return "BOOL";
    case LEVELS_STRING : return "STRING";
    default:
      levels_fail ("type2string: unsupported data type");
    }
  }

/*! Returns a C string describing the data type \a T. */
template<typename T> inline const char *type2typename ()
  { levels_fail(T::UNSUPPORTED_DATA_TYPE); }
template<> inline const char *type2typename<signed char> ()
  { return "signed char"; }
template<> inline const char *type2typename<unsigned char> ()
  { return "unsigned char"; }
template<> inline const char *type2typename<short> ()
  { return "short"; }
template<> inline const char *type2typename<unsigned short> ()
  { return "unsigned short"; }
template<> inline const char *type2typename<int> ()
  { return "int"; }
template<> inline const char *type2typename<unsigned int> ()
  { return "unsigned int"; }
template<> inline const char *type2typename<long> ()
  { return "long"; }
template<> inline const char *type2typename<unsigned long> ()
  { return "unsigned long"; }
template<> inline const char *type2typename<long long> ()
  { return "long long"; }
template<> inline const char *type2typename<unsigned long long> ()
  { return "unsigned long long"; }
template<> inline const char *type2typename<float> ()
  { return "float"; }
template<> inline const char *type2typename<double> ()
  { return "double"; }
template<> inline const char *type2typename<long double> ()
  { return "long double"; }
template<> inline const char *type2typename<bool> ()
  { return "bool"; }
template<> inline const char *type2typename<std::string> ()
  { return "std::string"; }

/*! mapping of "native" data types to integer constants */
enum NDT {
       NAT_CHAR,
       NAT_SCHAR,
       NAT_UCHAR,
       NAT_SHORT,
       NAT_USHORT,
       NAT_INT,
       NAT_UINT,
       NAT_LONG,
       NAT_ULONG,
       NAT_LONGLONG,
       NAT_ULONGLONG,
       NAT_FLOAT,
       NAT_DOUBLE,
       NAT_LONGDOUBLE,
       NAT_BOOL,
       NAT_STRING };

/*! Returns the \a NDT constant associated with \a T. */
template<typename T> inline NDT nativeType()
  { levels_fail(T::UNSUPPORTED_DATA_TYPE); }
template<> inline NDT nativeType<char>              () { return NAT_CHAR;      }
template<> inline NDT nativeType<signed char>       () { return NAT_SCHAR;     }
template<> inline NDT nativeType<unsigned char>     () { return NAT_UCHAR;     }
template<> inline NDT nativeType<short>             () { return NAT_SHORT;     }
template<> inline NDT nativeType<unsigned short>    () { return NAT_USHORT;    }
template<> inline NDT nativeType<int>               () { return NAT_INT;       }
template<> inline NDT nativeType<unsigned int>      () { return NAT_UINT;      }
template<> inline NDT nativeType<long>              () { return NAT_LONG;      }
template<> inline NDT nativeType<unsigned long>     () { return NAT_ULONG;     }
template<> inline NDT nativeType<long long>         () { return NAT_LONGLONG;  }
template<> inline NDT nativeType<unsigned long long>() { return NAT_ULONGLONG; }
template<> inline NDT nativeType<float>             () { return NAT_FLOAT;     }
template<> inline NDT nativeType<double>            () { return NAT_DOUBLE;    }
template<> inline NDT nativeType<long double>       () { return NAT_LONGDOUBLE;}
template<> inline NDT nativeType<bool>              () { return NAT_BOOL;      }
template<> inline NDT nativeType<std::string>       () { return NAT_STRING;    }

/*! Returns the size (in bytes) of the native data type \a type. */
inline int ndt2size (NDT type)
  {
  switch (type)
    {
    case NAT_CHAR      :
    case NAT_SCHAR     :
    case NAT_UCHAR     : return sizeof(char);
    case NAT_SHORT     :
    case NAT_USHORT    : return sizeof(short);
    case NAT_INT       :
    case NAT_UINT      : return sizeof(int);
    case NAT_LONG      :
    case NAT_ULONG     : return sizeof(long);
    case NAT_LONGLONG  :
    case NAT_ULONGLONG : return sizeof(long long);
    case NAT_FLOAT     : return sizeof(float);
    case NAT_DOUBLE    : return sizeof(double);
    case NAT_LONGDOUBLE: return sizeof(long double);
    case NAT_BOOL      : return sizeof(bool);
    default:
      levels_fail ("ndt2size: unsupported data type");
    }
  }

} // namespace levels

#endif
