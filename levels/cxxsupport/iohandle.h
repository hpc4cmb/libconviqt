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

/*! \file iohandle.h
 *  This file contains the declaration of the general I/O interface
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2004, 2005, 2006, 2009, 2010 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_IOHANDLE_H
#define PLANCK_IOHANDLE_H

#include <string>
#include <map>
#include <vector>
#include "safe_ptr.h"
#include "datatypes.h"
#include "arr.h"
#include "safe_cast.h"

/*! \defgroup iohandlegroup Abstract data I/O interface */
/*! \{ */

/*! Abstract interface class for I/O of LevelS data objects */
class iohandle
  {
  private:
    virtual void setKeyVoid
      (const std::string &key, const void *value, PDT type) = 0;
    virtual void getKeyVoid
      (const std::string &key, void *value, PDT type) const = 0;

    virtual void readColumnRawVoid
      (int colnum, void *data, PDT type, long num, uint64 offset=0) const = 0;
    virtual void writeColumnRawVoid
      (int colnum, const void *data, PDT type, long num, uint64 offset=0) = 0;
    virtual void appendColumnRawVoid
      (int colnum, const void *data, PDT type, long num) = 0;

  public:
    virtual ~iohandle() {}

    /*! Creates a new data object called \a name, of type \a type. */
    virtual void create (const std::string &name, const std::string &type) = 0;
    /*! Opens an existing data object called \a name for reading. */
    virtual void open (const std::string &name) = 0;
    /*! Opens an existing data object called \a name with type \a type
        for reading. \a type need not be the exact type of the object, but may
        be a parent type as well. */
    virtual void open (const std::string &name, const std::string &type) = 0;
    /*! Closes the object, potentially committing all changes. */
    virtual void close() = 0;
    /*! Tries to flush all changes to the object to disk/DB. */
    virtual void flush() = 0;
    /*! Tries to reserve an object called \a name with type \a type and returns
        the real name of the new object. */
    virtual std::string getID (const std::string &name, const std::string &type)
      const = 0;
    /*! Deletes the object called \a name. */
    virtual void deleteObject (const std::string &name) const = 0;

    /*! Returns the name of the object connected to the handle */
    virtual std::string objectName() const = 0;
    /*! Returns the exact data type of the object. */
    virtual std::string objectType() const = 0;
    /*! Returns the number of data columns of the object. */
    virtual int numColumns() const = 0;

    /*! Returns \a true if the meta data key \a key is present in the
        object, else \a false. */
    virtual bool keyPresent (const std::string &key) const = 0;
    /*! Sets the value of \a key to \a value. */
    template<typename T> void setKey (const std::string &key, const T &value)
      { setKeyVoid (key, &value, planckType<T>()); }
    /*! Sets the value of \a key to \a value, potentially converting \a value
        to another datatype. */
    virtual void setKeyTypeless (const std::string &key,
      const std::string &value) = 0;
    /*! Returns the value associated with \a key. */
    template<typename T> T getKey (const std::string &key) const
      { T val; getKeyVoid (key, &val, planckType<T>()); return val; }
    /*! Returns the value associated with \a key in \a value. */
    template<typename T> void getKey (const std::string &key, T &value) const
      { getKeyVoid (key, &value, planckType<T>()); }
    /*! Returns the value associated with \a key as a string. */
    virtual std::string getKeyTypeless (const std::string &key) const = 0;
    /*! Deletes the key \a key. */
    virtual void deleteKey (const std::string &key) = 0;
    /*! Returns a list of all (non-reserved) key names in \a keys. */
    virtual void getAllKeys (std::vector<std::string> &keys) const = 0;

    /*! Executes a query described by \a query_string and returns
        the results in \a result. */
    virtual void query (const std::string &query_string,
      std::vector<std::string> &result) const = 0;

    /*! Returns \a true if the requested column name exists,
        else \a false. */
    virtual bool columnPresent (const std::string &colname) const = 0;
    /*! Returns the number of the column called \a colname. */
    virtual int columnNumber (const std::string &colname) const = 0;
    /*! Returns the name of the column with number \a colnum. */
    virtual std::string columnName (int colnum) const = 0;
    /*! Returns the data type of the column with number \a colnum. */
    virtual PDT columnType (int colnum) const = 0;
    /*! Returns the data type of the column called \a colname. */
    PDT columnType (const std::string &colname) const
      { return columnType (columnNumber(colname)); }
    /*! Returns the unit of the column with number \a colnum. */
    virtual std::string columnUnit (int colnum) const = 0;
    /*! Returns the unit of the column called \a colname. */
    std::string columnUnit (const std::string &colname) const
      { return columnUnit (columnNumber(colname)); }
    /*! Returns the number of entries in the column with number \a colnum. */
    virtual uint64 columnLength (int colnum) const = 0;
    /*! Returns the number of entries in the column called \a colname. */
    uint64 columnLength (const std::string &colname) const
      { return columnLength (columnNumber(colname)); }
    /*! Returns a hint regarding the number of values that should be
        read/written from/to \a colnum in a single call for best performance. */
    virtual uint64 efficientChunkSize (int colnum) const = 0;
    /*! Returns a hint regarding the number of values that should be
        read/written from/to \a colname in a single call for best performance.*/
    uint64 efficientChunkSize (const std::string &colname) const
      { return efficientChunkSize (columnNumber(colname)); }

    template<typename T> void readColumnRaw
      (int colnum, T *data, long num, uint64 offset=0) const
      {
      if (num>0)
        readColumnRawVoid (colnum, data, planckType<T>(), num, offset);
      }
    template<typename T> void readColumnRaw
      (const std::string &colname, T *data, long num, uint64 offset=0) const
      { readColumnRaw (columnNumber(colname), data, num, offset); }
    /*! Fills the array \a data with entries from column \a colnum, starting
        from column entry \a offset. */
    template<typename T> void readColumn
      (int colnum, arr<T> &data, uint64 offset=0) const
      { readColumnRaw (colnum, &(data[0]), data.size(), offset); }
    /*! Fills the array \a data with entries from column \a colname, starting
        from column entry \a offset. */
    template<typename T> void readColumn
      (const std::string &colname, arr<T> &data, uint64 offset=0) const
      { readColumn (columnNumber(colname), data, offset); }
    /*! Fills the array \a data with entries from column \a colnum, starting
        from column entry \a offset. */
    template<typename T> void readColumn
      (int colnum, std::vector<T> &data, uint64 offset=0) const
      { readColumnRaw (colnum, &(data[0]), data.size(), offset); }
    /*! Fills the array \a data with entries from column \a colname, starting
        from column entry \a offset. */
    template<typename T> void readColumn
      (const std::string &colname, std::vector<T> &data, uint64 offset=0) const
      { readColumn (columnNumber(colname), data, offset); }
    /*! Returns the entry with index \a offset of column \a colnum in
        \a data. */
    template<typename T> void readColumn
      (int colnum, T &data, uint64 offset=0) const
      { readColumnRaw (colnum, &data, 1, offset); }
    /*! Returns the entry with index \a offset of column \a colname in
        \a data. */
    template<typename T> void readColumn
      (const std::string &colname, T &data, uint64 offset=0) const
      { readColumn (columnNumber(colname), data, offset); }
    /*! Returns the entire contents of the column \a colnum in \a data.
        \a data is resized accordingly. */
    template<typename T> void readEntireColumn
      (int colnum, arr<T> &data) const
      {
      data.alloc(safe_cast<tsize>(columnLength(colnum)));
      readColumn (colnum, data);
      }
    /*! Returns the entire contents of the column \a colname in \a data.
        \a data is resized accordingly. */
    template<typename T> void readEntireColumn
      (const std::string &colname, arr<T> &data) const
      { readEntireColumn (columnNumber(colname), data); }
    /*! Returns the entire contents of the column \a colnum in \a data.
        \a data is resized accordingly. */
    template<typename T> void readEntireColumn
      (int colnum, std::vector<T> &data) const
      {
      data.resize(safe_cast<tsize>(columnLength(colnum)));
      readColumn (colnum, data);
      }
    /*! Returns the entire contents of the column \a colname in \a data.
        \a data is resized accordingly. */
    template<typename T> void readEntireColumn
      (const std::string &colname, std::vector<T> &data) const
      { readEntireColumn (columnNumber(colname), data); }

    template<typename T> void writeColumnRaw
      (int colnum, const T *data, long num, uint64 offset=0)
      {
      if (num>0)
        writeColumnRawVoid(colnum, data, planckType<T>(), num, offset);
      }
    template<typename T> void writeColumnRaw
      (const std::string &colname, const T *data, long num, uint64 offset=0)
      { writeColumnRaw(columnNumber(colname), data, num, offset); }
    /*! Writes the contents of \a data to the column \a colnum, starting at
        index \a offset. */
    template<typename T> void writeColumn
      (int colnum, const arr<T> &data, uint64 offset=0)
      { writeColumnRaw (colnum, &data[0], data.size(), offset); }
    /*! Writes the contents of \a data to the column \a colname, starting at
        index \a offset. */
    template<typename T> void writeColumn
      (const std::string &colname, const arr<T> &data, uint64 offset=0)
      { writeColumn (columnNumber(colname), data, offset); }
    /*! Writes the contents of \a data to the column \a colnum, starting at
        index \a offset. */
    template<typename T> void writeColumn
      (int colnum, const std::vector<T> &data, uint64 offset=0)
      { writeColumnRaw (colnum, &data[0], data.size(), offset); }
    /*! Writes the contents of \a data to the column \a colname, starting at
        index \a offset. */
    template<typename T> void writeColumn
      (const std::string &colname, const std::vector<T> &data, uint64 offset=0)
      { writeColumn (columnNumber(colname), data, offset); }
    /*! Writes \a data to the column \a colnum at index \a offset. */
    template<typename T> void writeColumn
      (int colnum, const T &data, uint64 offset=0)
      { writeColumnRaw (colnum, &data, 1, offset); }
    /*! Writes \a data to the column \a colname at index \a offset. */
    template<typename T> void writeColumn
      (const std::string &colname, const T &data, uint64 offset=0)
      { writeColumn (columnNumber(colname), data, offset); }

    template<typename T> void appendColumnRaw
      (int colnum, const T *data, long num)
      {
      if (num>0)
        appendColumnRawVoid(colnum, data, planckType<T>(), num);
      }
    template<typename T> void appendColumnRaw
      (const std::string &colname, const T *data, long num)
      { appendColumnRaw(columnNumber(colname), data, num); }
    /*! Appends the contents of \a data to the column \a colnum. */
    template<typename T> void appendColumn
      (int colnum, const arr<T> &data)
      { appendColumnRaw (colnum, &data[0], data.size()); }
    /*! Appends the contents of \a data to the column \a colname. */
    template<typename T> void appendColumn
      (const std::string &colname, const arr<T> &data)
      { appendColumn (columnNumber(colname), data); }
    /*! Appends the contents of \a data to the column \a colnum. */
    template<typename T> void appendColumn
      (int colnum, const std::vector<T> &data)
      { appendColumnRaw (colnum, &data[0], data.size()); }
    /*! Appends the contents of \a data to the column \a colname. */
    template<typename T> void appendColumn
      (const std::string &colname, const std::vector<T> &data)
      { appendColumn (columnNumber(colname), data); }
    /*! Appends \a data to the column \a colnum. */
    template<typename T> void appendColumn
      (int colnum, const T &data)
      { appendColumnRaw (colnum, &data, 1); }
    /*! Appends \a data to the column \a colname. */
    template<typename T> void appendColumn
      (const std::string &colname, const T &data)
      { appendColumn (columnNumber(colname), data); }
  };

/*! Class which can handle several different types of \a iohandle. */
class iohandleManager
  {
  private:
    typedef iohandle *(*pfunc)();
    typedef std::map<std::string,pfunc> maptype;
    maptype dict;
    std::string deflt;

    void getTypeAndName (const std::string &fullname, std::string &type,
      std::string &name) const;

  public:
    void registerType (const std::string &type, pfunc func)
      { dict[type] = func; }
    void registerDefaultType (const std::string &type, pfunc func)
      { dict[type] = func; deflt=type; }

    /*! Returns the name of the default \a iohandle implementation. */
    std::string defaultType() const
      { return deflt; }

    /*! Returns a pointer to a new \a iohandle of type \a type. */
    iohandle *makeHandle(const std::string &type) const;
    /*! Returns a pointer to a new \a iohandle of default type. */
    iohandle *makeHandle() const;

    /*! Returns a pointer to a new \a iohandle of default type, which is
        connected to a newly-created object of type \a ddltype, with a name
        based on \a fullname. If \a fullname starts with "FITS:", "HFIDMC:",
        or "TOODI:", the respective type of \a iohandle is created, otherwise
        the default type is used. */
    iohandle *createObject(const std::string &fullname,
      const std::string &ddltype) const;

    /*! Returns a pointer to a new \a iohandle of default type, which is
        connected for reading to an object of type \a ddltype, with a name
        based on \a fullname. If \a fullname starts with "FITS:", "HFIDMC:",
        or "TOODI:", the respective type of \a iohandle is created, otherwise
        the default type is used. */
    iohandle *openObject(const std::string &fullname,
      const std::string &ddltype) const;
    /*! Returns a pointer to a new \a iohandle of default type, which is
        connected for reading to an object with a name based on \a fullname.
        If \a fullname starts with "FITS:", "HFIDMC:", or "TOODI:", the
        respective type of \a iohandle is created, otherwise the default type
        is used. */
    iohandle *openObject(const std::string &fullname) const;

    /*! Deletes an object with a name based on \a fullname.
        If \a fullname starts with "FITS:", "HFIDMC:", or "TOODI:", the
        respective type of \a iohandle is used, otherwise the default type
        is used. */
    void deleteObject(const std::string &fullname) const;

    /*! Obtains an ID for a new object of type \a ddltype, with a name based on
        \a fullname. If \a fullname starts with "FITS:", "HFIDMC:", or "TOODI:",
        the respective type of \a iohandle is used, otherwise the default type
        is used. */
    std::string getID(std::string &fullname, const std::string &ddltype);
  };

extern iohandleManager HandleManager;

/*! \} */

#endif
