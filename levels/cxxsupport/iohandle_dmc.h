/*
 *  This file contains the declaration of the DMC I/O interface
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2004-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

#ifndef PLANCK_IOHANDLE_DMC_H
#define PLANCK_IOHANDLE_DMC_H

#include <string>
#include "iohandle.h"
#include "paramfile.h"

class iohandle_dmc: public iohandle
  {
  private:
    bool connected;
    bool readonly;
    void *handle;
    std::string type_;
    arr<std::string> colname_, colunit_;
    arr<PDT> coltype_;

    void assert_connected() const;
    void assert_write() const;
    void assert_read() const;
    void clean_all();
    std::string getStringKey_connected (const std::string &key) const;
    void fill_cache();

    virtual void setKeyVoid
      (const std::string &key, const void *value, PDT type);
    virtual void getKeyVoid
      (const std::string &key, void *value, PDT type) const;

    virtual void readColumnRawVoid
      (int colnum, void *data, PDT type, long num, uint64 offset=0) const;
    virtual void writeColumnRawVoid
      (int colnum, const void *data, PDT type, long num, uint64 offset=0);
    virtual void appendColumnRawVoid
      (int colnum, const void *data, PDT type, long num);

  public:
    class Manager
      {
      private:
        std::map<std::string,std::string> paramMap;

      public:
        Manager (const std::string &properties);
        Manager (int argc, const char **argv);
        ~Manager();
        paramfile getParams (bool verbose=true) const;
      };

    static iohandle *makeHandle();
    static std::string Type();

    iohandle_dmc ();
    virtual ~iohandle_dmc();

    virtual void create (const std::string &name, const std::string &type);
    virtual void open (const std::string &name);
    virtual void open (const std::string &name, const std::string &type);
    virtual void close();
    virtual void flush();
    virtual std::string getID (const std::string &name, const std::string &type)
      const;
    virtual void deleteObject (const std::string &name) const;

    virtual std::string objectName() const;
    virtual std::string objectType() const;
    virtual int numColumns() const;

    virtual bool keyPresent (const std::string &key) const;
    virtual void setKeyTypeless (const std::string &key,
      const std::string &value);
    virtual std::string getKeyTypeless (const std::string &key) const;
    virtual void deleteKey (const std::string &key);
    virtual void getAllKeys (std::vector<std::string> &keys) const;
    virtual void query (const std::string &query_string,
      std::vector<std::string> &result) const;

    virtual bool columnPresent (const std::string &colname) const;
    virtual int columnNumber (const std::string &colname) const;
    virtual std::string columnName (int colnum) const;
    virtual PDT columnType (int colnum) const;
    virtual std::string columnUnit (int colnum) const;
    virtual uint64 columnLength (int colnum) const;
    virtual uint64 efficientChunkSize (int colnum) const;
  };

#endif
