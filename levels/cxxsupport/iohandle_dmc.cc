/*
 *  This file contains the implementation of the DMC I/O interface
 *  used by the Planck LevelS package.
 *
 *  Copyright (C) 2004-2011 Max-Planck-Society
 *  Author: Martin Reinecke
 */

// make sure exactly one DMC is specified
#if ((defined(USE_TOODI)) && (defined(USE_HFIDMC)))
#error support for both DMCs requested
#endif

#if ((defined(USE_TOODI)) || (defined(USE_HFIDMC)))

#include <cstdlib>
#include <cstring>
#include "iohandle_dmc.h"
#include "string_utils.h"
#ifdef USE_TOODI
#include "LLIOcommon_t.h"
#else
#include "HL2_COMMONLLIO/PIOCommonLLIO.h"
#endif

using namespace std;

//#define CHECK(call) { cerr << "doing " << #call << endl; int status = call; if (status!=0) { cerr << "DMC error " << status << endl; planck_fail("Stopped");} cerr << "done" << endl; }
//#define CHECK(call) { int status = call; if (status!=0) { cerr << "DMC error while doing " << endl << #call << endl << "exit status: " << status << endl; char *msg; dmcExceptionString(status,&msg); cerr << msg << endl; free(msg); planck_fail("Stopped");} }
#define CHECK(call) { int status = call; if (status!=0) { cerr << "DMC error while doing " << endl << #call << endl << "exit status: " << status << endl; planck_fail("Stopped");} }

#define HANDLE *(static_cast<ObjectHandle *> (handle))
#define PHANDLE (static_cast<ObjectHandle *> (handle))

namespace {

SessionHandle session;
ObjectHandle initHandle;

bool keyReserved (const string &key)
  {
#ifdef USE_TOODI
  const char *exclist[] = { "acceStat","authorId","frzeStat","objType",
    "creaDate","intrStat","authAffi","version","frstName","dbUserId",
    "dmcVersn","tranStat","name","ddlVersn","lastName","deleStat","MRun",
    "Refcount","accessCounter","cpu.per.node","cpu.type","dataRemoved",
    "machine.type","max.cputime","memory.size","moduleType","node.count",
    "os.name","os.type","parentMDO" };
#else
  const char *exclist[] = {};
#endif
  const tsize n_exc = sizeof(exclist)/sizeof(char *);

  for (tsize i=0; i<n_exc; ++i)
    if (equal_nocase(key,exclist[i])) return true;
  return false;
  }

void getSimParams (ObjectHandle &inp, map<string,string> &par)
  {
  int nkeys;
  dmcChar **keys;
  CHECK(dmcGetMetaKeys(inp, &nkeys, &keys));
  for (int m=0; m<nkeys; ++m)
    {
    if (!keyReserved(keys[m]))
      {
      dmcChar *result, *type=0;
      CHECK(dmcGetMetaValue(inp,keys[m],&result,&type));
      if (result)
        {
        par[keys[m]] = result;
        free(result);
        }
      else
        par[keys[m]] = "";

      if (type) free(type);
      }
    free(keys[m]);
    }
  free(keys);
  }

void noconv (PDT t1, PDT t2)
  {
  planck_fail (string("unsupported conversion from '")+type2string(t1)
               +string("' to '")+type2string(t2)+"'");
  }

} // unnamed namespace

iohandle_dmc::Manager::Manager (const string &properties)
  {
#ifdef USE_TOODI
// FIXME: this is not conforming to the CLLIO specification
  if (properties=="")
    CHECK(dmcInitializeLLIO("TOODI%file",&initHandle))
  else
    CHECK(dmcInitializeLLIO(properties.c_str(),&initHandle));
#else
  CHECK(dmcInitializeLLIO(properties.c_str(),&initHandle));
#endif

  CHECK(dmcOpenStores(initHandle,&session));
  CHECK(dmcBeginTransaction(session));
  if (properties!="")
    getSimParams (initHandle, paramMap);
  }
iohandle_dmc::Manager::Manager (int argc, const char **argv)
  {
  planck_assert(argc>=2,"incorrect command line format");
  if ((argc==2) && (string(argv[1]).find("=")==string::npos))
    {
    string properties=argv[1];
#ifdef USE_TOODI
// FIXME: this is not conforming to the CLLIO specification
    if (properties=="")
      CHECK(dmcInitializeLLIO("TOODI%file",&initHandle))
    else
      CHECK(dmcInitializeLLIO(properties.c_str(),&initHandle));
#else
    CHECK(dmcInitializeLLIO(properties.c_str(),&initHandle));
#endif
    CHECK(dmcOpenStores(initHandle,&session));
    CHECK(dmcBeginTransaction(session));
    if (properties!="")
      getSimParams (initHandle, paramMap);
    }
  else
    {
#ifdef USE_TOODI
// FIXME: this is not conforming to the CLLIO specification
    CHECK(dmcInitializeLLIO("TOODI%file",&initHandle))
#else
    CHECK(dmcInitializeLLIO("",&initHandle));
#endif
    parse_cmdline_equalsign(argc,argv,paramMap);
    CHECK(dmcOpenStores(initHandle,&session));
    CHECK(dmcBeginTransaction(session));
    }
  }
iohandle_dmc::Manager::~Manager ()
  {
  CHECK(dmcCommitTransaction(session));
  CHECK(dmcCloseStores(session));
  CHECK(dmcCloseLLIO(initHandle));
  }

paramfile iohandle_dmc::Manager::getParams (bool verbose) const
  { return paramfile(paramMap,verbose); }

void iohandle_dmc::assert_connected() const
  { planck_assert (connected, "not connected to an object"); }
void iohandle_dmc::assert_write() const
  {
  assert_connected();
  planck_assert (!readonly, "trying to write to a read-only object");
  }
void iohandle_dmc::assert_read() const
  {
  assert_connected();
  planck_assert (readonly, "reading not allowed from this object");
  }
void iohandle_dmc::clean_all()
  {
  if (!connected) return;
  CHECK(dmcCloseObject(HANDLE));
  if (!readonly) CHECK(dmcCommitAndResume(session));
  connected=false;
  }

//static
iohandle *iohandle_dmc::makeHandle()
  { return new iohandle_dmc; }

//static
string iohandle_dmc::Type()
#ifdef USE_TOODI
  { return "TOODI"; }
#else
  { return "HFIDMC"; }
#endif

iohandle_dmc::iohandle_dmc ()
  : connected(false), handle(new ObjectHandle) {}
// virtual
iohandle_dmc::~iohandle_dmc()
  { clean_all(); delete PHANDLE; }

void iohandle_dmc::fill_cache()
  {
  assert_connected();
  dmcChar *result;
#ifdef USE_TOODI
// FIXME: this is not conforming to the CLLIO specification
  type_=getStringKey_connected("objType");
#else
  CHECK(dmcGetObjectType(HANDLE,&result));
  type_ = result;
  free(result);
#endif
  int ncols;
  CHECK(dmcGetNumberOfColumns (HANDLE,&ncols));
  colname_.alloc(ncols);
  colunit_.alloc(ncols);
  coltype_.alloc(ncols);
  for (int m=0; m<ncols; ++m)
    {
    CHECK(dmcGetNameOfColumn (HANDLE,m,&result));
    colname_[m] = result;
    free(result);
    CHECK(dmcGetUnitOfColumn (HANDLE,m,&result));
    colunit_[m] = result;
    free(result);
    CHECK(dmcGetTypeOfColumn (HANDLE,m,&result));
    string stype (result);
    free (result);
    if (stype=="float64")
      coltype_[m] = PLANCK_FLOAT64;
    else if (stype=="float32")
      coltype_[m] = PLANCK_FLOAT32;
    else if (stype=="byte")
      coltype_[m] = PLANCK_INT8;
    else if (stype=="char")
      coltype_[m] = PLANCK_INT8;
    else if (stype=="int16")
      coltype_[m] = PLANCK_INT16;
    else if (stype=="int32")
      coltype_[m] = PLANCK_INT32;
    else if (stype=="int64")
      coltype_[m] = PLANCK_INT64;
    else if (stype=="string")
      coltype_[m] = PLANCK_STRING;
    else if (stype=="boolean")
      coltype_[m] = PLANCK_BOOL;
    else
      planck_fail ("Unsupported column type '"+stype+"'");
    }
  }

// virtual
void iohandle_dmc::create (const string &name, const string &type)
  {
  clean_all();

#ifdef USE_HFIDMC
// FIXME: this is not conforming to the CLLIO specification
  if (name=="UNKNOWN")
    {
    dmcChar *outname;
    CHECK(dmcGetID (session,1,type.c_str(),name.c_str(),&outname));
    cerr << "OUT ID " << outname << endl;
    CHECK(dmcOpenNewObject (session,type.c_str(),outname,PHANDLE));
    free(outname);
    }
  else
    CHECK(dmcOpenNewObject (session,type.c_str(),name.c_str(),PHANDLE));
#else
  CHECK(dmcOpenNewObject (session,type.c_str(),name.c_str(),PHANDLE));
#endif
  connected = true;
  readonly = false;
  fill_cache();
  }
// virtual
void iohandle_dmc::open (const string &name)
  {
  clean_all();
  CHECK(dmcOpenPersistentObject (session, "common.GenericCore", name.c_str(), PHANDLE));
  connected = true;
  readonly=true;
  fill_cache();
  }
// virtual
void iohandle_dmc::open (const string &name, const string &type)
  {
  clean_all();
  CHECK(dmcOpenPersistentObject (session,type.c_str(),name.c_str(),PHANDLE));
  connected = true;
  readonly=true;
  fill_cache();
  }
// virtual
void iohandle_dmc::close()
  { clean_all(); }

// virtual
void iohandle_dmc::flush()
  { CHECK(dmcCommitAndResume(session)); }

#ifdef USE_TOODI
//virtual
string iohandle_dmc::getID (const string &name, const string &type) const
  {
  char *result=0;
  CHECK(dmcGetID(session,true,type.c_str(),name.c_str(),&result));
  planck_assert(result,"error in dmcGetID()");
  string res = result;
  free(result);
  return res;
  }
#else
//virtual
string iohandle_dmc::getID (const string &, const string &) const
  {
  planck_fail ("(HFIDMC) iohandle_dmc::getID() called");
  }
#endif

// virtual
void iohandle_dmc::deleteObject(const string &name) const
  {
  ObjectHandle hnd;
  CHECK(dmcOpenPersistentObject (session,"common.GenericCore",name.c_str(),&hnd));
  CHECK(dmcDeletePersistent(hnd));
  CHECK(dmcCloseObject(hnd));
  CHECK(dmcCommitAndResume(session));
  }

// virtual
string iohandle_dmc::objectName() const
  {
  assert_connected();
#ifdef USE_TOODI
  return getStringKey_connected("objName");
#else
  return "unknown";
#endif
  }
// virtual
string iohandle_dmc::objectType() const
  {
  assert_connected();
  return type_;
  }
// virtual
int iohandle_dmc::numColumns() const
  {
  assert_connected();
  return colname_.size();
  }

// virtual
bool iohandle_dmc::keyPresent (const string &key) const
  {
  assert_read();
  dmcBoolean result;
  CHECK(dmcMetaKeyIsPresent(HANDLE,key.c_str(),&result));
  return result;
  }

//virtual
void iohandle_dmc::setKeyVoid (const string &key, const void *value, PDT type)
  {
  assert_write();

  if (keyReserved(key))
    {
    cerr << "Warning: key '"+key+"' is reserved and will not be written!"
         << endl;
    return;
    }

  switch (type)
    {
    case PLANCK_INT8:
      CHECK(dmcSetMetaValueAsByte(HANDLE, key.c_str(), *static_cast<const int8 *>(value)))
      break;
    case PLANCK_INT16:
      CHECK(dmcSetMetaValueAsInt16(HANDLE, key.c_str(), *static_cast<const int16 *>(value)))
      break;
    case PLANCK_INT32:
      CHECK(dmcSetMetaValueAsInt32(HANDLE, key.c_str(), *static_cast<const int32 *>(value)))
      break;
    case PLANCK_INT64:
      CHECK(dmcSetMetaValueAsInt64(HANDLE, key.c_str(), *static_cast<const int64 *>(value)))
      break;
    case PLANCK_FLOAT32:
      CHECK(dmcSetMetaValueAsFloat(HANDLE, key.c_str(), *static_cast<const float32 *>(value)))
      break;
    case PLANCK_FLOAT64:
      CHECK(dmcSetMetaValueAsDouble(HANDLE, key.c_str(), *static_cast<const float64 *>(value)))
      break;
    case PLANCK_BOOL:
      CHECK(dmcSetMetaValueAsBoolean(HANDLE, key.c_str(), *static_cast<const bool *>(value)))
      break;
    case PLANCK_STRING:
// FIXME: dmcSetMetaValueAsString() has a bad prototype, use dmcSetMetaValue()
// for now
      CHECK(dmcSetMetaValue(HANDLE, key.c_str(), static_cast<const string *>(value)->c_str(),"string"))
      break;
    default:
      planck_fail (string("unsupported data type '")+type2string(type)+"'");
    }
  }
//virtual
void iohandle_dmc::getKeyVoid (const string &key, void *value, PDT type) const
  {
  assert_read();
#ifdef USE_TOODI
//FIXME: should not be necessary and should be removed once the DMC fixes are in
  planck_assert(keyPresent(key),"Key '"+key+"' not present");
#endif
  switch (type)
    {
    case PLANCK_INT8:
      CHECK(dmcGetMetaValueAsByte(HANDLE, key.c_str(), static_cast<dmcByte *>(value)))
      break;
    case PLANCK_INT16:
      CHECK(dmcGetMetaValueAsInt16(HANDLE, key.c_str(), static_cast<dmcInt16 *>(value)))
      break;
    case PLANCK_INT32:
      CHECK(dmcGetMetaValueAsInt32(HANDLE, key.c_str(), static_cast<dmcInt32 *>(value)))
      break;
    case PLANCK_INT64:
      CHECK(dmcGetMetaValueAsInt64(HANDLE, key.c_str(), static_cast<dmcInt64 *>(value)))
      break;
    case PLANCK_FLOAT32:
      CHECK(dmcGetMetaValueAsFloat(HANDLE, key.c_str(), static_cast<dmcFloat *>(value)))
      break;
    case PLANCK_FLOAT64:
      CHECK(dmcGetMetaValueAsDouble(HANDLE, key.c_str(), static_cast<dmcDouble *>(value)))
      break;
    case PLANCK_BOOL:
      {
      dmcBoolean tmp;
      CHECK(dmcGetMetaValueAsBoolean(HANDLE, key.c_str(), &tmp))
      *static_cast<bool *>(value) = tmp;
      break;
      }
    case PLANCK_STRING:
      {
      dmcChar *tmp;
      CHECK(dmcGetMetaValueAsString(HANDLE, key.c_str(), &tmp))
      *static_cast<string *>(value) = tmp;
      free(tmp);
      break;
      }
    default:
      planck_fail (string("unsupported data type '")+type2string(type)+"'");
    }
  }

string iohandle_dmc::getStringKey_connected (const string &key) const
  {
  assert_connected();
  dmcChar *result, *type=0;
  CHECK(dmcGetMetaValue(HANDLE,key.c_str(),&result,&type));
  string result2 = "";
  if (result)
    {
    result2 = result;
    free(result);
    }
  if (type) free(type);
  return result2;
  }

//virtual
void iohandle_dmc::setKeyTypeless (const string &key, const string &value)
  {
  assert_write();
  if (keyReserved(key))
    {
    cerr << "Warning: key '"+key+"' is reserved and will not be written!"
         << endl;
    return;
    }
#ifdef USE_TOODI
  CHECK(dmcSetMetaValue(HANDLE, key.c_str(), value.c_str(),"voodoo"))
#else
  CHECK(dmcSetMetaValue(HANDLE, key.c_str(), value.c_str(),"string"))
#endif
  }
//virtual
string iohandle_dmc::getKeyTypeless (const string &key) const
  {
  assert_read();
  dmcChar *result, *type=0;
#ifdef USE_TOODI
//FIXME: should not be necessary and should be removed once the DMC fixes are in
  planck_assert(keyPresent(key),"Key '"+key+"' not present");
#endif
  CHECK(dmcGetMetaValue(HANDLE,key.c_str(),&result,&type));
  string result2 = "";
  if (result)
    {
    result2 = result;
    free(result);
    }
  if (type) free(type);
  return result2;
  }

// virtual
void iohandle_dmc::deleteKey (const string &key)
  {
  assert_write();
  CHECK(dmcRemoveMetaKey(HANDLE,key.c_str()));
  }

//virtual
void iohandle_dmc::getAllKeys (vector<string> &keys) const
  {
  assert_read();
  int nkeys;
  dmcChar **rawkeys;
  CHECK(dmcGetMetaKeys(HANDLE, &nkeys, &rawkeys));
  for (int m=0; m<nkeys; ++m)
    {
    if (!keyReserved(rawkeys[m]))
      keys.push_back(rawkeys[m]);
    free(rawkeys[m]);
    }
  free(rawkeys);
  }

//virtual
void iohandle_dmc::query (const string &query_string,
  vector<string> &result) const
  {
  int nresults;
  dmcChar **res;
  result.clear();
  CHECK(dmcQueryFor(session, query_string.c_str(), &nresults, &res));
  for (int m=0; m<nresults; ++m)
    {
    result.push_back(res[m]);
    free(res[m]);
    }
  free(res);
  }

//virtual
bool iohandle_dmc::columnPresent (const string &colname) const
  {
  assert_connected();
  for (int m=0; m<numColumns(); ++m)
    if (colname_[m]==colname) return true;
  return false;
  }

//virtual
int iohandle_dmc::columnNumber (const string &colname) const
  {
  assert_connected();
  for (int m=0; m<numColumns(); ++m)
    if (colname_[m]==colname) return m;
  planck_fail ("Column name '" + colname + "' not found.");
  }
// virtual
string iohandle_dmc::columnName (int colnum) const
  {
  assert_connected();
  return colname_[colnum];
  }
// virtual
PDT iohandle_dmc::columnType (int colnum) const
  {
  assert_connected();
  return coltype_[colnum];
  }
// virtual
string iohandle_dmc::columnUnit (int colnum) const
  {
  assert_connected();
  return colunit_[colnum];
  }
//  virtual
uint64 iohandle_dmc::columnLength (int colnum) const
  {
  assert_read();
  dmcInt64 result;
  CHECK(dmcGetSizeOfColumn (HANDLE,colnum,&result));
  return result;
  }
//  virtual
uint64 iohandle_dmc::efficientChunkSize (int /*colnum*/) const
  {
  assert_connected();
#ifdef USE_HFIDMC
  return 1000000000000ull; // avoid chunks at all costs for HFI
#else
  return 1024*1024; // FIXME: no better heuristics available at the moment
#endif
  }

// virtual
void iohandle_dmc::readColumnRawVoid
  (int colnum, void *data, PDT type, long num, uint64 offset) const
  {
  assert_read();
  switch (type)
    {
    case (PLANCK_FLOAT32):
      {
      dmcFloat *tdata = static_cast<dmcFloat *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_FLOAT32:
          CHECK(dmcGetFloatValues(HANDLE, colnum, offset, num, tdata));
          break;
        case PLANCK_FLOAT64:
          {
          arr<dmcDouble> tarr(num);
          CHECK(dmcGetDoubleValues(HANDLE, colnum, offset, num, &tarr[0]));
          for (int m=0; m<num; ++m) tdata[m]=tarr[m];
          break;
          }
        case PLANCK_INT32:
          {
          arr<dmcInt32> tarr(num);
          CHECK(dmcGetInt32Values(HANDLE, colnum, offset, num, &tarr[0]));
          for (int m=0; m<num; ++m) tdata[m]=tarr[m];
          break;
          }
        case PLANCK_INT64:
          {
          arr<dmcInt64> tarr(num);
          CHECK(dmcGetInt64Values(HANDLE, colnum, offset, num, &tarr[0]));
          for (int m=0; m<num; ++m) tdata[m]=tarr[m];
          break;
          }
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_FLOAT64):
      {
      dmcDouble *tdata = static_cast<dmcDouble *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_FLOAT64:
          CHECK(dmcGetDoubleValues(HANDLE, colnum, offset, num, tdata));
          break;
        case PLANCK_FLOAT32:
          {
          arr<dmcFloat> tarr(num);
          CHECK(dmcGetFloatValues(HANDLE, colnum, offset, num, &tarr[0]));
          for (int m=0; m<num; ++m) tdata[m]=tarr[m];
          break;
          }
        case PLANCK_INT32:
          {
          arr<dmcInt32> tarr(num);
          CHECK(dmcGetInt32Values(HANDLE, colnum, offset, num, &tarr[0]));
          for (int m=0; m<num; ++m) tdata[m]=tarr[m];
          break;
          }
        case PLANCK_INT64:
          {
          arr<dmcInt64> tarr(num);
          CHECK(dmcGetInt64Values(HANDLE, colnum, offset, num, &tarr[0]));
          for (int m=0; m<num; ++m) tdata[m]=tarr[m];
          break;
          }
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_INT8):
    case (PLANCK_UINT8):
      {
      dmcByte *tdata = static_cast<dmcByte *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT8:
        case PLANCK_UINT8:
          CHECK(dmcGetByteValues(HANDLE, colnum, offset, num, tdata));
          break;
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_INT16):
      {
      dmcInt16 *tdata = static_cast<dmcInt16 *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT16:
          CHECK(dmcGetInt16Values(HANDLE, colnum, offset, num, tdata));
          break;
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_INT32):
      {
      dmcInt32 *tdata = static_cast<dmcInt32 *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT32:
          CHECK(dmcGetInt32Values(HANDLE, colnum, offset, num, tdata));
          break;
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_INT64):
      {
      dmcInt64 *tdata = static_cast<dmcInt64 *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT64:
          CHECK(dmcGetInt64Values(HANDLE, colnum, offset, num, tdata));
          break;
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_BOOL):
      {
      bool *tdata = static_cast<bool *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_BOOL:
          {
          arr<dmcBoolean> tarr(num);
          CHECK(dmcGetBooleanValues(HANDLE, colnum, offset, num, &tarr[0]));
          for (int i=0; i<num; ++i)
            tdata[i] = tarr[i];
          break;
          }
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    case (PLANCK_STRING):
      {
      string *tdata = static_cast<string *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_STRING:
          {
          arr<dmcChar *> tdata2(num);
          CHECK(dmcGetStringValues(HANDLE,colnum,offset,num,&tdata2[0]));
          for (int m=0; m<num; ++m)
            {
            tdata[m] = tdata2[m];
            free (tdata2[m]);
            }
          break;
          }
        default:
          noconv (columnType(colnum),type);
        }
      break;
      }
    default:
      planck_fail (string("unsupported type '")+type2string(type)+"'");
    }
  }

// virtual
void iohandle_dmc::writeColumnRawVoid
  (int /*colnum*/, const void * /*data*/, PDT /*type*/, long /*num*/,
   uint64 /*offset*/)
  {
  planck_fail ("not yet implemented");
  }

// virtual
void iohandle_dmc::appendColumnRawVoid
  (int colnum, const void *data, PDT type, long num)
  {
  assert_write();
  switch (type)
    {
    case (PLANCK_FLOAT32):
      {
      const dmcFloat *tdata = static_cast<const dmcFloat *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_FLOAT32:
          CHECK(dmcAppendFloatValues(HANDLE, colnum, num, tdata, -1));
          break;
        case PLANCK_FLOAT64:
          {
          arr<dmcDouble> tarr(num);
          for (int m=0; m<num; ++m) tarr[m]=tdata[m];
          CHECK(dmcAppendDoubleValues(HANDLE, colnum, num, &tarr[0], -1));
          break;
          }
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_FLOAT64):
      {
      const dmcDouble *tdata = static_cast<const dmcDouble *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_FLOAT64:
          CHECK(dmcAppendDoubleValues(HANDLE, colnum, num, tdata, -1));
          break;
        case PLANCK_FLOAT32:
          {
          arr<dmcFloat> tarr(num);
          for (int m=0; m<num; ++m) tarr[m]=tdata[m];
          CHECK(dmcAppendFloatValues(HANDLE, colnum, num, &tarr[0], -1));
          break;
          }
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_INT8):
    case (PLANCK_UINT8):
      {
      const dmcByte *tdata = static_cast<const dmcByte *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT8:
        case PLANCK_UINT8:
          CHECK(dmcAppendByteValues(HANDLE, colnum, num, tdata, -1));
          break;
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_INT16):
      {
      const dmcInt16 *tdata = static_cast<const dmcInt16 *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT16:
          CHECK(dmcAppendInt16Values(HANDLE, colnum, num, tdata, -1));
          break;
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_INT32):
      {
      const dmcInt32 *tdata = static_cast<const dmcInt32 *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT32:
          CHECK(dmcAppendInt32Values(HANDLE, colnum, num, tdata, -1));
          break;
        case PLANCK_FLOAT32:
          {
          arr<dmcFloat> tarr(num);
          for (int m=0; m<num; ++m) tarr[m]=tdata[m];
          CHECK(dmcAppendFloatValues(HANDLE, colnum, num, &tarr[0], -1));
          break;
          }
        case PLANCK_FLOAT64:
          {
          arr<dmcDouble> tarr(num);
          for (int m=0; m<num; ++m) tarr[m]=tdata[m];
          CHECK(dmcAppendDoubleValues(HANDLE, colnum, num, &tarr[0], -1));
          break;
          }
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_INT64):
      {
      const dmcInt64 *tdata = static_cast<const dmcInt64 *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_INT64:
          CHECK(dmcAppendInt64Values(HANDLE, colnum, num, tdata, -1));
          break;
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_BOOL):
      {
      const bool *tdata = static_cast<const bool *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_BOOL:
          {
          arr<dmcBoolean> tarr(num);
          for (int i=0; i<num; ++i)
            tarr[i] = tdata[i];
          CHECK(dmcAppendBooleanValues(HANDLE, colnum, num, &tarr[0], -1));
          break;
          }
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    case (PLANCK_STRING):
      {
      const string *tdata = static_cast<const string *>(data);

      switch(columnType(colnum))
        {
        case PLANCK_STRING:
          {
          arr<dmcChar *> tdata2(num);
          for (int m=0; m<num; ++m)
            {
            tdata2[m] = new dmcChar[tdata[m].length()+1];
            strcpy (tdata2[m], tdata[m].c_str());
            }
          CHECK(dmcAppendStringValues(HANDLE,colnum,num,const_cast<const dmcChar **>(&tdata2[0]),-1));
          for (int m=0; m<num; ++m)
            delete[] tdata2[m];
          break;
          }
        default:
          noconv (type,columnType(colnum));
        }
      break;
      }
    default:
      planck_fail (string("unsupported type '")+type2string(type)+"'");
    }
// don't do this for TOODI to improve performance
#ifndef USE_TOODI
  CHECK(dmcCommitAndResume(session));
#endif
  }

namespace {

class DMC_Typecheck
  {
  public:
    DMC_Typecheck()
      {
      planck_assert (sizeof(dmcByte)  ==1, "bad DMC type size for dmcByte");
      planck_assert (sizeof(dmcInt16) ==2, "bad DMC type size for dmcInt16");
      planck_assert (sizeof(dmcInt32) ==4, "bad DMC type size for dmcInt32");
      planck_assert (sizeof(dmcInt64) ==8, "bad DMC type size for dmcInt64");
      planck_assert (sizeof(dmcFloat) ==4, "bad DMC type size for dmcFloat");
      planck_assert (sizeof(dmcDouble)==8, "bad DMC type size for dmcDouble");
      }
  };

DMC_Typecheck dmc_typechecker;

} // unnamed namespace

#endif
