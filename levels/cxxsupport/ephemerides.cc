#include "ephemerides.h"
#include "lsconstants.h"

using namespace std;

#ifndef USE_HFIDMC

class ephemerides_cllio: public ephemerides
  {
  private:
    arr<string> body, scalar, array;
    arr<double> data;
    tsize nsamp;
    double t0, inv_dt;
    safe_ptr<iohandle> inp;
    tsize o_x, o_y, o_z, o_angdiam, o_dsun, o_sto;

    tsize getScalarOffset (const string &name) const
      { return scalar.find(name); }

    double getScalar (tsize offset) const
      { return data[offset]; }

    double getScalar (const string &name) const
      { return getScalar(getScalarOffset(name)); }

    tsize getArrayOffset (const string &name) const
      { return scalar.size() + nsamp*array.find(name); }

    struct inter_struct
      {
      double frac;
      tsize idx;
      };

    inter_struct getArrayInterpol (double time) const
      {
      time -= t0;
      planck_assert(time>=0.,"requested time < t0");
      inter_struct inter;
      inter.frac = time*inv_dt;
      inter.idx = tsize(inter.frac);
      planck_assert(inter.idx<nsamp-1,"requested time too large");
      inter.frac -= double(inter.idx);
      return inter;
      }

    double getArrayValue (tsize offset, const inter_struct &inter) const
      {
      tsize ofs = offset + inter.idx;
      return (1.-inter.frac)*data[ofs] + inter.frac*data[ofs+1];
      }
    double getArrayValue (tsize offset, double time) const
      { return getArrayValue(offset,getArrayInterpol(time)); }
    double getArrayValue (const string &name, double time) const
      { return getArrayValue (getArrayOffset(name), time); }

  public:
    ephemerides_cllio (const string &obj)
      {
      inp = HandleManager.openObject(obj,"ephemeris.LS_ephemeris");
      nsamp = inp->getKey<int64>("nsamples");
      t0 = inp->getKey<double>("T0");
      inv_dt = 1./inp->getKey<double>("deltaT");
      inp->readEntireColumn("bodies",body);
      inp->readEntireColumn("scalars",scalar);
      inp->readEntireColumn("arrays",array);

      data.alloc (scalar.size() + nsamp*array.size());

      planck_assert(data.size()>0,"no data fields in ephemerides file");
      planck_assert(body.size()>0,"no bodies in ephemerides file");

      o_x = getArrayOffset("x_rel_sat_ecl/m");
      o_y = getArrayOffset("y_rel_sat_ecl/m");
      o_z = getArrayOffset("z_rel_sat_ecl/m");
      o_angdiam = getArrayOffset("angdiam_rel_sat/arcsec");
      o_dsun = getArrayOffset("d_sun/m");
      o_sto = getArrayOffset("angle_S_T_O/deg");
      }

    virtual const arr<string> &getBodies () const
      { return body; }

    virtual void loadBody (const string &name)
      { inp->readColumn("data",data,body.find(name)*data.size()); }

    virtual vec3 posRelSat_m (double time) const
      {
      inter_struct inter = getArrayInterpol (time);
      return vec3 (getArrayValue (o_x, inter),
                   getArrayValue (o_y, inter),
                   getArrayValue (o_z, inter));
      }
    virtual double dSun_m (double time) const
      { return getArrayValue (o_dsun, time); }
    virtual double angdiam_arcmin (double time) const
      { return 1./60.*getArrayValue (o_angdiam, time); }
    virtual double angle_S_T_O_rad (double time) const
      { return degr2rad*getArrayValue (o_sto, time); }
  };

ephemerides *getEphemerides (const string &obj)
  { return new ephemerides_cllio (obj); }

#else

#include <iomanip>
#include "pointing.h"
#include "HL2_DMC/PIODB.h"
#include "HL2_COMMONLLIO/PIOCommonLLIO.h"

#define toCharP(arg) const_cast<char *>(arg)
#define CHECK(call) { int status = call; if (status!=0) { cerr << "Error while doing " << endl << #call << endl << "exit status: " << status << endl; planck_fail("Stopped");} }

class ephemerides_hfi: public ephemerides
  {
  private:
    PIOGroup* ephgrp;
    arr<string> body, scalar, array, dmcobjname;
    arr<double> R, Delta, STO, ObsEcLon, ObsEcLat, AngDiam, o_x, o_y, o_z;
    tsize nsamp;
    double t0, inv_dt;

    struct inter_struct
      {
      double frac;
      tsize idx;
      };

    inter_struct getArrayInterpol (double time) const
      {
      if (time > 1.5e9) time -= t0; // prevents unwanted errors when asking for times already t0-corrected (1.5e9 is somewhere in 2005 in secs)
      else cerr << "WARNING: getArrayInterpol.time = " << setprecision(12) << time << " / ephemerides.t0 = " << setprecision(12) << t0 << endl;

      if (time < 0.0) {
        cerr << "getArrayInterpol.time = " << setprecision(12) << time << " / ephemerides.t0 = " << setprecision(12) << t0 << endl;
        planck_assert(time >= 0.0, "requested time < t0");
      }
      inter_struct inter;
      inter.frac = time*inv_dt;
      inter.idx = tsize(inter.frac);
      planck_assert(inter.idx<nsamp-1,"requested time too large");
      inter.frac -= double(inter.idx);
      return inter;
      }

    double getArrayValue (const arr<double> &a, const inter_struct &inter) const
      {
      return (1.0-inter.frac)*a[inter.idx] + inter.frac*a[inter.idx+1];
      }

    double getArrayValue (const arr<double> &a, double time) const
      {
      return getArrayValue(a, getArrayInterpol(time));
      }

    void readObjectColumn(const string &objname, const string &colname,
      arr<double> &coldata)
      {
      int objidx = body.find(objname),
          colidx = array.find(colname);
      PIOSTRING piocommand;
      void *dmcdata;

      cout << "reading column " << array[colidx] << " of object "
           << body[objidx] << endl;
      sprintf(piocommand, "tab=%d:%d,%d:%ld", colidx, colidx, 0, nsamp-1);
      PIOErr pioerr = PIOReadTAB2DObject(&dmcdata,
        toCharP(dmcobjname[objidx].c_str()), toCharP("PIODOUBLE"), piocommand,
        ephgrp);
      if (pioerr<0) {
        cerr << "error #" << pioerr << " in PIOReadTAB2DObject" << endl;
        planck_fail("error in PIOReadTAB2DObject");
      }
      planck_assert(tsize(pioerr)==nsamp,
        "wrong number of samples in PIOReadTAB2DObject");
      memcpy(static_cast<void *>(&coldata[0]), dmcdata, nsamp*sizeof(double));
      PIODeleteTAB2DTable(&dmcdata, ephgrp);
      }

  public:
    ephemerides_hfi (const string &obj)
      {

      PIOSTRING GrpName, Comment, keyval;

// **** open phemerides TAB2D group
      // inp = HandleManager.openObject(obj,"ephemeris.LS_ephemeris");
      dmcGetObjectGroup(toCharP(obj.c_str()), GrpName);
      cout << "Opening HFI ephemerides group: " << GrpName << endl;

      ephgrp = static_cast<PIOGroup*> (PIOOpenTAB2DGrp(GrpName, toCharP("r")));
      planck_assert(ephgrp!=0, "Unable to open ephemerides group");

// **** read Time of first data sample in seconds since Jan 1, 1958
      // t0 = inp->getKey<double>("T0");
      CHECK(PIOReadKeywordGrp(&t0, Comment, toCharP("T0_UT"), toCharP("PIODOUBLE"), ephgrp));
      cout << "t0: " << t0 << endl;

// **** read data sampling period in seconds
      // dt = inp->getKey<double>("deltaT");
      double dt;
      CHECK(PIOReadKeywordGrp(&dt, Comment, toCharP("deltaT"), toCharP("PIODOUBLE"), ephgrp));
      //dt = 3600.0;
      cout << "dt: " << dt << endl;
      inv_dt = 1./dt;

// **** get list of SSOs
      // inp->readEntireColumn("bodies",body);
      PIOSTRING  *FLGname = NULL;
      PIOINT      NbFlg = 0;
      PIOSTRING  *T2DObjType = NULL;
      PIOSTRING  *T2DObjName = NULL;
      PIOINT      NbT2DObj = 0;
      PIOSTRING   tempname;

      CHECK(PIOInfoTAB2DGrp(&FLGname, &NbFlg, &T2DObjType, &T2DObjName, &NbT2DObj, ephgrp));
      body.alloc(NbT2DObj);
      dmcobjname.alloc(NbT2DObj);
      for (int i=0; i<NbT2DObj; i++) {
        dmcobjname[i] = string(T2DObjName[i]);
        CHECK(PIOReadKeywordObject(static_cast<void *>(tempname), Comment, toCharP("SSOname"), toCharP("PIOSTRING"), toCharP(dmcobjname[i].c_str()), ephgrp));
        body[i] = string(tempname);
        cout << "reading SSO[" << i << "]: " << body[i] << " in " << dmcobjname[i] << endl;
      }
      CHECK(PIOFreeInfoTAB2DGrp(FLGname, T2DObjType, T2DObjName));

// **** get column names ("array")
      // inp->readEntireColumn("arrays",array);
      PIOLONG   NbKeys, NbCols, colnumber;
      PIOSTRING *keywords, *keytypes;

      NbKeys = PIOKeywordListGrp(&keywords, &keytypes, ephgrp);
      planck_assert(NbKeys>0, "Unable to read keywords from ephemerides group");
      NbCols = 0;
      for (int i=0; i<NbKeys; i++) {
        if (strncmp(keywords[i], "col_", 4) == 0) {
          NbCols++;
        }
      }

      array.alloc (NbCols);
      for (int i=0; i<NbKeys; i++) {
        if (strncmp(keywords[i], "col_", 4) == 0) {
          CHECK(PIOReadKeywordGrp(keyval, Comment, keywords[i], keytypes[i], ephgrp));
          colnumber = atol(keywords[i]+4) - 1;
          array[colnumber] = string(keyval);
          cout << "col[" << colnumber << "]: " << array[colnumber] << endl;
        }
      }

      _PIOFREE(keywords);
      _PIOFREE(keytypes);

// **** get number of samples per object
      // nsamp = inp->getKey<int64>("nsamples");
      PIOLONG begidx, endidx;

      begidx = PIOGetBeginObjectIdx(toCharP(dmcobjname[0].c_str()), ephgrp);
      cout << "begidx: " << begidx << endl;
      endidx = PIOGetEndObjectIdx(toCharP(dmcobjname[0].c_str()), ephgrp);
      cout << "endidx: " << endidx << endl;
      nsamp = (endidx - begidx + 1) / array.size();
      cout << "nsamp: " << nsamp << endl;

      // inp->readEntireColumn("scalars",scalar);
      R.alloc       (nsamp);
      Delta.alloc   (nsamp);
      STO.alloc     (nsamp);
      ObsEcLon.alloc(nsamp);
      ObsEcLat.alloc(nsamp);
      AngDiam.alloc (nsamp);
      o_x.alloc     (nsamp);
      o_y.alloc     (nsamp);
      o_z.alloc     (nsamp);

      planck_assert(array.size()>0,"no data fields in ephemerides file");
      planck_assert(body.size()>0,"no bodies in ephemerides file");
      }

    virtual ~ephemerides_hfi ()
      {
      planck_assert(PIOCloseTAB2DGrp(&ephgrp) == 0, "PIOCloseTAB2DGrp");
      }

    virtual const arr<string> &getBodies () const
      { return body; }

    virtual void loadBody (const string &name)
      {
      readObjectColumn(name, "r(au)",            R);
      readObjectColumn(name, "delta(au)",        Delta);
      readObjectColumn(name, "S-T-O(deg)",       STO);
      readObjectColumn(name, "ObsEcLon(deg)",    ObsEcLon);
      readObjectColumn(name, "ObsEcLat(deg)",    ObsEcLat);
      readObjectColumn(name, "Ang-diam(arcsec)", AngDiam);

//      pointing ptg;
      vec3 relpos;
      for (tsize i = 0; i < nsamp; i++) {
        relpos = pointing(halfpi - ObsEcLat[i] * degr2rad, ObsEcLon[i] * degr2rad);
        relpos *= astronomicalUnit * Delta[i];
        o_x[i] = relpos.x;
        o_y[i] = relpos.y;
        o_z[i] = relpos.z;
      }
      }

    virtual vec3 posRelSat_m (double time) const
      {
      if (0) {
        pointing ptg;
        inter_struct inter = getArrayInterpol (time);

        ptg.phi = getArrayValue (ObsEcLon, inter) * degr2rad;
        ptg.theta = halfpi - getArrayValue (ObsEcLat, inter) * degr2rad;
        vec3 relpos = ptg; // implicit conversion to x,y,z;
        relpos *= astronomicalUnit * getArrayValue (Delta, inter);
        return relpos;
      }

      inter_struct inter = getArrayInterpol (time);
      return vec3 (getArrayValue (o_x, inter),
                   getArrayValue (o_y, inter),
                   getArrayValue (o_z, inter));

      }

    virtual double dSun_m (double time) const
      { return astronomicalUnit * getArrayValue (R, time); }

    virtual double angdiam_arcmin (double time) const
      { return 1.0/60.0 * getArrayValue (AngDiam, time); }

    virtual double angle_S_T_O_rad (double time) const
      { return degr2rad * getArrayValue (STO, time); }
  };

ephemerides *getEphemerides (const string &obj)
  { return new ephemerides_hfi (obj); }

#endif
