#test
MPICXX := $(HPCP_MPICXX)

CXXFLAGS := $(HPCP_CXXFLAGS) $(HPCP_OMPFLAGS)

INCLUDE := -I/project/projectdirs/planck/modules/hopper/gnu/level_s-20130510-20131216/include $(shell toastconfig --cppflags)

LINK := $(shell toastconfig --mpiflibs) $(HPCP_LDFLAGS) $(HPCP_LIBS) \
	-L/project/projectdirs/planck/modules/hopper/gnu/level_s-20130510-20131216/lib \
	-lcxxmod -lhealpix_cxx -lcxxsupport -lfftpack -lc_utils -lcfitsio \
	-L/project/projectdirs/planck/software/prezeau/reijo/IScalm.a

OBJECTS := cmult.o cmult_module.o

all : cmult

default : cmult

cmult_module.o: cmult_module.cc
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) -c $<

cmult.o: cmult_module.o cmult.cc
	$(MPICXX) $(CXXFLAGS) $(INCLUDE) -c cmult.cc

cmult : $(OBJECTS)
	$(MPICXX) -o cmult $(OBJECTS) $(LINK)

clean:
	rm -f *~ *.o

