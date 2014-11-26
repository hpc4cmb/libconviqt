PKG:=libfftpack

SD:=$(SRCROOT)/$(PKG)
OD:=$(BLDROOT)/$(PKG)

FULL_INCLUDE+= -I$(SD)

HDR_$(PKG):=$(SD)/*.h
LIB_$(PKG):=$(LIBDIR)/libfftpack.a
BIN:=ffttest
OBJ:=fftpack.o bluestein.o ls_fft.o
ALLOBJ:=$(OBJ) ffttest.o
OBJ:=$(OBJ:%=$(OD)/%)
ALLOBJ:=$(ALLOBJ:%=$(OD)/%)

ODEP:=$(HDR_$(PKG)) $(HDR_c_utils)
BDEP:=$(LIB_$(PKG)) $(LIB_c_utils)

$(OD)/ffttest.o: $(SD)/ls_fft.c
$(OD)/fftpack.o: $(SD)/fftpack_inc.c

$(ALLOBJ): $(ODEP) | $(OD)_mkdir
$(LIB_$(PKG)): $(OBJ)
BIN:=$(BIN:%=$(BINDIR)/%)
$(BIN): $(BINDIR)/% : $(OD)/%.o $(BDEP)

all_hdr+=$(HDR_$(PKG))
all_lib+=$(LIB_$(PKG))
all_cbin+=$(BIN)
