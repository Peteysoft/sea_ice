BASEDIR = /home/petey

# processor type:
PROC = 64

LIBDIR = $(BASEDIR)/lib$(PROC)
INCLUDE_DIR = $(BASEDIR)/include

CC = g++
OPT_VER = #-g 
CFLAGS = -I$(INCLUDE_DIR) $(OPT_VER) -m$(PROC) -O2 
LDFLAGS = -L$(LIBDIR) -lgsl -lgslcblas

F77 = gfortran

mcc_ice_main: mcc_ice_main.cc mcc_ice_model.o mcc_ice_util.o inverf.o ckw.o bin_search.o
	$(CC) $(CFLAGS) -o mcc_ice_main mcc_ice_main.cc \
		mcc_ice_model.o mcc_ice_util.o inverf.o ckw.o bin_search.o $(LDFLAGS)

#ice_model_wrap.so: ice_model_wrap.cxx ice_model.o ice_model_util.o inverf.o ckw.o
#	$(CC) $(CFLAGS) -shared -o ice_model_wrap.so ice_model.o ice_model_util.o inverf.o ckw.o \
#		-L$(LIBDIR) -lpetey$(OPT_VER)$(PROC) -lgsl -lgslcblas

mcc_ice_idl_wrapper.so: mcc_ice_model.o mcc_ice_util.o inverf.o ckw.o bin_search.o mcc_ice_idl_wrapper.o
	ld -shared -o mcc_ice_idl_wrapper.so mcc_ice_idl_wrapper.o \
		mcc_ice_model.o mcc_ice_model_util.o inverf.o ckw.o bin_search.o \
		-L$(LIBDIR) -lpetey$(OPT_VER)$(PROC) -lgsl -lgslcblas

mcc_ice_idl_wrapper.o: mcc_ice_idl_wrapper.cc mcc_ice_model.h
	gcc $(CFLAGS) -c mcc_ice_idl_wrapper.cc

clean:
	rm -f *.o
	rm -f mcc_ice_main
	rm -f mcc_ice_idl_wrapper.so

mcc_ice_model.o: mcc_ice_model.cc mcc_ice_model.h
	$(CC) -c $(CFLAGS) mcc_ice_model.cc

mcc_ice_util.o: mcc_ice_util.cc mcc_ice_util.h
	$(CC) -c $(CFLAGS) mcc_ice_util.cc

bin_search.o: bin_search.cc
	$(CC) -c $(CFLAGS) bin_search.cc

inverf.o: inverf.c inverf.h
	$(CC) -c $(CFLAGS) inverf.c

ckw.o: ckw.f
	$(F77) -c $(CFLAGS) ckw.f

mcc_ice_model.h: mcc_ice_util.h

mcc_ice_util.h: bin_search.h

