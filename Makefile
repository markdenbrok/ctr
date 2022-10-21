SHELL = /bin/sh
LIBS =-lm  -L/z/mark/lib/ -lgsl -L/z/mark/lib/ -lgslcblas -lm  -L/z/mark/lib/ -lfftw3  -lcfitsio
PLOTLIBS = -L/usr/X11R6/lib -lX11 -lcpgplot -lpgplot -lpng -lz 
INCLUDEDIRS = -I ./
LIBDIRS =  
CC=gcc -g  -O0 -fno-inline -Wall $(INCLUDEDIRS)

.c.o: $(CC) -c -Wall $(INCLUDEDIRS) $<


CTR = kepler.o potential.o stuff.o line_model.o readcol.o read_mge.o\
	write_fits.o print_mge.o read_fits.o gaussfit.o calc_sigma_field.o

TEST_POT_MODEL = test_pot_model.o test_pot_model_2.o

TEST = test_tr_model.o


NGC3862_S = ngc3862_singlemodel.o

NGC3862_SER = ngc3862_series.o

NGC3862_MEM = ngc3862_memcee.o


ctr:	$(CTR)
	ar rcs libctr.a $(CTR)


ngc3862_singlemodel: $(NGC3862_S)
	make -s ctr
	$(CC) $(NGC3862_S) -o ngc3862_singlemodel $(LIBDIRS) -lctr $(LIBS) -L . -lpthread 

ngc3862_series: $(NGC3862_SER)
	make -s ctr
	$(CC) $(NGC3862_SER) -o ngc3862_series $(LIBDIRS) -lctr $(LIBS) -L . -lpthread 

ngc3862_memcee: $(NGC3862_MEM)
	make -s ctr
	$(CC) $(NGC3862_MEM) -o ngc3862_memcee $(LIBDIRS) -lctr -lmemcee $(LIBS) -L . -lpthread 

single_model: $(SINGLE_MODEL)
	make -s ctr
	$(CC) single_model.o -o single_model $(LIBDIRS) -lctr $(LIBS) -L . 

test:	$(TEST)
	make -s ctr
	$(CC) test_tr_model.o -o test_tr_model $(LIBDIRS) $(LIBS) -L . -lctr


test_pot_model:	$(TEST_POT_MODEL)
	make -s ctr
	$(CC) test_pot_model.o -o test_pot_model $(LIBDIRS) $(LIBS) -L . -lctr
	$(CC) test_pot_model_2.o -o test_pot_model_2 $(LIBDIRS) $(LIBS) -L . -lctr

clean: 
	rm *.o
	rm libctr.a
	rm nest/*

