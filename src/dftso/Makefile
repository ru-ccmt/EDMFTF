include ../Makefile.in

#FORT = gfortran
#PREPROC = cpp
#CC = cc
#FFLAGS     = -O3 -ffree-line-length-none #-fcheck=all # -DFFTW2 -FR -mp1 -w -prec_div -pc80 -pad -align -DINTEL_VML -traceback
#WLIBS   =  -framework Accelerate -lfftw3_omp -lfftw3 -lm -L/usr/local/lib #-mkl -L/usr/lib64 -lpthread -L/opt/intel/composer_xe_2013.1.117/mkl/lib/intel64/ -lfftw2xf_double_intel -lfftw2x_cdft_DOUBLE -lmkl_rt

FORT = $(PF90)
LDFLAGS= $(WFOPT) $(WFOPTL) 
FFLAGS  = $(WFOPT) 
###############################################################################

#..............................................................................
#
#  Routines which can be compiled unmodified
#
OBJS = reallocate.o modules.o abclm.o    couple.o   errclr.o  errflg.o   garadme.o gtfnam.o  \
        harmon2.o  atpar.o    hmsec.o   init.o     kptin.o   latgen2.o \
        outwin.o   rint13.o  lapwso.o   cali.o    sphbes.o vnsrint.o  \
        lomain.o   abc.o      cputim.o  dvbes1.o   bess.o    lagder.o  \
        rotate.o   ylm.o      harmon.o  sphbru.o   outerr.o  lap_bp.o \
        inouh.o    diracout.o inth.o    select.o   abc_r.o   rlomain.o \
        kptout.o  angle.o gaunt1.o gaunt2.o \
        dergl.o    hscalc.o   dergl2.o  vderiv.o  rotdef.o \
        symop.o garadorb.o hns.o hsocalc.o horbcalc.o vorblo.o euler.o

OBJSMOD = abc.o abc_r.o abclm.o angle.o atpar.o cali.o couple.o \
          dergl.o dergl2.o diracout.o garadme.o garadorb.o harmon.o \
          hns.o horbcalc.o hscalc.o hsocalc.o \
          init.o inouh.o kptin.o kptout.o lapwso.o latgen2.o lomain.o opnfs.o \
          outwin.o rint13.o rlomain.o rotdef.o stern.o vderiv.o vectest.o \
          vnsrint.o vorblo.o


#..............................................................................

lapwso : $(OBJS)
	$(FORT) -o lapwso $(OBJS) $(LDFLAGS) $(WLIBS)
#..............................................................................
#
#  All routines depend upon an include file (contains common PARAMETERS)
#
#$(OBJS):		Makefile
#$(OBJSMOD):		modules.o

clean:
	-rm -f *.o _tmp_.* *.P .real .complex *.mod *__genmod.f90 lapwso

#...................
# Pattern rules
#...................
%.o : %.f90
	$(FORT) $(FFLAGS) -c $< -o $@

%.o : %.FP
	$(preproc) -P  $(Fmpi_define) $< $*.f90
	$(FORT) $(FFLAGS) -c $*.f90
	rm $*.f90

%.o : %.c
	$(CC) -c $< -o $@

%.so : %.f90
	$(CMP) -c $< -m $*
