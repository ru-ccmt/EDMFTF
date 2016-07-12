include ../Makefile.in

all : gpoint.so yw_excor.so

gpoint.so : gpoint.f90
	$(CMP) $(F2PL) -c $? -m gpoint

yw_excor.so : yw_excor.f90
	$(CMP) $(F2PL) -c $? -m yw_excor


clean :
	rm -rf *.so *.o *.a *.mod *.pyc
	rm -rf *.dSYM
