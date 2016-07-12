# C++ = g++
# CFLAGS = -O3

CFLAGS = $(OFLAGS)
include ../../Makefile.in

all : broad

broad : broad.cc
	$(C++) $(CFLAGS) -o $@ $? $(LIBS)

clean :
	rm -f broad

