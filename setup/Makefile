# 
# Makefile for DMDsetup
#
include ../config.mk

OBJS = resLibrary.o potentials.o Structure.o IO.o paramSet.o
INCS = vars.h

#shorter exe name due to cygwin bug
EXEFILE = DMDSet
#EXEFILE = DMDSetup

DMDsetup: $(MODS) $(OBJS) $(LIBDIR) main.o
	$(F77) $(LFLAGS) $(ENDFLAG) -o $(EXEFILE) $(LOPTIONS) $(OBJS) main.o -lDiscrete 

main.o: vars.h

Structure.o: Structure.f 
	$(F77) -ffree-form $(OPTIONS) $(CFLAGSNOCHECK) Structure.f

.f.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
.F.o: $<
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<

clean: 
	rm *.o *.mod $(EXEFILE)
