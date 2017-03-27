include ../config.mk
#
OBJS = IO.o energy.o colisio.o DMD.o
MODS = paramSet.o Simulation.o random.o
EXEFILE = discrete
LIBFILE = $(LIBDIR)/libDiscrete.a

discrete: ${MODS} $(OBJS) $(LIBFILE) main.o
	$(F77) $(LFLAGS) $(ENDFLAG) -o $(EXEFILE) $(LOPTIONS) $(OBJS) $(MODS) main.o -lDiscrete 

main.o: vars.h	

.f.o: $< )
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
.F.o: $< 
	$(F77) -ffree-form $(OPTIONS) $(CFLAGS) $<
clean: 
	rm *.o *.mod $(EXEFILE)
