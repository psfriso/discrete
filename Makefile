V = 0.2.7
EXEDIR = ../../exe
#
.PHONY: lib setup discrete
#
all: setup discrete
	
lib: 
	cd lib; make

discrete: lib
	cd discrete; make

setup:	lib
	cd setup; make

install: all
	cp discrete/discrete $(EXEDIR)/discrete$(V)
	cp setup/DMDSet $(EXEDIR)/DMDSetup$(V)

clean:
	cd lib;make clean
	cd discrete; make clean
	cd setup; make clean
	rm $(EXEDIR)/discrete$(V) $(EXEDIR)/DMDSet$(V)
        	
