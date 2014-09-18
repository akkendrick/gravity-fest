# Makefile for GeoFEST Finite Element Code
#		Unified Version

  PROGRAM 	= GeoFEST

  CC 		= cc
# CC 		= gcc
# CFLAGS	= -g -Wall

  libs =-lm

  COMMONS=mydefs.h finel.h errcodes.h output_phase.h
  HEADERS=generat.h inphase.h main.h solver.h stiff.h strain.h utility.h grav.h quake.h
  OBJECTS=generat.o inphase.o main.o solver.o stiff.o strain.o utility.o grav.o quake.o
  SOURCES=generat.c inphase.c main.c solver.c stiff.c strain.c utility.c grav.c quake.c

all:	$(OBJECTS) $(PROGRAM)

clean:
	/bin/rm -f $(OBJECTS) $(PROGRAM)

$(PROGRAM):	$(OBJECTS)
	$(CC) -o $@ $(OBJECTS) $(libs)

$(OBJECTS):	$(SOURCES) $(HEADERS)	$(COMMONS)
	$(CC) -c $(SOURCES)


