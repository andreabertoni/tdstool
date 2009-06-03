TARGET = schro3D.x

#################################################################
#  stio3D_qtbm     time-independent Schroedinger solver with open
#  boundaries in 3D on a cuboidal domain.
#  Uses the Quantum Transmitting Boundary Method.
#  from stio2d_qtbm 3jun08 - 4jun08             by Andrea Bertoni
#################################################################
#   Makefile ifort 
#################################################################

## the command shell
SHELL = /bin/csh

## name of the fortran 90, fortran 77, C compilers
F90 = ifort 
F77 =
CC =

## flags to use for fortran 90, fortran 77, C compilation
F90FLAGS = -ipo -O2
F77FLAGS =
CFLAGS =

## name of the program to link the program units
LNK = $(F90)

## flags to use at link time
LFLAGS = $(F90FLAGS)

LLIBS = -L/usr/local/dislin/ -Llib -ldislin -lm -lfftw3

#################################################################
#      objects files
#################################################################

OBJS =  dislin.o            \
        mod_fftw3params.o   \
        mod_filewrite.o     \
        mod_bicubic.o       \
	mod_staticdata.o    \
	mod_strtopot.o      \
	mod_indata.o        \
	mod_splitstep.o     \
	mod_viewgraph.o     \
	mod_inputgui.o      \
	schro2.o

#################################################################
#      rules 
#################################################################

tdstool : $(OBJS)
	$(LNK) -o $(TARGET) $(LFLAGS) $^ $(LLIBS)

%.o : %.f90
	$(F90) $(F90FLAGS) -c $< -o $@


clean:
	rm -f $(OBJS)
	rm -f *.mod
	rm -f $(TARGET)

#################################################################
#      dependencies
#################################################################

schro.o : dislin.o          \
        mod_fftw3params.o   \
        mod_filewrite.o     \
        mod_bicubic.o       \
	mod_staticdata.o    \
	mod_strtopot.o      \
	mod_indata.o        \
	mod_splitstep.o     \
	mod_viewgraph.o     \
	mod_inputgui.o      \

mod_fftw3params.o: 

mod_filewrite.o:

mod_bicubic.o:

mod_staticdata.o:

mod_strtopot.o:		mod_staticdata.o

mod_indata.o:	        mod_staticdata.o      \
			mod_strtopot.o        \
			mod_bicubic.o

mod_mainalgo.o:	        mod_staticdata.o      \
			mod_indata.o          \
			mod_filewrite.o

mod_splitstep.o:	mod_staticdata.o      \
			mod_indata.o          \
			mod_filewrite.o       \
			mod_fftw3params.o

mod_viewgraph.o:	mod_staticdata.o      \
			mod_indata.o          \
			mod_filewrite.o       \
			dislin.o              \
			mod_bicubic.o

mod_inputgui.o:		mod_staticdata.o      \
			dislin.o              \
			mod_indata.o          \
			mod_strtopot.o        \
			mod_mainalgo.o        \
			mod_viewgraph.o
