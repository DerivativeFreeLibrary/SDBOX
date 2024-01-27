
FC = gfortran
CC = gcc
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g


OBJS = sdbox.o problem.o

all: sdbox_f

sdbox_f:  $(OBJS)
	$(FC) -o sdbox_f $(OBJS)
	mv sdbox_f ..

.SUFFIXES : .f90 .o

.f90.o: $* ; $(FC) $(FFLAGS) -c $*.f90

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) sdbox_f
	$(RM) sdbox_c

