
FC = gfortran
CC = gcc
RM = rm -f

FFLAGS = -fPIC --static --shared -O3
#FFLAGS = -g


OBJS = main.o sdbox.o problem.o
OBJSLIB = sdbox.o

all: sdbox_c lib

sdbox_c:  $(OBJS)
	$(CC) -o sdbox_c $(OBJS)  -lm
	mv sdbox_c ..

lib:  $(OBJSLIB)
	$(CC) -fPIC --shared -lm -o libsdbox.a $(OBJSLIB)
	mv libsdbox.a ../Julia_interface

.SUFFIXES : .c   .o

.c.o:   $* ; $(CC) $(FFLAGS) -c $*.c

clean: 
	$(RM) *.o
	$(RM) *.mod
	$(RM) sdbox_f
	$(RM) sdbox_c
	$(RM) libsdbox.a

