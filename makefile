
FC = gfortran
CC = gcc
RM = rm -f

#FFLAGS = -O3
FFLAGS = -g


all: fver cver

fver:
	(cd ./F90_src ; make ; cd ..)

cver:
	(cd ./C_src ; make ; cd ..)

clean: 
	(cd ./C_src ; make clean ; cd ..)
	(cd ./F90_src ; make clean ; cd ..)
	$(RM) sdbox_f
	$(RM) sdbox_c
	$(RM) PYTHON_src/*.pyc
	$(RM) libsdbox.a
	$(RM) Julia_interface/libsdbox.a

