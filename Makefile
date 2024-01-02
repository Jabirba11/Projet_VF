EXE = main
F90 = gfortran
#OPT = -O0 -pedantic -g -ffpe-trap=invalid,zero,overflow -fbounds-check -fcheck=all -Wall
OPT = -O2
OBJ = mod_parametres.o mod_fonctions.o mod_schemas.o mod_output.o mod_test.o $(EXE).o

$(EXE): $(OBJ)
	$(F90) $(OPT) -o $(EXE) $^

%.o: %.f90
	$(F90) $(OPT) -c $<

%.mod: %.f90
	$(F90) $(OPT) -c $<

clean:
	rm *.o *.mod $(EXE) *.vtk sortie.dat

exe: $(EXE)
	./$(EXE)

plot_error: error.gp error.dat
	gnuplot -p $<
	
error.dat: $(EXE) convergence.sh
	./convergence.sh $<