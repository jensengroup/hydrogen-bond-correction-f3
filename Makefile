
All:
	@gfortran -c src/f3.f90 -o f3.o
	@gfortran -c src/f3_program.f90 -o f3_program.o
	@gfortran f3.o f3_program.o -o f3_exe

clean:
	rm *.o
	rm *.mod
