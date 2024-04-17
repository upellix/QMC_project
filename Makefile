PROGRAM_NAME = exe
FC = gfortran

$(PROGRAM_NAME): code.F90 hamiltonian.F90 utility.F90 montecarlo.F90
	$(FC) -o $@ $^

clean:
	rm -f $(PROGRAM_NAME)

