# QUANTUM MONTE CARLO Project

The ```code.F90``` is a program to compute energies for simple hydrogen or helium-based system using Monte Carlo algorithms. Specifically, the program implements the Variational Monte Carlo and the Pure Diffusion Monte Carlo algorithms.
The code is capable of computing energy for systems made by hydrogen and Helium with at least 2 electrons, such as H, He, H<sub>2</sub><sup>+</sup>, H<sub>2</sub> and H<sub>3</sub><sup>+</sup>.
For more details, please refer to ```Document``` file.

## Download and installation
To get a copy of the program, clone the repository from <https://github.com/upellix/QMC_project.git> or download the ```.zip``` file from GitHub.
The program is written in Fortran 90 so a fortran compiler (*gfortran*, *pgfortran*, *ifort*, ...) is required to compile the code.

## Run the program  
Afetr installation, the program can be compiled by executing:
```
make
```
Esure this command is run in the same directory of the ```Makefile```.

To run the code, two file are needed:
- ### Input file ```mc_input.inp```
  This file contains all the parameters for the calculation,allowing selection between 
the Variational Monte Carlo (var) algorith and the Pure Diffusion Monte Carlo (dif) 
algorithm. An example is provided in the repository.
- ### Geometry file ```[syestem].xyz```
  This file contain the geometry information of the system, including the number and the position of the atoms. Several exaples of tested system are available in the repository.

Once these files are prepared in the same directory as the program, execute it by running:
```
./exe
```

## Output
The program gives the **energy** and the **acceptance ratio** togheter with their respective errors directly to the terminal.
Additionally, it generates an output file ```energy.dat``` containing the energy per walker in a format suitable for plotting.

## License
The code is open source and is licensed under GPL License. Please refer to ```License``` for more information.
