# Compilation Instructions
To compile the code, please follow these steps:

1) Configure the Makefile:
Adapt the Makefile by specifying your local paths (particularly for GetFEM and other dependencies).

2) Execute compilation:

```bash
make
```

## Additional options:
-j4: Compile in parallel using 4 processors

parallel: to compile with the necessary parallel libraries to have MUMPS solver

russian: to compile from the main_russian_doll, and have the PV-PLC structure. it has to be paired with an input file with the same structure as data_russian.txt

2>out: Redirect compilation output to a file named 'out' for detailed review

to print all the detailed info execute

```bash
make info
```

# Main Executables Description
main_coupled: Implements a fully coupled monolithic formulation of the problem, solving all equations within a single linear system.

(if we compile through make russian, the name will be main_coupled_russain)

## Execution Instructions
After successful compilation, execute the program with:

```bash
./main_coupled -f ./input/data_filename
```
## data files
here are the different data_filename inside the input folder, that can be passed as input
-data_ideal_parameter: force applied on top, using all physical parameters as 1
-data_real_parameter: same as ideal parameter, but using real experimental physical parameters
-data_shear: rotating stress applied to the principal axe
-data_local_force: horizonal displacement applied only on a small portion of the bone
-data_2d_annulus: 2d test with a annular mesh
-data_3d_cube: 3d test with an unit cube
-data_2d_square_convergence: this datafile has been used paired with the convergence_test.py to do the convergence test

-data_russian: data to use with the main_russian_doll.cc file, in order to have the correct input structure for the PV-PLC problem

## Output Specifications
The program generates output in the following format:

Results are exported to the output_vtk directory

Output includes pressure and displacement fields at each time step

Corresponding exact solution values are also exported (can be ignored if no exact solution has been configured)

## Additional Notes
For optimal performance, ensure that:

- All dependency paths are correctly specified in the Makefile

- Sufficient memory is available for the monolithic solver approach

- Output directory permissions allow file creation and writing

This implementation provides a comprehensive framework for solving coupled poroelastic problems using a monolithic finite element approach. The VTK output format allows for convenient visualization of results using standard scientific visualization tools.

