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

2>out: Redirect compilation output to a file named 'out' for detailed review


# Main Executables Description
main_coupled: Implements a fully coupled monolithic formulation of the problem, solving all equations within a single linear system.

## Execution Instructions
After successful compilation, execute the program with:

```bash
./main -f data_filename
```
## data files
here are the different data_filename that can be passed as input
-data_ideal_parameter: force applied on top, using all physical parameters as 1
-data_real_parameter: same as ideal parameter, but using real experimental physical parameters
-data_shear: rotating stress applied to the principal axe
-data_local_force: horizonal displacement applied only on a small portion of the bone


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

