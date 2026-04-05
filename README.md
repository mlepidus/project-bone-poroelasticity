# Coupled Poroelastic Finite Element Solver

This code solves time‑dependent poroelasticity problems (e.g., fluid flow in deformable porous media) using a **monolithic finite element method**. It supports 2D/3D meshes, various boundary conditions, and outputs VTK files for visualisation.

**Key features**  
- Fully coupled displacement‑pressure formulation  
- Parallel solver support (MUMPS via `make parallel`)  
- Russian‑doll multiscale structure option  
- Convergence testing utilities

## Compilation

### Prerequisites
- [GetFEM](http://getfem.org/) – set its installation path before compiling:
  ```bash
  export GETFEM_PREFIX=/path/to/getfem
  ```

- Execute compilation:

```bash
make
```

### Additional options:
-j4: Compile in parallel using 4 processors

parallel: to compile with the necessary parallel libraries to have MUMPS solver

russian: to compile from the main_russian_doll, and have the PV-PLC structure. it has to be paired with an input file with the same structure as data_russian.txt

2>out: Redirect compilation output to a file named 'out' for detailed review

to print all the detailed current flags and active options:
```bash
make info
```

to print all the possible options and requirements:
```bash
make help
```

## Main Executables Description
main_coupled: Implements a fully coupled monolithic formulation of the problem, solving all equations within a single linear system.

(if we compile through make russian, the name will be main_coupled_russain)

## Usage
After successful compilation, execute the program with:

```bash
./main_coupled -f ./input/data_filename
```
replace data_filename with one of the following files:
### input files
paired with main_coupled:

- data_ideal_parameter: force applied on top, using all physical parameters as 1

- data_real_parameter: same as ideal parameter, but using real experimental physical parameters  

- data_shear: rotating stress applied to the principal axis 

- data_local_force: horizontal displacement applied only on a small portion of the bone  

- data_2d_annulus: 2d test with an annular mesh  

- data_3d_cube: 3d test with a unit cube  

- data_2d_square_convergence: this datafile has been used paired with the convergence_test.py to do the convergence test  


paired with main_russian_doll:
- data_russian: same as the data_ideal_parameter, but with the correct input structure for the Russian Doll model

## Output

- **Directory**: `output_vtk/` (created automatically if it doesn’t exist)  
- **Files**: One `.vtk` file per time step, containing pressure and displacement fields.  
- **Visualization**: Open with ParaView, VisIt, or any VTK‑compatible viewer.

## Requirements & Notes

- **Memory** – The monolithic solver can be memory‑intensive. For large 3D problems, 16+ GB RAM is recommended.  
- **Permissions** – Ensure you have write access to the `output_vtk/` directory.  
- **Paths** – Double‑check the `GETFEM_PREFIX` and any other library paths in the Makefile before compiling.
