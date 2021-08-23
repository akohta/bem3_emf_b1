# bem3_emf_b1
This is the three-dimensional electromagnetic field analysis program for arbitrary object irradiated by arbitrary beams. 
This is based on boundary element method, the own developed numerical solution is used. 
This is the full vector field three-dimensional analysis, the corner probrem free. 
Intel Math Kernel Library is required. 
Gmsh is used to create a mesh data for object. 
The electromagnetic field analysis program "multi_fbeam" is used to analyze incident field. 

## Usage of example code

1. type 'make' command to compile.  
   The executable d3b1_bv_solver, example1.out, example2.out are created. 
   The executable d3b1_bv_solver is the main solver of boundary integral equations. 
   The example1.out is the executable of source code example1.c, it shows a simplest example of usage. 
   The example2.out is the execubable of source code example2.c, it shows a example of electromagnetic field intensity analysis. 
  
2. type './d3b1_bv_solver' with arguments of medium datafile name, mesh datafile name and output dafafile name.  
   For example, './d3b1_bv_solver medium_data.txt sphere_m1.msh ex.dat'. 
   The medium_data.txt is the sample of medium data file, two mediums are defined in it. The domain numbers are assigned to the medium from 1 in order. 
   The sphere_m1.msh is the example of mesh datafile, it is a sphere consists of two hemispheres.
   It was created by Gmsh geometry file sphere_m1.geo in mesh_sample folder. 
   The sphere_m1_image.pdf is the visulalization result of the sphere_m1.msh ( using Gmsh ).
   The d3b1_bv_solver solves boundary integral equations with specified data, outputs them to binary file with the output datafile name.
   The ipw.txt is the sample of incident field data file, a plane-wave is defined in it. Please refer to the "multi_fbeam" for detail. 
   This sample requires about 32GB of memory to solve boundary integral equations ( to run d3b1_bv_solver ). 
   The d3b1_bv_solver has optional arguments for rotation and translation of object. 
   When the vector defining rotation axis is (rx, ry, rz), the rotation angle is theta, the translation vector is (tx, ty, tz), 
   the arguments are './d3b1_bv_solver medium_data.txt sphere_m1.msh ex.dat rx ry rz theta tx ty tz'. 
   Rodrigues' rotation formula is used.  

3. type './example1.out' with a argument of datafile name outputed by d3b1_bv_solver.  
   For example, './example1.out ex.dat'. 
   This executable calculates electromagnetic field, radiation force and torque. 

4. type './example2.out' with a argument of datafile name outputed by d3b1_bv_solver.   
   For example, './example2.out ex.dat'. 
   This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.pdf is the visualization result of intensity distributions, created by the Gnuplot script gscript_example2.plt.
   
Please see 'd3b1_src/bem3_emf_b1.h' for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable 'OMP_NUM_THREADS'.  


## About mesh file

This code can use quadrangular ( bi-linear ) and triangular ( linear triangular ) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh data is in the folder mesh_sample. 
The file with extension '.geo' is the Gmsh geometry file. 
The file with extension '.msh' is the mesh file created from Gmsh geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 xxxx.geo' in command line ( xxxx.geo is a geometry file). 
The domain number ( Physical Surface ) 99 is assigned to the open region in Gmsh geometry file, becase Gmsh can't use the number 0 ( assigned to open region in the code). 
Please refre to the manual of Gmsh for detail of geometry file.  


## Verifications

### Verification 1  

The verification results using 'emf_mie_mmls' are in the folder 'verification1'.
This is the analysis result of plane wave scattering by the two-layered sphere.
The sphere_m2_image.pdf is the visualization result of mesh data. 
The I_example2.pdf is the visualization result of intensity distributions ( outputs of 'example2.out' ).
The results of 'emf_mie_mmls' is in the folder 'emf_mie_mmls_result'.

### Verification 2  

The verification results using 'emf_mie_ms' are in the folder 'verification2'.
This is the analysis result of plane wave scattering by three spheres.
The sphere_m3_image.pdf is the visualization result of mesh data. 
The I_example2.pdf is the visualization result of intensity distributions ( outputs of 'example2.out' ).
The results of 'emf_mie_ms' is in the folder 'emf_mie_ms_result'.


## Example of analysis

The additional example of analysis is in the folder 'analysis_sample'.
This is the analysis result of plane wave scattering by cone shaped metal.
The cone_m1_image.pdf is the visualization result of mesh data. 
The I_example2.pdf is the visualization result of intensity distributions ( outputs of 'example2.out' ), using Gnuplot script gscript_example2_logcb.plt.  


## References
1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)
2. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)
3. The electromagnetic field analysis program [multi_fbeam](https://github.com/akohta/multi_fbeam/) 
4. The electromagnetic field analysis program [emf_mie_mmls](https://github.com/akohta/emf_mie_mmls/)
5. The electromagnetic field analysis program [emf_mie_ms](https://github.com/akohta/emf_mie_ms/)
