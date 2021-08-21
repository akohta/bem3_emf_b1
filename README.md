# bem_emf_b1
This is the electromagnetic field analysis program for arbitrary object irradiated by arbitrary beams. This is based on boundary element method, the own developed numerical solution is used. This is the full vector field analysis, the corner probrem free. Intel Math Kernel Library is required. Gmsh is used to create a mesh data for object. The electromagnetic field analysis program "multi_fbeam" is used to analyze incident field. 

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

3. type './example1.out' with a argument of datafile name outputed by d3b1_bv_solver.  
   For example, './example1.out ex.dat'.  
   This executable calculates electromagnetic field, radiation force and torque. 

4. type './example2.out' with a argument of datafile name outputed by d3b1_bv_solver.   
   For example, './example2.out ex.dat'.   
   This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.pdf is the visualization result of intensity distributions, created by the Gnuplot script gscript_example2.plt.
   
