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
   The sphere_m1.msh is the example of mesh data file, it is a sphere consists of two hemispheres.
   It was created by Gmsh geometry file sphere_m1.geo in mesh_sample folder.  

3. type './example1.out' with a argument of datafile name.  
   For example, '

