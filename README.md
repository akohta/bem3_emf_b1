# bem3_emf_b1
This is the three-dimensional electromagnetic field analysis program for arbitrary object irradiated by arbitrary beams. 
This is based on boundary element method, the own developed numerical solution is used. 
This is the full vector field three-dimensional analysis, the corner probrem free. 
Intel Math Kernel Library and libpng are required. 
Gmsh is used for create a mesh data of object. 
The electromagnetic field analysis program "multi_fbeam" is used for analyze incident field. 


## Usage of example code  

1. type 'make' command to compile.  
   The executable d3b1_bv_solver, example1.out, example2.out, example3.out are created. 
   The executable d3b1_bv_solver is the main solver of boundary integral equations. 
   The example1.out is the executable of source code example1.c, it shows a simplest example using "bem3_emf_b1". 
   The example2.out is the execubable of source code example2.c, it shows a example of electromagnetic field intensity analysis.
   The example3.out is the executable of source code example3.c, it shows a example of outputting the instantaneous value of electromagnetic field as an image. 
  
2. type './d3b1_bv_solver' with arguments of medium datafile name, mesh datafile name and output dafafile name.  
   For example, './d3b1_bv_solver medium_data.txt sphere_m1.msh ex.dat'. 
   The medium_data.txt is the sample of medium datafile, two mediums are defined in it. The domain numbers are assigned to the medium from 1 in order. 
   The sphere_m1.msh is the example of mesh datafile, it is a sphere consists of two hemispheres.
   It was created by Gmsh geometry file sphere_m1.geo in the mesh_sample folder. 
   The sphere_m1_image.png is the visulalization result of the sphere_m1.msh.
   The d3b1_bv_solver solves boundary integral equations with the specified datafiles, outputs the results to binary file with the output datafile name.
   The ipw.txt is the sample of incident field datafile, a plane-wave is defined in it. Please refer to "multi_fbeam" for detail. 
   This sample requires about 32GB of memory to solve boundary integral equations (required to run d3b1_bv_solver). 
   The d3b1_bv_solver has optional arguments for rotation and translation of the object. 
   When the vector defining rotation axis is (rx, ry, rz), the rotation angle is theta, the translation vector is (tx, ty, tz), 
   the arguments are './d3b1_bv_solver medium_data.txt sphere_m1.msh ex.dat rx ry rz theta tx ty tz'. 
   Rodrigues' rotation formula is used. 
   As a simple representation of the analysis model, the nodes used for the surface integral are output as point cloud data. 
   In this example, the file ex.particles is output, and the visualization result is ex_particle.png.
   The image was created using gnuplot script "gscript_particles.plt" and converted eps to png by using ImageMagick.

3. type './example1.out' with an argument of datafile name output by d3b1_bv_solver.  
   For example, './example1.out ex.dat'. 
   This executable calculates electromagnetic field, radiation force and torque.  

4. type './example2.out' with an argument of datafile name output by d3b1_bv_solver.   
   For example, './example2.out ex.dat'. 
   This executable calculates electromagnetic field intensity distributions, outputs them to text files. 
   The I_example2.png is the visualization result of intensity distributions, created by using gnuplot script gscript_example2.plt.  
   
5. type './example3.out' with an argument of datafile name output by d3b1_bv_solver.  
   For example, './example3.out ex.dat'. 
   This executable calculates instantaneous value of the electromagnetic fields, outputs them to png image files. 
   The image files are output to the folder which has a name adding "images" to the datafile name specified in the argument (file-extension is excluded). 
   Each image file has a name that indicates the cross section, field component, and number of time steps (ex. xz_Ex_014.png). 
   The color bar is output as color_bar.png in the same folder. 
   The range of color bar in each cross section is output to the info.txt file (ex. xy_info.txt for z=0 plane). 
   The xz_Ex.gif, yz_Ex.gif and xy_Ex.gif are animated gifs that concatenate the png files created by using the shell script gif_animation.sh.  

Please see d3b1_src/bem3_emf_b1.h for detail of functions. 
The main parts of the code are parallelized by using OpenMP. 
The number of threads is controlled by the environment variable OMP_NUM_THREADS. 
The additional analysis examples are in the folder analysis_sample1 ~ analysis_sample4.  

![mesh 0](sphere_m1_image.png "mesh image of the object (sphere_m1_image.png)")  
![point cloud data 0](ex_particles.png "nodes for surface integral (ex_particles.png)") 
![intensity distributions 0](I_example2.png "intensity distributions (I_example2.png)")
![xz_Ex.gif](xz_Ex.gif "instantaneous value of the E_x on y=0 plane (xz_Ex.gif)")![yz_Ex.gif](yz_Ex.gif "instantaneous value of the E_x on x=0 plane (yz_Ex.gif)")  
![xy_Ex.gif](xy_Ex.gif "instantaneous value of the E_x on z=0 plane (xy_Ex.gif)")  


## Analysis sample 1 ( in the folder analysis_sample1 )  

This is the analysis result of plane wave scattering by the cone shaped metal.

![mesh 1](analysis_sample1/cone_m1_image.png "mesh image of the cone (analysis_sample1/cone_m1_image.png)")  
![point cloud data 1](analysis_sample1/ex1_particles.png "nodes for surface integral (analysis_sample1/ex1_particles.png)") 
![intensity distributions 1](analysis_sample1/I_example2.png "intensity distributions (analysis_sample1/I_example2.png)")  
![xz_Ex.gif 1](analysis_sample1/xz_Ex.gif "instantaneous value of the E_x on y=0 plane (analysis_sample1/xz_Ex.gif)")![yz_Ex.gif 1](analysis_sample1/yz_Ex.gif "instantaneous value of the E_x on x=0 plane (analysis_sample1/yz_Ex.gif)")  
![xy_Ex.gif 1](analysis_sample1/xy_Ex.gif "instantaneous value of the E_x on z=0 plane (analysis_sample1/xy_Ex.gif)")  


## Verifications  

### Verification 1    

The verification result using "emf_mie_mmls" is in the folder verification1.
It is the analysis result of plane wave scattering by the two-layered sphere.
The sphere_m2_image.png is the visualization result of mesh datafile. 
The I_example2.png is the visualization result of intensity distributions.
The result of "emf_mie_mmls" is in the folder emf_mie_mmls_result.  

![mesh v1](verification1/sphere_m2_image.png "mesh image of the two-layered sphere (verification1/sphere_m2_image.png)")  

### Verification 2  

The verification result using "emf_mie_ms" is in the folder verification2.
It is the analysis result of plane wave scattering by the three arranged spheres.
The sphere_m3_image.png is the visualization result of mesh datafile. 
The I_example2.png is the visualization result of intensity distributions.
The result of "emf_mie_ms" is in the folder emf_mie_ms_result.  

![mesh v2](verification2/sphere_m3_image.png "mesh image of the three arranged spheres (verification2/sphere_m3_image.png)")  


## About mesh file

This code can use quadrangular ( bi-linear ) and triangular ( linear triangular ) elements. 
I recommend using quadrangular element for reduce required memory. 
The samples of mesh data are in the folder mesh_sample. 
The file with extension .geo is the Gmsh geometry file. 
The file with extension .msh is the mesh datafile created by Gmsh geometry file. 
These mesh files are created by the command 'gmsh -2 -tol 1.0e-15 -format msh2 xxxx.geo' in command line ( xxxx.geo is a geometry file). 
The msh file format must be version 2. 
The domain number ( Physical Surface ) 99 is assigned to the open region in Gmsh geometry file, because Gmsh can't use the number 0 (assigned to open region in the code). 
Please refer to the manual of Gmsh for detail of geometry file.  


## System of units  

This program use the own defined system of units (OSU), optimized for optics. 
The system of units is defined as <img src="https://latex.codecogs.com/gif.latex?c_0=1"> ( speed of light in vacuum ), 
<img src="https://latex.codecogs.com/gif.latex?\mu_0=1"> ( permeability of vacuum ). 
For the conversion from OSU to MKSA system of units, the unit of length in OSU is defined as 
<img src="https://latex.codecogs.com/gif.latex?1\times10^{-6}"> [m] in MKSA, the unit of power in OSU is defined as
<img src="https://latex.codecogs.com/gif.latex?1\times10^{-3}"> [W] in MKSA. The conversions of base unit are follows.  
<img src="https://latex.codecogs.com/gif.latex?a=1\times10^{-6}">,  
<img src="https://latex.codecogs.com/gif.latex?b=1\times10^{-3}">,  
<img src="https://latex.codecogs.com/gif.latex?a\,\mathrm{[m]}=1\,\mathrm{[L]}">,  
<img src="https://latex.codecogs.com/gif.latex?\frac{ab}{c_0^3}\,\mathrm{[kg]}=1\,\mathrm{[M]}">,  
<img src="https://latex.codecogs.com/gif.latex?\frac{a}{c_0}\,\mathrm{[s]}=1\,\mathrm{[T]}">,  
<img src="https://latex.codecogs.com/gif.latex?\sqrt{\frac{b}{c_0\mu_0}}\,\mathrm{[A]}=1\,\mathrm{[I]}">.  
Please see com_src/osu_mksa.h and com_src/osu_mksa.c for detail of conversions.  


## References  

1. Intel Math Kernel Library [MKL](https://software.intel.com/mkl)  
2. The official PNG reference library [libpng](http://www.libpng.org/pub/png/libpng.html)  
3. Three-dimensional mesh generator [Gmsh](https://gmsh.info/)  
4. The command-line driven graphing utility [gnuplot](http://www.gnuplot.info/)  
5. The utilities for manipulating images [ImageMagick](https://imagemagick.org/)  
6. The electromagnetic field analysis program [multi_fbeam](https://github.com/akohta/multi_fbeam/)   
7. The electromagnetic field analysis program [emf_mie_mmls](https://github.com/akohta/emf_mie_mmls/)  
8. The electromagnetic field analysis program [emf_mie_ms](https://github.com/akohta/emf_mie_ms/)  
