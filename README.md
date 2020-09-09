# shallow_water
2D Shallow Water Equation Solver
we git fork form https://github.com/b2kay/shallow_water2D, and we have done some
efforts on refactoring and optimizing 2D Shallow Water Equation Solver on the 
Sunway TaihuLight supercomputer.
## Compilation and Run
 apt-get install libnetcdf libnetcdf-dev
 apt install libnetcdf-c++4-dev libnetcdf-cxx-legacy-dev
 $make
 
 $./main.x
 
or

 $sh Allrun.sh

## plot and save png
 ncview output.nc
 
or 
edit x and y spacing in animation.m. eta*.png files will be saved in png folder)

 $octave animation.m
 
 $convert -delay 10 ./png/*.png anim.gif

## Parameters file
1st row : x spacing

2nd row : y spacing

3rd : timestep

4th : total time

