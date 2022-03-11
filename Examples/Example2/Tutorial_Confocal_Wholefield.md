# Tutorial 2: Mesh generation from confocal z-stacks and running simulations on confocal images with wholefield excitation

This tutorial will describe how to run a simulation using a z-stack as input. Here we will be using wholefield laser excitation.

## Parameter set
Parameters have been defined in the 'params.txt' file.

## Loading images
Open the "BcLOVModelwithExcitation_3D_FELICITY.m" file and run the first two cells, termed "Set Parameters" and "Load Image". You should see something like this pop up:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/Segmentation.png?raw=true" align="center" width="450" ></a>
<br/>
Note the different outlines regions here. The yellow marks the nucleus, the red marks the cytoplasm, and the blue marks the structure illumination regions.

## Mesh and excitation volume generation
Now run the next cell, which will extrapolate a hemi-ellipsoid geometry from the 2D cell image. Run the next code cell. Something like this should pop-up:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/Mesh.png?raw=true" align="center" width="450" ></a>
<br/>


This shows the resultant geometry and mesh, as well as iso-surfaces of the resultant excitation volume. Check that these are correct before proceeding.

## Running code
The next cell "Build FEM matrices with FELICITY and solve," will compile the FEM matrices and run the simulation. Start it now. After it finishes running you should get a graphical output of the results. It will appear something like this (graph shown is for a different cell):

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/ExcitationROI.png?raw=true" align="center" width="450" ></a>

<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/OutsideROI.png?raw=true" align="center" width="450" ></a>
<br/>

One of these graphs displays the average cytosolic concentration and membrane density over the simulation time across the entire cell. The other displays the average within the excitation region. Note that all the raw simulation data is available within the "Soln" variable for any further custom processing you wish to do.

## Image generation
You may want to generate cross-sections through the 3D result to visualize what it would look like under the microscope. You can do this by running the next cell, which interpolates the mesh results onto a 2D grid at the z-position of interest. In this case, we set the z-location as z = 5 um, i.e.a bit above the focal plabe, to generate the following images (shown t time = 1 s, 51 s, and 200 s).

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/TimePoints/1.jpg?raw=true" align="left" width="250" ></a>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/TimePoints/51.jpg?raw=true" align="left" width="250" ></a>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/TimePoints/201.jpg?raw=true" align="left" width="250" ></a>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
