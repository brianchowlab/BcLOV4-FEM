# Tutorial 3: Spatial resolution simulation using different excitation volumes with Virtual Cell integration

This tutorial demonstrates how to use the Virtual Cell integration capabilities of this toolbox to simulate the achievable BcLOV4 spatial resolution using different stimulation approuches, specifically 1P (laser-scanning confocal), 2P, and TIRF stimulation

For the purpsoses of this tutorial, we have including Virtual Cell simulation files (.vcml) for 1P, 2P, and TIRF files. Here, we will focus on 2P stimulation (2P_Final.vcml), but the protocol is generalizable between different simulations.

## Setting up Virtual Cell Simulation

The process to setting up a Virtual Cell simulation includes:  
1) Setting up a reaction diagram. 
2) Defining the reactions. 
3) Loading in a MATLAB geometry by conversion to an .STL file. 
4) Setting up initial conditions. 
5) Defining time- and space-varying excitation light. 
6) Running the simulation. 
7) Porting data back to MATLAB for analysis. 

Let's focus on the first 2 steps in this section.

### Reaction Diagram
When you open "2P_Final.vcml" you will see the reaction diagram already setup. This reflects the following reactions:  
1) Interconversion between lit-state membrane protein and dark-state membrane protein. 
2) Interconversion between lit-state cytosolic protein and dark-state cytosolic protein. 
3) Interconversion between lit-state membrane protein and lit-state cytosolic protein. 
4) Interconversion between dark-state membrane protein and dark-state cytosolic protein. 
5) An additional state ("light") is included to allow for a time- and space-varying excitation intput that drives reactions (1) and (2), i.e. the dark->lit conversion.

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/FrontPage.png?raw=true" align="center" width="450" ></a>
<br/>

### Definining reactions
We next define all of our reactions according to the biophysical parameter values for the protein. Navigate to the following page:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Reactions.png?raw=true" align="center" width="450" ></a>
<br/>
<br/>

We define the parameters by clicking on each reaction. The parameters here are equivalent to our inputs into the "params" file in our MATLAB toolbox. For example, a reaction would appear like this:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Reaction-ex.png?raw=true" align="center" width="450" ></a>
<br/>

## Porting mesh from MATLAB to VCell
Virtual Cell has built-in mesh generating capabilities, but this toolbox gives more flexibility in mesh generation, so you may want to port over the toolbox mesh to Virtual Cell. We have created a MATLAB script to allow for this. To do this, open "Mesh2STL.m," load your toolbox mesh into MATLAB, and then run the script. This will create a STL file called 'optics_geom.stl" which defines the mesh in a Virtual Cell compatible format. Then go to the geometry section in Virtual Cell and go to the option to "Replace Geometry." Import the STL file into Virtual Cell. You should not see see something like this:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Geom-3D.png?raw=true" align="center" width="450" ></a>
<br/>

You now need to define which part of the mesh cooresponds to th ecytoplasm, nucleus, plasma membrane, etc. This is also done in the geometry tab and should look like this:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Geom.png?raw=true" align="center" width="450" ></a>
<br/>



## Setting up time- and space-varying excitation
If we have time- and space- varying illumination (i.e. pulsed structures illumination, such as a laser-scanning confocal excitation ROI), we will need to input that into Virtual Cell. Virtual Cell does not natively support such inputs, but we can make it work by defining a time- and space- varying initial condition in the initial conditions tab:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Initial.png?raw=true" align="center" width="450" ></a>
<br/>

We have incldued a MATLAB Script to define the light function. Load the "boolean_VCell.m" script, which gives options to generate rastered 1P, 2P, or TIRF elimination. You may need to edit this scripts for your specific needs. Runt he cooresponding cell, which will generate a TXT file that gives the inout function. Copy the input function into Virtual Cell. It will look somethign like this:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Eq.png?raw=true" align="center" width="450" ></a>
<br/>

Do not forget to define the rest of your initial conditions, e.g. BcLOV4 initial concentration.

## Running simulation
You are now ready to run the simulation in the simulations tab. You can edit the runtime. It will look something like this:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Simul.png?raw=true" align="center" width="450" ></a>
<br/>

Start your simulation and let it run.

## Analyzing Results and converting back to MATLAB
Once the simulation has finished running, you may want to analyze it. This can be done back in MATLAB. First export it to a VTK file. Then use the "read_vtk.m" script to convert it back to a ".mat" file. The result should be something like the "vcell_2P_result.mat."

You can plot it using the "VCell_Plots.m" script. The result should look something like this:

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Plot.png?raw=true" align="center" width="450" ></a>
<br/>

