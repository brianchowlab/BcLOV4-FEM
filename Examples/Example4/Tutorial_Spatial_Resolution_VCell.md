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

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Geom.png?raw=true" align="center" width="450" ></a>
<br/>

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Geom-3D.png?raw=true" align="center" width="450" ></a>
<br/>

## Setting up time- and space-varying excitation
<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Initial.png?raw=true" align="center" width="450" ></a>
<br/>

<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Eq.png?raw=true" align="center" width="450" ></a>
<br/>

## Running simulation
<br/>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example4/Images/Simul.png?raw=true" align="center" width="450" ></a>
<br/>

## Analyzing Results and concerting bac to MATLAB
