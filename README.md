# BcLOV4-FEM
Code used for 3D finite element analysis for BcLOV4 spatiotemporal dynamics.

## Introduction:

Dynamic membrane recruitment of a cytosol-sequestered protein is one of the most generalizable and ubiquitously applied approaches in optogenetics to control intracellular signaling. While the photosensory domains that functionally enable the optogenetic tools are often envisioned as binary signaling switches for the fused effector, they themselves are complex signaling proteins whose respective photochemical and biophysical properties give rise to varied spatiotemporal dynamics that have physiological consequences.  Here, we present FEM code that describes the cell-wide 3D spatiotemporal dynamics of the BcLOV4 membrane recruitment system in response to pulsatile and spatially patterned illumination. 

## Capabilities:
This model is designed to be relatively flexible and easily extendable beyond the implementation shown here. Currently, this repository support 3D spatio-temporal modeling of BcLOV4 expressing cells that are excited with wholefield or spatially patterned (e.g. DMD) light. Current supported microscopy setups are widefield and confocal. The main input data is either a 2D micrograph of a cell or a 3D z-stack of a cell. This informaiton is then used to construct a 3D mesh of the cell. In the case of 2D image data this is done by hemi-ellipsoid projection versus in the case of a z-stack this is done by segementationa nd interpolation between segmented frames. See Tutorial 2 for more informationr egarding mesh generation from z-stacks. Stay tuned for our upcoming paper regarding applications of this toolbix.

## Dependencies:
  
-MATLAB with PDE Toolbox. Development was done with MATLAB R2019b, may work for earlier versions. <br/>
-FELICITY toolbox for MATLAB (https://www.mathworks.com/matlabcentral/fileexchange/31141-felicity). Packaged in this GitHub repository. <br/>
-convolution3D_FFTdomain (https://www.mathworks.com/matlabcentral/fileexchange/35613-3d-convolution-in-the-fft-domain?s_tid=srchtitle). Used for fast 3D convolution of model output with microscope PSF. Packaged in ths GitHub repository. <br/>
-interparc (https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc?s_tid=srchtitle). Used for interpolation along the membrane. Packaged in ths GitHub repository. <br/>
-ndSparse (https://www.mathworks.com/matlabcentral/fileexchange/29832-n-dimensional-sparse-arrays?s_tid=srchtitle). Used for efficient memory handling. Packaged in ths GitHub repository. <br/>
-readImageJROI (https://www.mathworks.com/matlabcentral/fileexchange/32479-readimagejroi?s_tid=srchtitle). Used for reading in ROIs from imageJ. Packaged in this GitHub repository. <br/>
-mesh_xsections (https://www.mathworks.com/matlabcentral/fileexchange/70238-mesh_xsections?s_tid=srchtitle). Used for membrane interpolation. Packaged in this GitHub repository. <br/>


## Instructions for Use:

### Organization:
The main driver script is called "BcLOVModelwtihExcitation_3D_FELICITY.m," while the associated functions are placed within the "Functions" folder.

### Minimal steps to ensure functionality:
Download the repository and open up the driver script. Add the "Functions" folder to your path.

### Parameter files:
The parameter file governs the behavior of the simulation and is where you input the cell image, seed the concentration, set the simulation conditions (time-step, mesh size, etc), and define the parameter values for the model. Currently, the following things are editable:

| Entry | Value| Description|
| ------------- | ------------- | ------------- |
|im_file|XYZ.tif|Image file name. <br/>|
|il_roi_file| XYZ.roi |Illumination ROI, can be NA.  <br/>|
|nucleus_roi_file| XYZ.roi        |    Nucleus ROI, can be NA. <br/>|
|cyto_roi_file| XYZ.roi          |     Cytosol ROI, can be NA. <br/>|
|concentration_im_file| XYZ.tif |      Image from which to extract concentration values from. <br/>|
|mesh| XYZ.mat|                        For loading in a pre-generated mesh (i.e. from confocal z-stacks). <br/>|
|z_resolution| 0.15|                   Z-stack interval in um. Relevant for interpolation. <br/>|
|D| 8.932|                             Cytosolic diffusivity in um^2/s. <br/>|
|S| 6000|                              Number of membrane binding sites/um^2.<br/> |
|D_m|0.028          |                 Lateral membrane diffusivity in um^2/s. <br/>|
|k_off_p|0.0541    |                  Photoreversion rate in s-1. <br/>|
|k_on_d|1126     |                   Membrane on rate in dark state in M-1 s-1. <br/>|
|k_on_l|31901    |                    Membrane on rate in lit state in M-1 s-1. <br/>|
|k_off_d|0.0225 |                     Membrane off rate in dark state in s-1. <br/>|
|k_off_l|0.038 |                      Membrane off rate in lit state in s-1. <br/>|
|power_density| 1    |                 Irradiance at the focal plane in W/cm^2. <br/>|
|excitation_type| 1|                   Excitation type: 0 for wholefield, 1 for structured. <br/>|
|quantum_yield_signaling_state |0.726| Quantum yield of transition from signaling state formation (i.e. absorbs photon and transitions from dark to lit state). <br/>|
|duty_cycle |1                |       Duty cycle of stimulation in %. <br/>|
|min_element_size| 1         |         Minimum tetrahedral element volume in um^3. <br/>|
|max_element_size |5        |          Maximum tetrahedral element volume in um^3. <br/>|
|downsample |1                     |   Downsampling for interpolation. Settign to 2 downsamples by a factor of 2. <br/>|
|scale_len_x| 0.1                   |    Image pixel size in um. <br/>|
|scale_len_y| 0.1                   |    Image pixel size in um. <br/>|
|hn |4                           |     Maximum height of nucleus in um in the + and - direction (i.e. top-to-bottom height will be 2x this value). <br/>|
|h |6                          |       Maximum height of cytoplasm in um in the + and - direction (i.e. top-to-bottom height will be 2x this value). <br/>|
|extrude |30                    |      Governs number of points along the arc of the hemi-ellipsoid extrusion when extrapolating a volume. <br/>|
|excitation_wavelength|450   |        Excitation wavelength in nm. <br/>|
|extinction_coeff|12500       |       Extinction coefficient of flavin at excitation wavelength in M-1 cm-1. <br/>|
|NA |1.4 |                             Objective lens numerical aperture. <br/>|
|mag |63  |                            Objective lens magnification. <br/>|
|immersion_n| 1.52 |                   Objective lens refractive index of immersion medium. <br/>|
|dt| 1e-1            |                 Time step for ODE/PDE solver in s. <br/>|
|num_steps |600 |                      Number of time steps. <br/>|
|store_interval| 1            |        How many time steps between stored solutions. <br/>|
|interpolation_interval| 1   |         For interpolation steps, interpolate every nth stored solution. <br/>|
|period |10     |                      Periodicity of excitation (s). <br/>|
|theta|0.2929   |                     Value of theta for solver.<br/> |
|conc_ratio |500|                      Ratio of cellular fluorescence to concentration (uM). <br/>|
|offset |0        |                    Background fluorescence. <br/>|
|tol |1e-6      |                      Tolerance for Newton's Method. <br/>|
|alpha_radius |4 |                     For meshing. Higher values give more stable meshing, but decrease mesh concavity. <br/>|
|plot| 1 |                             Plot results?  <br/>|
|debug |0  |                           Generates more intermediate output. <br/>|


## Example:
For this example, we will use an idealized circular cell with concentric nuclear and plasma membrane borders. The nuclear radius is 15 um and cytoplasmic radius is 25 um. This cell can be found in ./Example/OpticsModelCell.

### Step 1: Loading parameters
Open up the main driver script and set the first line as:<br/>

filename = './Example/params';<br/>

Then run the first code block. Do not forget to add the "Functions" folder to your path. This will load the parameter file which can be found at "./Exmaples/params.txt".<br/>

### Step 2: Loading cell image and selecting ROIs
Run: <br/>

[contours,mask_il,I] = LoadImages(param); <br/>

This will load up the cell image. Draw a box in the center of the nucleus to seed the activecontour model so that it can find the nucleus. Let the code run for a few seconds and you should get a result like this, where the cytoplasm, nucleus, and excitation ROI are all highlighted: <br/>

<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example3/Outputs/Cell_regions.png?raw=true" align="left" width="450" ></a>

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
<br/>

### Step 3: Meshing
Run: <br/>

[mesh_c,poly,shp_n] = GenMesh(contours,param);<br/>

This will generate the mesh using the prameters set in the parameter file and using a geometry generated by hemi-ellipsoid projection of the 2D cell image given above. The result should appear like this: <br/>


<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example3/Outputs/Mesh.png?raw=true" align="left" width="350" ></a>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example3/Outputs/Mesh_alt.png?raw=true" align="left" width="350" ></a>

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


### Step 4: Illumination ROI
Run: <br/>

[photo_on_scale,idx_excited] = ExcitationROI(mesh_c,mask_il,poly,param);

This will predict the excitation volume from the illumination ROI. The result should look like this:  

<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example3/Outputs/Il.png?raw=true" align="left" width="450" ></a>

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
<br/>
<br/>
<br/>

### Step 5: Run simulation
Running the next code block will start the simulation.


### Step 6: Output
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/?raw=true" align="left" width="450" ></a>
<a href="url"><img src="https://github.com/brianchowlab/BcLOV4-FEM/blob/main/Examples/Example2/Images/?raw=true" align="left" width="450" ></a>        

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
<br/><br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>
<br/>


## Tutorials:

### Tutorial 1: Widefield input with structured illumination
<a href='https://github.com/brianchowlab/BcLOV4-FEM/edit/main/Examples/Example1/Tutorial_Widefield_DMD.md'>This tutorial</a> demonstrates how to set up an experiment using an actual cell. Here, the images are assumed to be acquired under widefield, so z-stack data is not available and the resultant mesh has to be generated by hemi-ellipsoid projection. Here the excitation is structured.

### Tutorial 2: Confocal z-stack input with wholefield illumination
<a href='https://github.com/brianchowlab/BcLOV4-FEM/edit/main/Examples/Example2/Tutorial_Confocal_Wholefield.md'>This tutorial</a> demonstrates how to set up an experiment where the cell is imaged visa confocal and z-stack data is available. The resultant mesh is constructured by segmentaiton of z-stack images and interpolation in the z-plane. Here, the excitation is assumed to be wholefield, rather than structured as above.
