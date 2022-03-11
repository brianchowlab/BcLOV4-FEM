# Tutorial: Running simulations on widefield images with structured illumination

This tutorial will describe how to run a simulation for a cell imaged using widefield microscopy and stimulated using structured illumination.

## Parameter set
Parameters have been defined in the 'params.txt' file:

| Entry | Value| Description|
| ------------- | ------------- | ------------- |
|im_file|./Examples/Example1/Slice.tif|Image file name. <br/>|
|il_roi_file| /Examples/Example1/ExcitationROI.roi |Illumination ROI, can be NA.  <br/>|
|nucleus_roi_file| X./Examples/Example1/Nucleus.roi       |    Nucleus ROI, can be NA. <br/>|
|cyto_roi_file| ./Examples/Example1/Cytosol.roi          |     Cytosol ROI, can be NA. <br/>|
|concentration_im_file| ./Examples/Example1/Slice.tif |      Image from which to extract concentration values from. <br/>|
|mesh| NA|                        For loading in a pre-generated mesh (i.e. from confocal z-stacks). <br/>|
|z_resolution| 0.15|                   Z-stack interval in um. Relevant for interpolation. <br/>|
|D| 8.932|                             Cytosolic diffusivity in um^2/s. <br/>|
|S| 6000|                              Number of membrane binding sites/um^2.<br/> |
|D_m|0.028          |                 Lateral membrane diffusivity in um^2/s. <br/>|
|k_off_p|0.0541    |                  Photoreversion rate in s-1. <br/>|
|k_on_d|1126     |                   Membrane on rate in dark state in M-1 s-1. <br/>|
|k_on_l|31901    |                    Membrane on rate in lit state in M-1 s-1. <br/>|
|k_off_d|0.0225 |                     Membrane off rate in dark state in s-1. <br/>|
|k_off_l|0.038 |                      Membrane off rate in lit state in s-1. <br/>|
|power_density| 0.012    |                 Irradiance at the focal plane in W/cm^2. <br/>|
|excitation_type| 0|                   Excitation type: 0 for wholefield, 1 for structured. <br/>|
|quantum_yield_signaling_state |0.726| Quantum yield of transition from signaling state formation (i.e. absorbs photon and transitions from dark to lit state). <br/>|
|duty_cycle |1                |       Duty cycle of stimulation in %. <br/>|
|min_element_size| 1         |         Minimum tetrahedral element volume in um^3. <br/>|
|max_element_size |5        |          Maximum tetrahedral element volume in um^3. <br/>|
|downsample |2                    |   Downsampling for interpolation. Settign to 2 downsamples by a factor of 2. <br/>|
|scale_len_x| 0.1                   |    Image pixel size in um. <br/>|
|scale_len_y| 0.1                   |    Image pixel size in um. <br/>|
|hn |3                           |     Maximum height of nucleus in um in the + and - direction (i.e. top-to-bottom height will be 2x this value). <br/>|
|h |65                         |       Maximum height of cytoplasm in um in the + and - direction (i.e. top-to-bottom height will be 2x this value). <br/>|
|extrude |30                    |      Governs number of points along the arc of the hemi-ellipsoid extrusion when extrapolating a volume. <br/>|
|excitation_wavelength|405   |        Excitation wavelength in nm. <br/>|
|extinction_coeff|7000       |       Extinction coefficient of flavin at excitation wavelength in M-1 cm-1. <br/>|
|NA |1.4 |                             Objective lens numerical aperture. <br/>|
|mag |63  |                            Objective lens magnification. <br/>|
|immersion_n| 1.515 |                   Objective lens refractive index of immersion medium. <br/>|
|dt| 1e-1            |                 Time step for ODE/PDE solver in s. <br/>|
|num_steps |2000 |                      Number of time steps. <br/>|
|store_interval| 10            |        How many time steps between stored solutions. <br/>|
|interpolation_interval| 1   |         For interpolation steps, interpolate every nth stored solution. <br/>|
|period |10     |                      Periodicity of excitation (s). <br/>|
|theta|0.2929   |                     Value of theta for solver.<br/> |
|conc_ratio |500|                      Ratio of cellular fluorescence to concentration (uM). <br/>|
|offset |0        |                    Background fluorescence. <br/>|
|tol |1e-6      |                      Tolerance for Newton's Method. <br/>|
|alpha_radius |2 |                     For meshing. Higher values give more stable meshing, but decrease mesh concavity. <br/>|
|plot| 1 |                             Plot results?  <br/>|
|debug |0  |                           Generates more intermediate output. <br/>|


## Loading images
Open the "BcLOVModelwithExcitation_3D_FELICITY.m" file and run the first two cells, termed "Set Parameters" and "Load Image". You should see something like this pop up:

![Segmentation] (Images/Segmentation.png)

Note the different outlines regions here. The yellow marks the nucleus, the red marks the cytoplasm, and the blue marks the structure illumination regions.

## Mesh and excitation volume generation
Now run the next cell, which will extrapolate a hemi-ellipsoid geometry from the 2D cell image. Run the next code cell. Something like this should pop-up:

![Segmentation] (Images/Mesh.png)

This shows the resultant geoemtry and mesh, as well as iso-surfaces of the resultant excitation volume. Check that these are correct before proceeding.

## Running code
The next cell "Build FEM matrices with FELICITY and solve," will compile the FEM matrices and run the simulation. Start it now. After it finishes running you should get a graphical output of the results. It will appear something like this (graph shown is for a different cell):
![Segmentation] (../Example2/Images/ExcitationROI.png)
![Segmentation] (../Example2/Images/OutsideROI.png)

One of these graphs displays the average cytosolic concentration and membrane density over the simulation time across the entire cell. The other displays the average within the excitation region. Note that all the raw simulation data is available within the "Soln" variable for any further custom processing you wish to do.

## Image generation
You may want to generate cross-sections through the 3D result to visualize what it would look like under the microscope.


