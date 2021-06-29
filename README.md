# BcLOV4-FEM
Code used for 3D finite element analysis for BcLOV4 spatiotemporal dynamics.

<b>Introduction:</b>

Dynamic membrane recruitment of a cytosol-sequestered protein is one of the most generalizable and ubiquitously applied approaches in optogenetics to control intracellular signaling. While the photosensory domains that functionally enable the optogenetic tools are often envisioned as binary signaling switches for the fused effector, they themselves are complex signaling proteins whose respective photochemical and biophysical properties give rise to varied spatiotemporal dynamics that have physiological consequences.  Here, we present FEM code that describes the cell-wide 3D spatiotemporal dynamics of the BcLOV4 membrane recruitment system in response to pulsatile and spatially patterned illumination. 

<b>Dependencies:</b>
  
-MATLAB with PDE Toolbox. Development was done with MATLAB R2019b, may work for earlier versions. <br/>
-FELICITY toolbox for MATLAB (https://www.mathworks.com/matlabcentral/fileexchange/31141-felicity). Packaged in this GitHub repository. <br/>
-convolution3D_FFTdomain (https://www.mathworks.com/matlabcentral/fileexchange/35613-3d-convolution-in-the-fft-domain?s_tid=srchtitle). Used for fast 3D convolution of model output with microscope PSF. Packaged in ths GitHub repository. <br/>
-interparc (https://www.mathworks.com/matlabcentral/fileexchange/34874-interparc?s_tid=srchtitle). Used for interpolation along the membrane. Packaged in ths GitHub repository. <br/>
-ndSparse (https://www.mathworks.com/matlabcentral/fileexchange/29832-n-dimensional-sparse-arrays?s_tid=srchtitle). Used for efficient memory handling. Packaged in ths GitHub repository. <br/>
-readImageJROI (https://www.mathworks.com/matlabcentral/fileexchange/32479-readimagejroi?s_tid=srchtitle). Used for reading in ROIs from imageJ. Packaged in this GitHub repository. <br/>
-mesh_xsections (https://www.mathworks.com/matlabcentral/fileexchange/70238-mesh_xsections?s_tid=srchtitle). Used for membrane interpolation. Packaged in this GitHub repository. <br/>


<b>Instructions for Use:</b>

Organization:
The main driver script is called "BcLOVModelwtihExcitation_3D_FELICITY.m," while the associated functions are placed within the "Functions" folder.

Minimal steps to ensure functionality:
Download the repository and open up the driver script. Add the "Functions" folder to your path.

Parameter files:

im_file XYZ.tif                     Image file name. <br/>
il_roi_file XYZ.roi                 Illumination ROI, can be NA.  <br/>
nucleus_roi_file XYZ.roi            Nucleus ROI, can be NA. <br/>
cyto_roi_file XYZ.roi               Cytosol ROI, can be NA. <br/>
concentration_im_file XYZ.tif       Image from which to extract concentration values from. <br/>
mesh XYZ.mat                        For loading in a pre-generated mesh (i.e. from confocal z-stacks). <br/>
z_resolution 0.15                   Z-stack interval in um. Relevant for interpolation. <br/>
D 8.932                             Cytosolic diffusivity in um^2/s. <br/>
S 6000                              Number of membrane binding sites/um^2.<br/> 
D_m 0.028                           Lateral membrane diffusivity in um^2/s. <br/>
k_off_p 0.0541                      Photoreversion rate in s-1. <br/>
k_on_d  1126                        Membrane on rate in dark state in M-1 s-1. <br/>
k_on_l 31901                        Membrane on rate in lit state in M-1 s-1. <br/>
k_off_d 0.0225                      Membrane off rate in dark state in s-1. <br/>
k_off_l 0.038                       Membrane off rate in lit state in s-1. <br/>
power_density 1                     Irradiance at the focal plane in W/cm^2. <br/>
excitation_type 1                   Excitation type: 0 for wholefield, 1 for structured. <br/>
quantum_yield_signaling_state 0.726 Quantum yield of transition from signaling state formation (i.e. absorbs photon and transitions from dark to lit state). <br/>
duty_cycle 1.                       Duty cycle of stimulation in %. <br/>
min_element_size 1                  Minimum tetrahedral element volume in um^3. <br/>
max_element_size 5                  Maximum tetrahedral element volume in um^3. <br/>
downsample 1                        Downsampling for interpolation. Settign to 2 downsamples by a factor of 2. <br/>
scale_len 0.1                       Image pixel size in um. <br/>
hn 4                                Maximum height of nucleus in um in the + and - direction (i.e. top-to-bottom height will be 2x this value). <br/>
h 6                                 Maximum height of cytoplasm in um in the + and - direction (i.e. top-to-bottom height will be 2x this value). <br/>
extrude 30                          Governs number of points along the arc of the hemi-ellipsoid extrusion when extrapolating a volume. <br/>
excitation_wavelength 450           Excitation wavelength in nm. <br/>
extinction_coeff 12500              Extinction coefficient of flavin at excitation wavelength in M-1 cm-1. <br/>
NA 1.4                              Objective lens numerical aperture. <br/>
mag 63                              Objective lens magnification. <br/>
immersion_n 1.52                    Objective lens refractive index of immersion medium. <br/>
dt 1e-1                             Time step for ODE/PDE solver in s. <br/>
num_steps 600                       Number of time steps. <br/>
store_interval 1                    How many time steps between stored solutions. <br/>
interpolation_interval 1            For interpolation steps, interpolate every nth stored solution. <br/>
period 10                           Periodicity of excitation (s). <br/>
theta 0.2929                        Value of theta for solver.<br/> 
conc_ratio 500                      Ratio of cellular fluorescence to concentration (uM). <br/>
offset 0                            Background fluorescence. <br/>
tol 1e-6                            Tolerance for Newton's Method. <br/>
alpha_radius 4                      For meshing. Higher values give more stable meshing, but decrease mesh concavity. <br/>
plot 1                              Plot results?  <br/>
debug 0                             Generates more intermediate output. <br/>


