# BcLOV4-FEM
Code used for 3D finite element analysis for BcLOV4 spatiotemporal dynamics.

<b>Introduction:</b>

Dynamic membrane recruitment of a cytosol-sequestered protein is one of the most generalizable and ubiquitously applied approaches in optogenetics to control intracellular signaling. It is particularly well-suited to creating inducible peripheral membrane protein systems with the requisite spatiotemporal precision for systematically probing the signaling dynamics that govern how cells respond to external stimuli and regulate their cytoskeletal architecture. The recruitment can be mediated by a protein-protein interaction (PPI) as with optogenetic heterodimerization systems, or by a light-regulated protein-lipid interaction (PLI) between a photoreceptor and the plasma membrane itself, as we have previously shown with BcLOV4, a LOV (light-oxygen-voltage) protein from B. cinerea.

While the photosensory domains that functionally enable the optogenetic tools are often envisioned as binary signaling switches for the fused effector, they themselves are complex signaling proteins whose respective photochemical and biophysical properties give rise to varied spatiotemporal dynamics that have physiological consequences.  For example, if the photocycle or lifetime of the photoactive signaling state is long, the tool beneficially requires less light for sustained signaling, but at the expense of temporal precision and spatial resolution due to persistent activation and lateral diffusion after stimulus removal. A validated kinetic model that accurately describes the non-equilibrium dynamics across an entire three-dimensional cell advances our biophysical understanding of the optogenetic membrane recruitment process and guides the forward-design of spatiotemporally complex signaling outputs. 

Here, we present FEM code that describes the cell-wide 3D spatiotemporal dynamics of the BcLOV4 membrane recruitment system15 in response to pulsatile and spatially patterned illumination. Briefly, a mesh is autogenerated from an image or volumetric image stack of a cell unique in its membrane contour and then seeded with the initial empirical subcellular protein distribution profile. Using experimentally derived constants for every biophysical parameter, the model can then be run to predict the spatiotemporal behavior of the system for any arbitrary light input.  Model performance was analyzed by comparison to corresponding single-cell video data after convolution  with the point-spread-function (PSF) of the imaging system. Importantly, the model mirrors single-cell data as opposed to fitting an idealized cell to population averages, corrects for the point-spread-function (PSF) of the imaging system for more accurate correlation between model output and single-cell data, and represents the global minimum solution since every biophysical parameter has been experimentally determined. 

The model recapitulates salient features observed experimentally, including near diffusion-limited membrane association, protein depletion from distal non-illuminated regions, and dependency of signaling efficiency on optical stimulation intensity and timing. This work greatly informs how optical inputs shape optogenetic signaling outputs and provides a foundation upon which to further construct more complex spatiotemporally accurate models of peripheral membrane protein signaling.

<b>Dependencies:</b>
  
MATLAB with PDE ToolBox. Development was done with MATLAB R2019b, may work for earlier version.

<b>Instructions for Use:</b>
