# sherwood_p3d
Measurements of clustering of Lyman alpha and halos in the Sherwood simulations. 

### Content

This repository contains the measured power spectra (flux auto, halo auto, and flux-halo cross) from the Sherwood suite of simulations, as presented in Givans et al. (2022).

The actual measurements are stored under the data/ folder, in FITS format. (Ignore the pickled_data/ folder, will be removed before publication).

The py/ folder contain some simple scripts to load the data and make some basic plots, including:
 - sherwood_simulation.py: basic objects describing a particular simulation box and 3D grid: box size, resolution, redshift, etc.
 - measured_power.py: book-keeping functions to find a particular power spectrum measurement from the data/ folder, and return a dictionary.
 - model_density.py: simple object to predict the linear power spectrum corresponding to the simulation, at a given redshift.
 - plot_data.py: example plotting script, reproducing Fig 2 of Givans et al. (2022).
 - plot_ratio.py: another plotting script, dividing the measured power by the linear power spectrum.
 - unpickle.py  (Ignore, will be removed soon)

### Install

The only dependency should be fitsio, that one can easily obtain from pip.

You will need to define an environmental variable SHERWOOD pointing to the local copy of the repository.
