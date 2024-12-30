# Spectral Indices
Spectral indices to measure the depth of the main features in the mid-infrared spectra of ultracool objects (see figure below).

Collection of Spitzer IRS mid-infrared spectra in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract) highlighting the main features at different spectral types. This repository allows you to measure the strength of the indicated water, methane, amonnia, and silicates features.

[//]: # (This is a comment.)

<p align="center">
    <img src="aux/Spitzer_IRS_spectra.png" title="Spitzer IRS spectra of ultracool objects" alt="Spitzer IRS spectra of ultracool objects" width="500">
</p>

Function included to measure spectral indices:
* silicate_index: silicate index for the silicate absorption feature at 9.3 microns as defined in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract) and improved in [Suárez & Metchev (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523.4739S/abstract).
* water_index: water index for the water absorption feature at 6.25 microns as defined in [Cushing et al. (2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...648..614C/abstract) and modified in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract).
* methane_index: methane index for the methane absorption feature at 7.65 microns as defined in [Cushing et al. (2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...648..614C/abstract) and modified in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract).
* ammonia_index: ammonia index for the ammonia absorption feature at 10.5 microns as defined in [Cushing et al. (2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...648..614C/abstract) and modified in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract).

## Tutorial
[This tutorial](https://github.com/suarezgenaro/spectral_indices/blob/main/aux/example_index_measurements.ipynb) explains how to obtain silicate, water, methane, and ammonia indices for mid-infrared spectra of brown dwarfs, as shown below:

<p align="center">
    <img src="aux/silicate_index_measurement.png" title="Silicate index" width="400">
    <img src="aux/water_index_measurement.png" title="Water index" width="400">
    <img src="aux/methane_index_measurement.png" title="Methane index" width="400">
    <img src="aux/ammonia_index_measurement.png" title="Ammonia index" width="400">
</p>

## Silicate Index

- [x] #739
- [ ] https://github.com/suarezgenaro/seda/issues/19
- [ ] Add delight to the experience when all tasks are complete :tada:

- [x] #739
- [ ] https://github.com/octo-org/octo-repo/issues/740
- [ ] Add delight to the experience when all tasks are complete :tada:
