# Spectral Indices
Spectral indices to measure the depth of the main features in the mid-infrared spectra of ultracool objects (see figure below):
* silicate_index: silicate index for the silicate absorption feature at 9.3 microns as defined in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract) and improved in [Suárez & Metchev (2023)](https://ui.adsabs.harvard.edu/abs/2023MNRAS.523.4739S/abstract).
* water_index: water index for the water absorption feature at 6.25 microns as defined in [Cushing et al. (2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...648..614C/abstract) and modified in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract).
* methane_index: methane index for the methane absorption feature at 7.65 microns as defined in [Cushing et al. (2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...648..614C/abstract) and modified in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract).
* ammonia_index: ammonia index for the ammonia absorption feature at 10.5 microns as defined in [Cushing et al. (2006](https://ui.adsabs.harvard.edu/abs/2006ApJ...648..614C/abstract) and modified in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract).

Collection of Spitzer IRS mid-infrared spectra in [Suárez & Metchev (2022)](https://ui.adsabs.harvard.edu/abs/2022MNRAS.513.5701S/abstract) highlighting the main features at different spectral types. This repository allows you to measure the strength of the indicated water, methane, amonnia, and silicates features.

[//]: # (This is a comment.)
![alt text](https://github.com/suarezgenaro/spectral_indices/blob/main/Spitzer_IRS_spectra.png)

## Silicate Index
	Description:
	------------
		Measure the strength of the mid-infrared silicate absorption considering the defined silicate index in Suárez & Metchev (2022,2023).

	Parameters:
	-----------
	- wl : array
		Spectrum wavelengths in microns.
	- flux : array
		Spectrum fluxes in Jy.
	- eflux : array
		Spectrum flux errors in Jy.
	- silicate_peak : float, optional
		Wavelength reference to indicate the center of silicate absorption. 
		Default value is 9.3 microns (Suárez & Metchev 2023).
	- silicate_window : float, optional
		Wavelength window around ``silicate_peak`` used to calculate the median flux at the absorption.
		Default value is 0.6 microns (Suárez & Metchev 2023).
	- silicate_cont1 : float, optional
		Wavelength reference to indicate the short-wavelength continuum of silicate absorption. 
		Default value is 7.45 microns (Suárez & Metchev 2023).
	- silicate_window_cont1 : float, optional
		Wavelength window around ``silicate_cont1`` used to calculate the median flux at the short-wavelength absorption continuum.
		Default value is 0.5 microns (Suárez & Metchev 2023).
	- silicate_cont2 : float, optional
		Wavelength reference to indicate the long-wavelength continuum of silicate absorption. 
		Default value is 13.5 microns (Suárez & Metchev 2023).
	- silicate_window_cont2 : float, optional
		Wavelength window around ``silicate_cont2`` used to calculate the median flux at the long-wavelength absorption continuum.
		Default value is 1.0 microns (Suárez & Metchev 2023).
	- default : {``SM22``, ``SM23``}, optional (default ``SM23``)
		Reference to set default parameters to measure the silicate index.
		``SM22`` for Suárez & Metchev (2022) or ``SM23`` (default) for Suárez & Metchev (2023).
	- plot : {``True``, ``False``}, optional (default ``False``)
		Plot (``True``) or do not plot (``False``) the silicate index measurement.
	- plot_title : str, optional
		Plot title (default ``'Silicate Index Measurement'``.
	- plot_xrange : list or array
		Wavelength range (in microns) of the plot (default [5.2, 14] um).
	- plot_yrange : list or array
		Flux range (in Jy) of the plot (default [5.2, 14] um).
	- plot_save : {``True``, ``False``}, optional (default ``True``)
		Save (``'True'``) or do not save (``'False'``) the resulting plot.

	Returns:
	--------
	- Dictionary 
		Dictionary with silicate index parameters:
			- ``'silicate'`` : Silicate index
			- ``'esilicate'`` : Silicate index uncertainty
			- ``'flux_silicate_peak'`` : flux at the absorption feature
			- ``'eflux_silicate_peak'`` : flux error at the absorption feature
			- ``'cont_silicate'`` : flux at the continuum of the absorption
			- ``'econt_silicate'`` : flux uncertainty at the continuum of the absorption
			- ``'slope'`` : Slope of the linear fit to the log(flux)-log(wavelength) space
			- ``'eslope'`` : Slope uncertainty
			- ``'constant'`` : Constant or intercept of the linear fit
			- ``'econstant'`` : Constant uncertainty
			- ``'silicate_peak'`` : input ``silicate_peak``
			- ``'silicate_window'`` : input ``silicate_window``
			- ``'silicate_con1'`` : input ``silicate_con1``
			- ``'silicate_window_con1'`` : input ``silicate_window_con1``
			- ``'silicate_con2'`` : input ``silicate_con2``
			- ``'silicate_window_con2'`` : input ``silicate_window_con2``
			- ``'wl'`` : input ``wl``
			- ``'flux'`` : input ``flux``
			- ``'eflux'`` : input ``eflux``
	- Plot of the silicate index measurement that will be stored if ``save`` with the name ``out_file``.
