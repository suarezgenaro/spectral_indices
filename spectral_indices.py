import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from sys import exit

##########################
def silicate_index(wl, flux, eflux, silicate_peak=None, silicate_window=None, silicate_con1=None, 
	               silicate_window_con1=None, silicate_con2=None, silicate_window_con2=None, default='SM23',
	               plot=False, plot_title=None, plot_xrange=None, plot_yrange=None, plot_save=True):
	'''
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
	- Plot of the silicate index measurement that will be stored if ``plot_save``.

	Author: Genaro Suárez
	'''

	# avoid nan values
	mask_nan = ~np.isnan(eflux)
	wl = wl[mask_nan]
	flux = flux[mask_nan]
	eflux = eflux[mask_nan]

	# use default values if optional values (peak and continuum regions) are not provided
	if default=='SM23':
		# values from Suárez & Metchev (2023)
		# center of the absorption and window region
		if (silicate_peak is None): silicate_peak = 9.30 # um
		if (silicate_window is None): silicate_window = 0.6 # um
		# continuum region definition
		# where CH4 and NH3 are not expected for mid-L type
		if (silicate_con1 is None): silicate_con1 = 7.45 # minimum in the blue side of the peak (um)
		if (silicate_window_con1 is None): silicate_window_con1 = 0.5 # window in um for the short-wavelength continuum region
		if (silicate_con2 is None): silicate_con2 = 13.5 # minimum in the red side of the peak (um)
		if (silicate_window_con2 is None): silicate_window_con2 = 1.0 # window in um for the long-wavelength continuum region
	elif default=='SM22':
		# values from Suárez & Metchev (2022)
		# center of the absorption and window region
		silicate_peak = 9.00 # um
		silicate_window = 0.6 # um
		# continuum region definition
		silicate_con1 = 7.5 # minimum in the blue side of the peak (um)
		silicate_window_con1 = 0.6 # window in um for the short-wavelength continuum region
		silicate_con2 = 11.5 # minimum in the red side of the peak (um)
		silicate_window_con2 = 0.6 # window in um for the long-wavelength continuum region
	else: raise Exception(f'"{default}" is not a recognized reference to set default parameters to measure the silicate index. \nTry "SM22" or "SM23".')

	# measure silicate index
	ind_silicate_peak = (wl>=(silicate_peak-silicate_window/2)) & (wl<=(silicate_peak+silicate_window/2))
	# using the mean
	flux_silicate_peak = np.mean(flux[ind_silicate_peak]) # average flux at the bottom of the feature
	eflux_silicate_peak = np.std(flux[ind_silicate_peak]) / np.sqrt(flux[ind_silicate_peak].size) # flux error as the standard error of the mean
	## using the median
	#flux_silicate_peak = np.median(flux[ind_silicate_peak]) # average flux at the bottom of the feature
	#eflux_silicate_peak = 1.2533 * np.std(flux[ind_silicate_peak]) / np.sqrt(flux[ind_silicate_peak].size) # standard error of the median = 1.2533 times the standard error of the mean
	##eflux_silicate_peak = (np.percentile(flux[ind_silicate_peak], 84, interpolation = 'linear') - np.percentile(flux[ind_silicate_peak], 16, interpolation = 'linear')) / np.sqrt(flux[ind_silicate_peak].size) # 68% confidence interval
	##eflux_silicate_peak = flux_silicate_peak * np.median(eflux[ind_silicate_peak]/flux[ind_silicate_peak]) # error to keep the fractional uncertainties

#	# data points in the first continuum region
#	ind_silicate_con1 = (wl>=(silicate_con1-silicate_window_con1/2)) & (wl<=(silicate_con1+silicate_window_con1/2))
#	# data points in the second continuum region
#	ind_silicate_con2 = (wl>=(silicate_con2-silicate_window_con2/2)) & (wl<=(silicate_con2+silicate_window_con2/2))

	# fit a line to the continuum points
	log_log_space = True # True or False to consider the log(flux)-log(lambda) space
	#--------------------------------------
	# fit a line to the continuum points in linear space
#	if (log_log_space=='no'):
#		ind_silicate_con1_con2 = np.concatenate((ind_silicate_con1, ind_silicate_con2)) # indices of the continuum data points on both regions
#		fit, cov_fit = np.polyfit(wl[ind_silicate_con1_con2,i], flux[ind_silicate_con1_con2,i], 1, w=1./eflux[ind_silicate_con1_con2,i], cov=True) # weigh by the inverse of the error
#		slope = fit[0] # slope of the histogram without any correction
#		eslope = np.sqrt(cov_fit[0,0]) # error of the slope
#		constant = fit[1]
#		econstant = np.sqrt(cov_fit[1,1]) 
#		# continuum at the point in of the absorption
#		cont_silicate = slope*silicate_peak + constant
#		econt_silicate = (((slope+eslope)*silicate_peak + constant) - ((slope-eslope)*silicate_peak + constant)) / 2

	
	# fit a line to the continuum points in the log(flux)-log(lambda) space
	if log_log_space:
		ind_silicate_con1_con2 = ((wl>=(silicate_con1-silicate_window_con1/2)) & (wl<=(silicate_con1+silicate_window_con1/2))) | \
		                         ((wl>=(silicate_con2-silicate_window_con2/2)) & (wl<=(silicate_con2+silicate_window_con2/2))) # data points in both continuum regions
		if (flux[ind_silicate_con1_con2].min()<=0): 
			print('Warning: Logarithm of negative and/or zero fluxes in the continuum windows')
			print('   Only positive fluxes will be used')
			ind_silicate_con1_con2 = (((wl>=(silicate_con1-silicate_window_con1/2)) & (wl<=(silicate_con1+silicate_window_con1/2))) | \
			                         ((wl>=(silicate_con2-silicate_window_con2/2)) & (wl<=(silicate_con2+silicate_window_con2/2)))) & \
			                         (flux>0) # positive fluxes in both continuum regions
		logwl = np.log10(wl[ind_silicate_con1_con2])
		logflux = np.log10(flux[ind_silicate_con1_con2])
		#elogflux = (1/np.log(10)) * eflux[ind_silicate_con1_con2]/flux[ind_silicate_con1_con2]
		elogflux = eflux[ind_silicate_con1_con2] # these errors are used only for weighing the fit. Thus, assigning the log errors equal to the linear error will weight by the error bars directly plotted in the flux vs wavelength plot. Otherwise, when properly propagating the linear errors to log errors will give more weight (lower errors) to higher fluxes.
	
		fit, cov_fit = np.polyfit(logwl, logflux, 1, w=1./elogflux, cov=True) # weigh by the inverse of the error
		slope = fit[0] # slope 
		eslope = np.sqrt(cov_fit[0,0]) # slope error 
		constant = fit[1] # constant
		econstant = np.sqrt(cov_fit[1,1]) # constant error
		# the slope error has to be scaled somehow to be used in linear plots
		fit_linear, cov_fit_linear = np.polyfit(wl[ind_silicate_con1_con2], flux[ind_silicate_con1_con2], 1, w=1./eflux[ind_silicate_con1_con2], cov=True) # weigh by the inverse of the error
		slope_linear = fit_linear[0] # slope 
		eslope_linear = np.sqrt(cov_fit_linear[0,0])
		# slope error normalization
		eslope = eslope / np.log10(eslope/eslope_linear) # this is an empirical correction I deduced
	
		# continuum at the absorption peak
		cont_silicate = 10**(slope*np.log10(silicate_peak) + constant)
		econt_silicate = (10**(((slope+eslope)*np.log10(silicate_peak)+constant)) - 10**((slope-eslope)*np.log10(silicate_peak)+constant)) / 2
		#econt_silicate = (10**(((slope+eslope)*np.log10(silicate_peak)+(constant-econstant))) - 10**((slope-eslope)*np.log10(silicate_peak)+(constant+econstant))) / 2
	
	## fit an exponential curve to the continuum points in the flux vs. lambda plot (this should be equivalent to a linear fit in the log(flux)-log(lambda) space)
	# ATTEMPT
	#if not log_log_space:
	#	#def exponential(x, a, b, c):
	# 	#   return a * np.exp(b * x) + c
	#	def exponential(x, a, b):
	# 	   return a*np.exp(b * x)
	#	ind_silicate_con1_con2 = np.concatenate((ind_silicate_con1, ind_silicate_con2)) # indices of the continuum data points on both regions
	#	fit_exp, cov_fit_exp = curve_fit(exponential, wl[ind_silicate_con1_con2,i], flux[ind_silicate_con1_con2,i])#, sigma=eflux[ind_silicate_con1_con2,i])

	# silicate index
	silicate = cont_silicate/flux_silicate_peak
	esilicate = silicate * np.sqrt((econt_silicate/cont_silicate)**2 + (eflux_silicate_peak/flux_silicate_peak)**2)

	# output dictionary
	out = {'silicate': silicate, 'esilicate': esilicate, 'flux_silicate_peak': flux_silicate_peak, 
	       'eflux_silicate_peak': eflux_silicate_peak, 'cont_silicate': cont_silicate, 'econt_silicate': econt_silicate}
	# add parameter of the slope fit
	out['slope'] = slope
	out['eslope'] = eslope
	out['constant'] = constant
	out['econstant'] = econstant
	# add parameters used to measure the index
	out['silicate_peak'] = silicate_peak
	out['silicate_window'] = silicate_window
	out['silicate_con1'] = silicate_con1
	out['silicate_window_con1'] = silicate_window_con1
	out['silicate_con2'] = silicate_con2
	out['silicate_window_con2'] = silicate_window_con2
	#out['log_log_space'] = log_log_space
	out['wl'] = wl
	out['flux'] = flux
	out['eflux'] = eflux

	# visualize how the silicate index was measured
	if plot: plot_silicate_index(out, plot_xrange, plot_yrange, plot_title, plot_save)

	return out

#++++++++++++++++++++++++++++++++++++++
# plot the silicate index measurement
def plot_silicate_index(out_silicate_index, plot_xrange=None, plot_yrange=None, plot_title=None, plot_save=True):

	# read parameters of interest
	wl = out_silicate_index['wl']
	flux = out_silicate_index['flux']
	eflux = out_silicate_index['eflux']
	si_index = out_silicate_index['silicate']
	esi_index = out_silicate_index['esilicate']
	cont_silicate = out_silicate_index['cont_silicate']
	econt_silicate = out_silicate_index['econt_silicate']
	# continuum fit
	slope = out_silicate_index['slope']
	eslope = out_silicate_index['eslope']
	constant = out_silicate_index['constant']
	econstant = out_silicate_index['econstant']
	# silicate feature peak
	silicate_min = out_silicate_index['silicate_peak'] # flux peak (um)
	silicate_window = out_silicate_index['silicate_window'] # window in um for the absorption
	# first window region for the continuum
	silicate_con1 = out_silicate_index['silicate_con1'] # center (in um) of the short-wavelength continuum region
	silicate_window_con1 = out_silicate_index['silicate_window_con1'] # window in um for the short-wavelength continuum region
	# second window region for the continuum
	silicate_con2 = out_silicate_index['silicate_con2'] # center (in um) of the long-wavelength continuum region
	silicate_window_con2 = out_silicate_index['silicate_window_con2'] # window in um for the long-wavelength continuum region

	#++++++++++++++++++++++++++++++
	# plot silicate index measurement
	fig, ax = plt.subplots()
	
	#xmin, xmax = wl[wl>0].min(), wl[wl>0].max()
	#ymin, ymax = np.percentile(flux[flux!=0], scale_percent), np.nanpercentile(flux[flux!=0], 100-scale_percent)
	#xmin, xmax = 0.90*silicate_con1-silicate_window_con1/2., 1.05*silicate_con2+silicate_window_con2/2.
	if plot_xrange is None: xmin, xmax = 5.2, 14
	else: xmin, xmax = plot_xrange
	if plot_yrange is None:
		#scale_percent = 0.1
		mask = (wl>=xmin) & (wl<=xmax)
		ymin, ymax = flux[mask].min(), flux[mask].max()
	else: ymin, ymax = plot_yrange
	
	# indicate the regions to measure the silicate index
	ax.fill([silicate_min-silicate_window/2, silicate_min-silicate_window/2, 
	        silicate_min+silicate_window/2, silicate_min+silicate_window/2], 
	        [ymin, ymax, ymax, ymin], facecolor='silver', linewidth=1, zorder=2)
	ax.fill([silicate_con1-silicate_window_con1/2, silicate_con1-silicate_window_con1/2, 
	        silicate_con1+silicate_window_con1/2, silicate_con1+silicate_window_con1/2], 
	        [ymin, ymax, ymax, ymin], facecolor='gainsboro', linewidth=1, zorder=2)
	ax.fill([silicate_con2-silicate_window_con2/2, silicate_con2-silicate_window_con2/2, 
	        silicate_con2+silicate_window_con2/2, silicate_con2+silicate_window_con2/2], 
	        [ymin, ymax, ymax, ymin], facecolor='gainsboro', linewidth=1, zorder=2)
	
	# plot spectrum
	mask = wl>0
	plt.plot(wl[mask], flux[mask])
	
	# linear fit
	xc_silicate = np.linspace(silicate_con1-silicate_window_con1/2., silicate_con2+silicate_window_con2/2., 100)
	yc_silicate = 10**(slope*np.log10(xc_silicate) + constant)
	ycu_silicate = 10**((slope+eslope)*np.log10(xc_silicate) + constant)
	ycd_silicate = 10**((slope-eslope)*np.log10(xc_silicate) + constant)
	ax.fill(np.append(xc_silicate, xc_silicate[::-1]), np.append(ycd_silicate, ycu_silicate[::-1]), 
	        facecolor='gray', edgecolor='gray', linewidth=0, alpha=0.20)
	ax.plot(xc_silicate, yc_silicate, '--', color='gray', linewidth=0.5)
	
	# fluxes in the continuum regions
	ind_silicate_con1 = np.where((wl>=(silicate_con1-silicate_window_con1/2)) & (wl<=(silicate_con1+silicate_window_con1/2)))[0]
	ind_silicate_con2 = np.where((wl>=(silicate_con2-silicate_window_con2/2)) & (wl<=(silicate_con2+silicate_window_con2/2)))[0]
	ax.errorbar(wl[ind_silicate_con1], flux[ind_silicate_con1], yerr=eflux[ind_silicate_con1], 
	            fmt='o', c='gray', markersize=3.0, linewidth=1.0, capsize=3, capthick=1)
	ax.errorbar(wl[ind_silicate_con2], flux[ind_silicate_con2], yerr=eflux[ind_silicate_con2], 
	            fmt='o', c='gray', markersize=3.0, linewidth=1.0, capsize=3, capthick=1)
	# interpolated continuum at the absorption
	ax.errorbar(silicate_min, cont_silicate, yerr=econt_silicate, fmt='o', c='black', 
	            markersize=3.0, linewidth=1.0, capsize=3, capthick=1)
	
	# fluxes in the absorption
	ind_silicate_min = np.where((wl>=(silicate_min-silicate_window/2)) & (wl<=(silicate_min+silicate_window/2)))[0]
	ax.errorbar(wl[ind_silicate_min], flux[ind_silicate_min], yerr=eflux[ind_silicate_min],  
	            fmt='o', c='gray', markersize=3.0, linewidth=1.0, capsize=3, capthick=1)
	# mean flux of the absorption 
	flux_silicate_min = np.mean(flux[ind_silicate_min])
	eflux_silicate_min = np.std(flux[ind_silicate_min]) / np.sqrt(flux[ind_silicate_min].size)
	ax.errorbar([silicate_min], [flux_silicate_min], yerr=eflux_silicate_min, 
	            fmt='o', c='red', markersize=3.0, linewidth=1.0, capsize=3, capthick=1, zorder=3)
	
	plt.text(0.55, 0.8, 'Silicate index='+str(round(si_index,2))+'$\pm$'+str(round(esi_index,2)), transform=fig.transFigure, size=12)
	
	plt.xlim(xmin, xmax)
	plt.ylim(ymin, ymax)
	ax.xaxis.set_minor_locator(AutoMinorLocator())
	ax.yaxis.set_minor_locator(AutoMinorLocator())
	ax.grid(True, which='both', color='gainsboro', alpha=0.5)
	
	plt.xlabel(r'$\lambda\ (\mu$m)', size=12)
	plt.ylabel(r'$F_\nu\ (Jy)$', size=12)
	if plot_title is not None: plt.title(plot_title, size=12)
	else: plt.title('Silicate Index Measurement', size=12)
	
	if plot_save: plt.savefig('silicate_index_measurement.pdf', bbox_inches='tight')
	plt.close()

	return
